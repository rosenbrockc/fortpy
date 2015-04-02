import os
import sys
from . import settings
import hashlib
from stat import S_ISDIR
import errno

def coderelpath(coderoot, relpath):
    """Returns the absolute path of the 'relpath' relative to the specified code directory."""
    from os import chdir, getcwd, path
    cd = getcwd()
    chdir(coderoot)
    result = path.abspath(relpath)
    chdir(cd)
    return result

paramiko = None

class FileSupport(object):
    """This tramp module provides support for having intellisense support
    over SSH when emacs uses tramp for editing the files on the super computer.
    """
    def __init__(self):
        self.config = sys.modules["config"]
        self.ssh = None
        self.ftp = None
        self._pkey = None

    def _setup_ssh(self):
        """Initializes the connection to the server via SSH."""
        global paramiko
        if paramiko is none:
            import paramiko

        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh.connect(self.server, username=self.user, pkey=self.pkey)

    def _setup_ftp(self):
        """Initializes an SFTP session for copying files from the 
        remote server."""
        self._check_ssh()
        self.ftp = self.ssh.open_sftp()
        
    def _check_ssh(self):
        """Makes sure an ssh connection is available. If it isn't start one."""
        if self.ssh is None:
            self._setup_ssh()

    def _check_ftp(self):
        """Makes sure an ssh connection is available. If it isn't start one."""
        if self.ftp is None:
            self._setup_ftp()

    def is_ssh(self, filepath):
        """Returns true if the specified filepath is a tramp path for
        accessing a file over SSH."""
        return filepath[:4] == "/ssh"
    
    def abspath(self, path):
        """Returns the absolute path to the specified relative or user-relative
        path. For ssh paths, just return the full ssh path."""
        if self.is_ssh(path):
            return path
        else:
            return os.path.abspath(path)

    def dirname(self, path):
        """Returns the full path to the parent directory of the specified
        file path."""
        if self.is_ssh(path):
            remotepath = self._get_remote(path)
            remotedir = os.path.dirname(remotepath)
            return self._get_tramp_path(remotedir)
        else:
            return os.path.dirname(path)
        
    def touch(self, filepath):
        """Touches the specified file so that its modified time changes."""
        if self.is_ssh(filepath):
            self._check_ssh()
            remotepath = self._get_remote(filepath)
            stdin, stdout, stderr = self.ssh.exec_command("touch {}".format(remotepath))
            stdin.close()
        else:
            os.system("touch {}".format(filepath))
        
    def expanduser(self, filepath, ssh=False):
        """Replaces the user root ~ with the full path on the file system.
        Works for local disks and remote servers. For remote servers, set
        ssh=True."""
        if ssh:
            self._check_ssh()
            stdin, stdout, stderr = self.ssh.exec_command("cd; pwd")
            stdin.close()
            remotepath = filepath.replace("~", stdout.read().split()[0])
            return self._get_tramp_path(remotepath)
        else:
            return os.path.expanduser(filepath)

    def exists(self, filepath):
        """Determines if the specified file/folder exists, even if it
        is on a remote server."""
        if self.is_ssh(filepath):
            self._check_ftp()
            remotepath = self._get_remote(filepath)
            try:
                self.ftp.stat(remotepath)
            except IOError as e:
                if e.errno == errno.ENOENT:
                    return False
            else:
                return True
        else:
            return os.path.exists(filepath)
    
    def getmtime(self, filepath):
        """Gets the last time that the file was modified."""
        if self.is_ssh(filepath):
            self._check_ftp()
            source = self._get_remote(filepath)
            mtime = self.ftp.stat(source).st_mtime
        else:
            mtime = os.path.getmtime(filepath)

        return mtime                        

    def read(self, filepath):
        """Returns the entire file as a string even if it is on a remote
        server."""
        target = self._read_check(filepath)

        if os.path.isfile(target):
            with open(target) as f:
                string = f.read()

            #If we got this file via SSH, delete it from the temp folder
            if self.is_ssh(filepath):
                os.remove(target)
        else:
            string = ""

        return string

    def readlines(self, filepath):
        """Returns the list of strings extracted from the file using either
        open() if it is a local file, or by copying over the network using
        sftp over SSH."""
        target = self._read_check(filepath)

        if os.path.isfile(target):
            with open(target) as f:
                lines = f.readlines()

            #If we got this file via SSH, delete it from the temp folder
            if self.is_ssh(filepath):
                os.remove(target)
        else:
            lines = []

        return lines

    def _read_check(self, filepath):
        """Returns the path of a file on the *local* system that can be read
        from. If the filepath is on a remote server, the file is first copied
        locally."""
        if self.is_ssh(filepath):
            self._check_ftp()
            #First we need to generate a file path on the local system to
            #copy the file to.
            source = self._get_remote(filepath)
            target = self._get_hashed_path(filepath)
            self.ftp.get(source, target)
            #Now we can just read it with the normal python commands.
        else:
            target = filepath
        
        return target

    def walk(self, dirpath):
        """Performs an os.walk on a local or SSH filepath."""
        if self.is_ssh(dirpath):
            self._check_ftp()
            remotepath = self._get_remote(dirpath)
            return self._sftp_walk(remotepath)
        else:
            return os.walk(dirpath)

    def _sftp_walk(self, remotepath):
        """Performs the same function as os.walk but over the SSH channel."""
        #Get all the files and folders in the current directory over SFTP.
        path = remotepath # We need this instance of path for the yield to work.
        files=[]
        folders=[]
        for f in self.ftp.listdir_attr(remotepath):
            if S_ISDIR(f.st_mode):
                folders.append(self._get_tramp_path(f.filename))
            else:
                files.append(f.filename)
        #We use yield so that if there are really large folders with 
        #complicated structures, we don't wait forever while the SFTP
        #keeps traversing down the tree.
        yield path, folders, files

        #Now call this method recursively for each of sub-directories.
        for folder in folders:
            new_path = os.path.join(remotepath, folder)
            for x in self._sftp_walk(new_path):
                yield x

    def _get_tramp_path(self, remotepath):
        """Creates a tramp path of the specified absolute path on the remote
        server."""
        return "/ssh:{}@{}:{}".format(self.user, self.server, remotepath)

    def _get_remote(self, filepath):
        """Gets the full path to the file on the remote server by parsing an
        emacs tramp style buffer-file-name."""
        return filepath.split(":")[-1]

    def _get_hashed_path(self, path):
        """Returns an md5 hash for the specified file path."""
        return self._get_path('%s.pkl' % hashlib.md5(path.encode("utf-8")).hexdigest())

    def _get_path(self, filename):
        """Creates the cache directory if it doesn't already exist. Returns the
        full path to the specified file inside the cache directory."""
        tempdir = settings._temp_directory
        if not os.path.exists(tempdir):
            os.makedirs(tempdir)
        return os.path.join(tempdir, filename)

    def _get_pkey(self):
        """Gets an RSAKey object for the private key file so that we can
        copy files without logging in with user/password."""
        keypath = self.config.server["pkey"]
        with open(os.path.expanduser(keypath)) as f:
            pkey = paramiko.RSAKey.from_private_key(f)

        return pkey

    @property
    def server(self):
        """Returns the name of the SSH server to connect to."""
        return self.config.server["name"]

    @property
    def user(self):
        """Returns the username for logging in to the SSH server."""
        return self.config.server["username"]

    @property
    def pkey(self):
        """Returns the private key for quick authentication on the SSH server."""
        if self._pkey is None:
            self._pkey = self._get_pkey()

        return self._pkey
