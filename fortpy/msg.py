from termcolor import cprint

#This module handles writing to the terminal or a log file with support
#for coloring for warnings, errors, etc.
def warn(msg):
    """Prints the specified message as a warning; prepends "WARNING" to
    the message, so that can be left off.
    """
    cprint("WARNING: " + msg, "yellow")

def err(msg):
    """Prints the specified message as an error; prepends "ERROR" to
    the message, so that can be left off.
    """
    cprint("ERROR: " + msg, "red")

def info(msg):
    """Prints the specified message as information."""
    cprint(msg, "cyan")

def okay(msg):
    """Prints the specified message as textual progress update."""
    cprint(msg, "green")

def gen(msg):
    """Prints the message as generic output to terminal."""
    cprint(msg, "blue")

def blank(n=1):
    """Prints a blank line to the terminal."""
    for i in range(n):
        print("")
