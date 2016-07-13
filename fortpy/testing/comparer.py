from .. import msg
import xml.etree.ElementTree as ET
import re
import os
from .results import CompareResult, BodyResult, BlockResult, CompatibilityResult, ListResult
from .templates import FileTemplate, LineGroup

class FileRepresentation(object):
    """The python-valued representation of a templated file.

    :arg lines: the original lines of the file as strings.
    :arg template: a template.FileTemplate used to intepret the lines.
    :arg version: the version of the file lines extracted from the first line.
    :arg path: the path to the file in the file system.
    
    :attr clean: the cleaned (decommented) lines.
    :attr simple: True if there is no template to intepret the file.
    :attr preamble: a dictionary of values extracted from the preamble. Keys
      are the line identifiers, values is a LIST of LineValue objects.
    :attr stored: the dictionary of 'global' variables used by the body and
      other lines in the preamble.
    :attr body: a list of body blocks extracted from the file. They are ordered
      by appearance in the file.
    :attr extracted: specified whether the template could find all the data it
      needed during extraction.
    """
    def __init__(self, lines, template, version, path):
        self.lines = lines
        self.template = template
        self.version = version
        self.path = path

        self.simple = False
        self.clean = None

        self.preamble = {}
        self.stored = {}
        self.body = []

        self._clean_values()
        self.extracted = self._extract()

    def _extract(self):
        """Extracts the data in python types line-by-line."""    
        success = True

        if self.template is not None:
            #We need to use the file lines to parse each line, ignoring any comments.
            l, failed = self._extract_preamble()
                
            #The rest of the output file is the body. The template decides how to
            #process the body.
            if not failed:
                success = self._extract_body(l)
            else:
                success = False
        else:
            #We are just comparing text values from each file, return the cleaned list.
            self.simple = True
            self.values = self.clean

        return success

    def _extract_preamble(self):
        """Extracts data using the template for the preamble."""
        #The preamble allows a template to specify some variables that are used
        #to extract data in the body of the output file. These are put in self.stored.
        l = 0
        failed = False

        for line in self.template.preamble:
            count = self._get_line_count(line, self.stored)
            if isinstance(line, LineGroup):
                lcounts = [self._get_line_count(gline, self.stored) for gline in line.lines]
                line.update_counts(lcounts, count)
                count = count*sum(lcounts)
            
            lineres = []
            for i in range(count):
                #If the output file doesn't have enough lines, we can't compare the files.
                if l + i < len(self.clean):
                    if isinstance(line, LineGroup):
                        values = line.parse(self.clean[l + i], i)
                    else:
                        values = line.parse(self.clean[l + i])
                    lineres.append(values)
                    #If the result stored any values, we need to make those publicly
                    #available to future lines.
                    self.stored.update(values.stored)
                else:
                    failed = True
                    break
                
            if failed:
                break                    

            #Associate all the lines with the identifier from the template
            self.preamble[line.identifier] = lineres
            l += count
            
        return (l, failed)

    def _extract_body(self, l):
        """Extracts data from the body using the preamble."""
        #Any values stored inside a body declaration only
        #apply to the current iteration of the body template
        count = self._get_line_count(self.template, self.stored)
        success = True
        if count == 0:
            return True

        if count > 1:
            #We will just process the body template count times
            for i in range(count):
                run = self._read_body_template(l, self.template.body)
                l = run[0]
                self.body.append(run[1])                    
                if not run[2]:
                    #It failed, stop the loop and exit
                    success = False
                    break
        elif self.template.stop == "EOF":
            #We will process the body until we run out of items
            while l < len(self.clean):
                run = self._read_body_template(l, self.template.body)
                l = run[0]
                self.body.append(run[1])                  
                if not run[2]:
                    #It failed, stop the loop and exit
                    success = False
                    break
        else:
            #The body is only being read through once.
            run = self._read_body_template(l, self.template.body)
            self.body.append(run[1])                  
            if not run[2]:
                #It failed, stop the loop and exit
                success = False

        return success

    def _read_body_template(self, l, body):
        """Reads the body template in once from the current (l) position in the lines."""
        #Each time we read the body, we get a dictionary of identifiers and values.
        result = {}
        success = True

        #The body stores results at this local level, they aren't available outside of
        #an iteration.
        stored = {}
        for line in body:
            count = self._get_line_count(line, stored)
            if isinstance(line, LineGroup):
                lcounts = [self._get_line_count(gline, stored) for gline in line.lines]
                line.update_counts(lcounts, count)
                count = count*sum(lcounts)
                
            if count == 0:
                #We don't want to read anything in, just keep moving along.
                continue

            lineres = []
            for i in range(count):
                #If the output file doesn't have enough lines, we can't compare the files.
                if l + i < len(self.clean):
                    if isinstance(line, LineGroup):
                        values = line.parse(self.clean[l + i], i)
                    else:
                        values = line.parse(self.clean[l + i])
                    lineres.append(values)
                    #If the result stored any values, we need to make those publicly
                    #available to future lines.
                    stored.update(values.stored)
                else:
                    success = False
                    break

            if not success:
                break

            #Associate all the lines with the identifier from the template
            result[line.identifier] = lineres
            l += count

        return (l, result, success)

    def _get_line_count(self, line, store):
        """Determines the number of lines that the specified template should be used for."""
        if isinstance(line.count, int):
            count = line.count
        else:
            #The count refers to a variable that should be stored from an
            #earlier line.
            if line.count in store:
                count = store[line.count]
            elif "$" in line.count:
                #We just use eval() on the attribute value they specified. Any variable
                #names need to be replaced by the appropriate value in the store dict.
                import operator
                rx = re.compile(r"(?<=[(\s,\-+*/%])(?<![\[.])[A-Za-z]+[\w\d_]*")
                nl = line.count.replace("$", "")
                for varname in rx.findall(line.count):
                    if varname not in ["if", "then", "else", "operator"]:
                        nl = nl.replace(varname, 'store["{}"]'.format(varname))
                count = eval(nl)
            else:
                msg.err("reference to stored value that does not exist '{}'".format(line.count))
                msg.gen("STORE: {}".format(list(store.keys())))
                exit(1)

        return count

    def _clean_values(self):
        """Extracts the values from the specified lines using a template."""
        #First extract any comments from the lines
        self.clean = []
        comments = self._get_comment_char()
        for line in self.lines:
            decommented = line.split(comments)[0].strip()
            if decommented != "":
                self.clean.append(decommented)

    def _get_comment_char(self):
        """Gets the character that was used to designate a comment in the output file."""
        default = "#"
        if self.template is not None:
            if self.template.comments != "":
                return self.template.comments
            else:
                return default
        else:
            return default

def compare_representations(rep1, rep2, mode = "default"):
    """Determines how similar the specified representations are.

    :arg mode: the level/strictness of comparisons to make.
    """
    #We will compare the preamble and body separately and see how closely
    #related they are. We have to take compatibility specs into account
    #and what can actually be compared if the versions are different.
    if not rep1.simple or not rep2.simple:
        #Make sure the mode specifed for comparisons exists
        _compare_validate(rep1, mode)
        _compare_validate(rep2, mode)

        preamble = _compare_preamble(rep1, rep2, mode)
        body = _compare_body(rep1, rep2, mode)
        result = CompareResult(preamble, body, rep1.path, rep2.path)
        result.mode = mode
        return result
    elif rep2.simple:
        #This is a really simple file comparison of straight values.
        #Since all the other comparers need template information, we
        #will just do a line by line comparison here.               
        result = ListResult(rep1.clean, rep2.clean)
        elcount = min([len(rep1.clean), len(rep2.clean)])
        for i in range(elcount):
            if rep1.clean[i] != rep2.clean[i]:
                result.different.append((rep1.clean[i], rep2.clean[i]))
            else:
                result.common += 1

        return result
    else:
        msg.err("a templated representation cannot be compared to a simple one.")
        
def _compare_validate(rep, mode):
    """Validates the contents of the representation."""
    if not mode in rep.template.comparisons:
        msg.err("the mode '{}' specified for comparisons does not exist.".format(mode))
        exit(1)
    if not mode in rep.template.outcomes:
        msg.err("the mode '{}' specified for outcomes does not exist.".format(mode))
        exit(1)

def _compare_preamble(rep1, rep2, mode):
    """Compares the contents of the preamble to the rep. preamble."""
    #The preamble is really just like a single body block once we get to
    #the comparison stage.
    return _compare_block(rep1, rep2, rep1.preamble, rep2.preamble, mode, None, None, False)

def _compare_body(rep1, rep2, mode):
    """Compares the contents of the specified body list to this rep. body."""
    #if a key was specified for the body, we generate a hashtable so that
    #each line will only be matched to a corresponding one that has
    #the same value for the key.
    result = BodyResult(rep1, rep2)
    if len(rep1.body) == 0 and len(rep2.body) == 0:
        #Don't bother computing anything, we have two empty sets.
        return result

    if rep1.template.key is not None:
        #Loop through the list of body elements and re-store them in a dictionary
        #by their key value
        thisd = _get_keyed_dict(rep1.template, rep1.body)
        thatd = _get_keyed_dict(rep2.template, rep2.body)

        for bkey in thisd:
            if bkey in thatd:
                result.blocks[bkey] = _compare_block(rep1, rep2, thisd[bkey], thatd[bkey], mode, bkey)
            else:
                result.only1[bkey] = thisd[bkey]

        #Do the same thing from the other body's perspective
        for bkey in thatd:
            if bkey not in thisd:
                result.only2[bkey] = thatd[bkey]
    else:
        #Loop through by index and compare each indexed body to the other.
        for i in range(len(rep1.body)):
            thisb = rep1.body[i]
            if i < len(rep2.body):
                thatb = rep2.body[i]
                result.blocks[i] = _compare_block(rep1, rep2, thisb, thatb, mode, None, i)
            else:
                result.only1[i] = rep1.body[i]

        #See if we had more entries in body2 than body1
        if len(rep2.body) > len(rep1.body):
            for i in range(len(rep1.body), len(rep2.body)):
                result.only2[i] = rep2.body[i]

    return result                


def _compare_block(rep1, rep2, thisb, thatb, mode, bkey = None, index = None, isbody = True):
    """Compares two blocks that are supposed to be similar.

    :arg rep1, rep2: FileRepresentations that the body blocks came from.
    :arg thisb, thatb: dictionaries that represent a single pass through
      the body template definition. Keys are the line identifiers, values
      are lists of LineValue objects, one for each read of a line, or
      for 'lines' tags, one for each pass with the same line template.
    :arg key: if the body blocks are keyed by value, the value that both
      blocks are being compared on.
    :arg index: if the body blocks are compared sequentially, the current
      index of the blocks in both parent lists.
    :arg isbody: specifies whether the block is being compared as part of
      the file body (True), or whether it is just the preamble 'block'.
    """
    result = BlockResult(rep1.template.outcomes[mode], bkey, index, rep1.template)
    #First look at those identifiers that are either common to both or
    #have compatibility specified in the first block's templates.
    for key in thisb:            
        #If the values for a line were named, we need to compare by name
        #otherwise, we just compare by index. Thisv and thatv are LISTS
        #of line value objects.
        thisv = thisb[key]
        thatv = None
        if key in thatb:
            thatv = thatb[key]
        else:
            #Check for compatibility setting to link the two versions
            if rep2.version in _get_compat(rep1.template, key, isbody):
                compat = _get_compat(rep1.template, key, isbody)[rep2.version]
                result.results[key] = _compare_compat(compat, thatb, thisv, key, rep1.template.outcomes[mode])

        if thatv is not None:
            if len(thisv) != len(thatv):
                loop = min([len(thisv), len(thatv)])
                msg.warn("line values extracted for key" + 
                         " '{}' in the two files have different".format(key) + 
                         " numbers of elements: {} vs. {}.".format(len(thisv), len(thatv)), 2)
            else:
                loop = len(thisv)

            for j in range(loop):
                if len(list(thisv[j].named.keys())) != 0:
                    res = rep1.template.comparisons[mode].compare_d(thisv[j].named,
                                                                    thatv[j].named, key,
                                                                    rep1.template.outcomes[mode])
                else:
                    res = rep1.template.comparisons[mode].compare_l(thisv[j].values,
                                                                    thatv[j].values, key,
                                                                    rep1.template.outcomes[mode])

                #We have one such result for each element in the list under this
                #line identifier (key). We will make a composite key that shows
                #the list items as being related to the parent line, but unique.
                if loop > 1:
                    result.results["{}.{}".format(key, j)] = res
                else:
                    result.results[key] = res
               
    #Handle compatibility settings that are in the second template pointing
    #to the first one.
    for key in thatb:
        thatv = thatb[key]
        if rep1.version in _get_compat(rep2.template, key, isbody):
            compat = _get_compat(rep2.template, key, isbody)[rep1.version]
            result.results[key] = _compare_compat(compat, thatv, thisb, key, rep2.template.outcomes[mode])

    return result

def _get_compat(template, key, isbody):
    """Gets the compatibility dictionary from the specified template by key
    and whether it is the body compatibility section."""
    if isbody:
        return template._body[key].compatibility
    else:
        return template._preamble[key].compatibility

def _compare_compat(compat, compatv, thatb, label, outcomes):
    """Does a compatibility comparison between the block the specified the
    compatibility mapping and some other block for a certain identifier
    in the compatibility block.

    :arg compat: Compat is a dictionary with target id as key and
      value as a dictionary of mappings from one name to another.
    :arg compatv: the block of line identifiers and LineValue objects that
      come from the same parent as the compatibility dictionary compat.
    :arg thatb: the other block whose values will be compared to those in
      the compatibility block.
    :arg label: the line identifier for this compatibility comparison.
    :arg outcomes: a TemplateOutcomes from the parent of the compatibility
      dictionary of mappings.
    """
    #We just need to go through all the
    #target mappings and do the comparisons
    result = CompatibilityResult(label, outcomes)
    for tkey in compat:
        if tkey in thatb:
            thatv = thatb[target]
            mappings = compat[tkey]
            #We only deal with named values in compatibility mode.
            for namekey in mappings:
                result.total += 1
                if namekey in thatv.named and namekey in compatv.named:
                    if compatv.named[namekey] == thatv.named[mappings[namekey]]:
                        result.common += 1
                    else:
                        result.different += 1
                else:
                    result.key_errors.append("{}.{}".format(tkey, namekey))
        else:
            result.key_errors.append(tkey)
                        
    return result
  
def _get_keyed_dict(template, body):
    """Returns a dictionary of body block dictionaries set by key value."""
    result = {}
    for block in body:
        if isinstance(template.key, list):
            keyval = _get_key_value_list(template, block)
        else:
            keyval = _get_key_value_single(template.key, block)
        if keyval is not None:
            result[keyval] = block
    return result

def _get_key_value_list(template, bodyblock):
    """Extracts a multi-valued key hash from the body block."""
    #If we don't find a value to list the key by, we return none.
    result = None
    values = [_get_key_value_single(k, bodyblock) for k in template.key]

    #We are going to string join the values in a bar-separated list and then
    #compute a hash on the string is the value to return.
    svalues = []
    for oneval in values:
        if isinstance(oneval, list):
            strval = ",".join([str(o) for o in oneval])
        else:
            strval = str(oneval)
        svalues.append(strval)

    result = hash("|".join(svalues))
    return result    

def _get_key_value_single(key, bodyblock):
    """Extracts the value of the key from the specified body block."""
    #If we don't find a value to list the key by, we return none.
    result = None
    if "." in key:
        keys = key.split(".")
        if len(keys) == 2 and keys[0] in bodyblock:
            #The lineval will be a list of line values. In the case where a lines
            #tag specifies multiple values for the same identifier, we will only
            #use the first elements named value.
            lineval = bodyblock[keys[0]][0]
            if keys[1] in lineval.named:
                result = lineval.named[keys[1]]
    else:
        if key in bodyblock:
            result = bodyblock[key].values[0]

    if result is None:
        msg.err("block {} did not return a valid key value for '{}'.".format(bodyblock, key))
        exit(1)

    return result    

class FileComparer(object):
    """Class for comparing different versions of output files.

    :arg template_folder: the path to the folder that contains the
    XML templates for comparing files across versions.
    """
    def __init__(self, fortpy_templates=None, template_folder=None):
        self.fortpy_templates = None
        """The pull path to the directory that holds the templates shipping with fortpy.
        """
        from fortpy.utility import set_fortpy_templates
        set_fortpy_templates(self, fortpy_templates)
        self.templates = {}
        """Dict of FileTemplate instances with the XML file name as key. 
        """
        from os import path
        self.folder = path.abspath(path.expanduser(template_folder))
        """The path to the folder that houses the XML template to use in comparing.
        """
        self.template = None
        """This is the main template used for all the versions of the files that
        will be compared.
        """
        
    def compare(self, source, target, template=None, mode="default", reps=False):
        """Compares the two files using an XML template if one exists.

        :arg source: the path to the first file to compare.
        :arg target: the path to the other file to compare.
        :arg template: the name of XML file to use as template for the files.
        :arg mode: the comparison mode to use defined in the template.
        :arg reps: when True, the file representations used in the comparison are 
          also returned from the function.
        """
        #Try to load the template specified by the testspec.
        self._load_from_xml(template)

        svalues = self.get_representation(source, template)
        tvalues = self.get_representation(target, template)

        #Since we re-use the comparer for multiple unit tests, we need to reset
        #the currently active template.
        self.template = None

        #We can't compare representations that don't exist...
        if svalues is not None and tvalues is not None:
            result = compare_representations(svalues, tvalues, mode)
            if not reps:
                return result
            else:
                return (result, svalues, tvalues)            

    def get_representation(self, path, template=None):
        """Creates a file representation for the specified file path.

        :arg path: the full path to file to get a templated representation of.
        :arg template: the name of the template XML file to use. If unspecified,
          the <fortpy> tag in the file should specify it.
        """
        source = os.path.expanduser(path)
        if not os.path.exists(source):
            msg.err("can't create representation for {}.".format(source) + 
                    " File does not exist.")            
            return None

        with open(source) as f:
            slines = f.readlines()

        #The first line in a file to be compared contains version information
        if len(slines) == 0:
            msg.err("The file {} is empty; can't create representation.".format(source))
            return None

        sf = self._get_fortpy(slines[0], source)
        sv = self._get_file_version(sf)

        #By convention, if the two files have different filenames, we only
        #load a template based on the source name and use it for both
        if self.template is None:
            self.template = self._load_template(source, sf)
        if (self.template is None and template is not None and
            template.lower() in self.templates):
            self.template = self.templates[template.lower()]

        stemplate = self._get_file_template(self.template, sv)

        #Get python-valued representations and compare them. If the first line
        #doesn't have a comment, we obviously have data starting the very first
        #line.
        if "#" in slines[0]:
            svalues = FileRepresentation(slines[1::], stemplate, sv, source)
        else:
            svalues = FileRepresentation(slines, stemplate, sv, source)

        if not svalues.extracted:
            msg.err("ouput file does not have the same format as template.\n{}".format(source))
            return None
        else:
            return svalues

    def _get_file_template(self, template, version):
        """Gets the file template for the specified file version if it exists."""
        if template is not None and version in template.contents:
            return template.contents[version]
        else:
            return None

    def _get_fortpy(self, line1, source):
        """Extracts the fortpy tag from the first line of the file.

        :arg source: the path to the file from which the line was extracted.
        """
        if "#" in line1:
            from fortpy.utility import XML_fromstring
            try:
                lxml = XML_fromstring(line1.split("#")[1], source)
                if lxml.tag == "fortpy":
                    return lxml
                else:
                    return None
            except ET.ParseError as err:
                msg.warn(err.msg, 2)
                msg.info("no version information found in the file. Assuming version 1.", 2)
                return None
        else:
            msg.info("no version information found in the file. Assuming version 1.", 2)
            return None

    def _get_file_version(self, fortpyxml):
        """Extracts file version information from the fortpy xml element."""
        if fortpyxml is not None and "version" in fortpyxml.attrib:
            return int(fortpyxml.attrib["version"])
        else:
            return 1
            
    def _load_from_xml(self, filename):
        """Loads the XML template for the specified XML file name if it exists."""
        if filename is None or filename.lower() in self.templates:
            return

        target = os.path.join(self.folder, filename)
        if os.path.isfile(target):
            self.templates[filename.lower()] = FileTemplate(target)
        else:
            #We could try the default templates directory for fortpy.
            target = os.path.join(self.fortpy_templates, filename)
            if os.path.isfile(target):
                self.templates[filename.lower()] = FileTemplate(target)            

    def _load_template(self, filepath, fortpyxml):
        """Tries to load the XML template for the file at the specified path."""
        #A template file name can be specified in the fortpy tag at the top of
        #the output file
        if fortpyxml is not None and "template" in fortpyxml.attrib:
            template = fortpyxml.attrib["template"]
            from fortpy.utility import get_dir_relpath
            xmlpath = get_dir_relpath(self.folder, template)
            if xmlpath[-4:] != ".xml":
                xmlpath += ".xml"
            xmlname = os.path.split(xmlpath)[1]
            if not os.path.isfile(xmlpath):
                #This could be one of fortpy's built-in templates
                xmlpath = os.path.join(self.fortpy_templates, xmlname)
        elif filepath != "" and filepath is not None:
            #The xml template file will have the same name as the files do
            #but an EXTRA extension of .xml - J.1.out -> J.1.out.xml
            #Look if we have one in the templates folder.
            xmlname = filepath.split("/")[-1] + ".xml"
            if self.folder != "":
                if xmlname not in self.templates:
                    xmlpath = os.path.join(self.folder, xmlname)
        else:
            #We can't get a template with the information we have, stop.
            return None

        if xmlname in self.templates:
            return self.templates[xmlname]
        if os.path.exists(xmlpath):
            msg.info("Using {} as the XML file template for comparison.".format(xmlpath), 2)
            self.templates[xmlname] = FileTemplate(xmlpath)    
            return self.templates[xmlname]
        else:
            msg.warn("Could not find the XML file template for {}.".format(filepath), 2)
            return None
