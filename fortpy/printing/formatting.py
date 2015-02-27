"""Module with functions for formatting the output of code or other data
in pretty-print for human-reading.
"""

def present_params(paramlist, spacing = 0, maxchars=90, linecont=", &"):
    """Creates the (paramlist) for a method call formatted nicely for calls
    with lots of parameters."""
    #The +2 is spacing is for the tab indent at the start of the line.
    #The +3 is for indent and the extra parenthesis at the start of the call.
    line = []
    length = 0
    result = []

    for param in paramlist:
        extra = len(list(param))
        if length + extra + 2 + spacing > maxchars:
            result.append(", ".join(line) + linecont)
            line = [ param ]
            length = extra + 2
        else:
            line.append(param)
            length += extra + 2

    #Add on the remaining bits of the line
    result.append(", ".join(line))

    return "\n{}".format(" ".join([ "" for i in range(spacing + 3)])).join(result)
