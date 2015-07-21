
"""
Contains custom exceptions and other tools for handling errors,
particularly with cython's unique error-handling tools.
"""


class PermissionException(Exception):
    """
    Thrown when permissions are likely the preventing culprit
     for a given operation.
    """
    pass


class ImproperArgumentError(ValueError):
    """
    Thrown when an argument is the correct type but the object
    is not a suitable target for it. e.g., laying out an unmapped
    read.
    """
    def __init__(self, message=""):
        super(ValueError, self).__init__()
        self.message = message
        if(message != ""):
            print("Improper Argument Error: %s" % message)


class UnsetRequiredParameter(ValueError):
    def __init__(self, message=""):
        super(ValueError, self).__init__()
        self.message = message
        if(message != ""):
            print(message)


class ConsideredHarmful(ValueError):
    def __init__(self, message=""):
        super(ValueError, self).__init__()
        self.message = message
        if(message != ""):
            print(message)


class IllegalArgumentError(ValueError):
    def __init__(self, message=""):
        super(ValueError, self).__init__()
        self.message = message
        if(message != ""):
            print("Illegal Argument Error: %s" % message)


class AbortMission(Exception):
    """
    Thrown when one wants to prematurely escape due to something unexpected.
    """
    def __init__(self, message=""):
        super(Exception, self).__init__()
        self.message = message


class MissingExternalTool(Exception):
    """
    Thrown when one wants to prematurely escape due to something unexpected.
    """
    def __init__(self, message=""):
        try:
            print("Repr of Chapman object: %s" % repr(globals()['Chapman']))
        except KeyError:
            print("Note: Chapman object not in globals. Globals:"
                  " %s" % repr(globals()))
        super(Exception, self).__init__()
        self.message = message


class ThisIsMadness(Exception):
    """
    Thrown when something just doesn't seem right and existing exceptions
    don't have enough panache.
    """
    def __init__(self, message=""):
        print(FPStr)
        # Call the base class constructor with the parameters it needs
        super(Exception, self).__init__(message)
        self.message = message


class ThisIsHKMadness(Exception):
    """
    Thrown when something just doesn't seem right and existing exceptions
    don't have enough panache.
    """
    def __init__(self, message=""):
        print(HKStr)
        # Call the base class constructor with the parameters it needs
        super(Exception, self).__init__(message)
        self.message = message


class FunctionCallException(Exception):
    """
    Thrown when errors arise when calling functions outside of the
    scope of CalledProcessError.
    """
    def __init__(self, call, message, shell):
        super(Exception, self).__init__(message)
        if(shell is True):
            call = call.split("#")[0]
        else:
            call = call.split("#")[1]
        print("Failed call: %s" % call)
        print("Message: %s" % message)
        print("Was shell: %s" % shell)


class MissingGlobalVariable(Exception):
    """
    Thrown when a global variable one expects doesn't exist.
    """
    def __init__(self, variable, message="Missing global variable ???"):
        if(message == "Missing global variable ???"):
            message = "Missing global variable %s ???" % variable
        self.variable = variable
        super(Exception, self).__init__(message)


FPStr = ("............................................________ \n............"
         "........................,.-'\"...................``~. \n..........."
         "..................,.-\"...................................\"-., \n."
         "........................,/........................................."
         "......\":, \n.....................,?..............................."
         "......................., \n.................../...................."
         ".......................................,} \n................./....."
         ".................................................,:`^`..} \n......."
         "......../...................................................,:\"..."
         "....../ \n..............?.....__..................................."
         "......:`.........../ \n............./__.(.....\"~-,_..............."
         "...............,:`........../ \n.........../(_....\"~,_........\"~,"
         "_....................,:`........_/ \n..........{.._$;_......\"=,_.."
         ".....\"-,_.......,.-~-,},.~\";/....} \n...........((.....*~_......."
         "\"=-._......\";,,./`..../\"............../ \n...,,,___.`~,......\"~"
         ".,....................`.....}............../ \n............(....`=-"
         ",,.......`........................(......;_,,-\" \n............/.`~"
         ",......`-...................................../ \n.............`~.*"
         "-,.....................................|,./.....,__ \n,,_.........."
         "}.>-._...................................|..............`=~-, \n..."
         "..`=~-,__......`,................................. \n.............."
         ".....`=~-,,.,............................... \n...................."
         "............`:,,...........................`..............__ \n...."
         ".................................`=-,...................,%`>--==`` "
         "\n........................................_..........._,-%.......` "
         "\n...................................,")


HKStr = ("\n"
         "   ( )---( )   \n"
         "   /       \   \n"
         "  +  0  o 0 +   ~~ This is madness!  \n"
         "   + _____ +   \n"
         "     / ^ -L   \n"
         "   (~( \/ | )   \n")
