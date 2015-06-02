
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


class IllegalArgumentError(ValueError):
    pass


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
        raise Exception(message)


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