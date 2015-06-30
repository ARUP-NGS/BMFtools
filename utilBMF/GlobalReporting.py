from utilBMF.ErrorHandling import IllegalArgumentError, MissingGlobalVariable
from copy import copy as ccopy


class LockedDictionary(dict):

    KeysetDict = {}

    def __init__(self):
        if(isinstance(self.__class__.KeysetDict, set)):
            self.KeysetDict = {i: None for i in
                               self.__class__.KeysetDict}
        elif(isinstance(self.__class__.KeysetDict, dict)):
            self.KeysetDict = ccopy(self.__class__.KeysetDict)
        else:
            try:
                # Attempt to coerce to iterable.
                self.KeysetDict = {i: None for i in
                                   list(self.__class__.KeysetDict)}
            except TypeError:
                raise IllegalArgumentError("KeysetDict not iterable...")

    def __setitem__(self, key, value):
        try:
            b = self.KeysetDict[key]  # Check if key is initialized.
            super(LockedDictionary, self).__setitem__(key, value)
        except KeyError:
            raise KeyError("A LockedDictionary can only hold keys with which "
                           "it was initialized. Keyset: "
                           "%s" % ",".join(self.keys()))

    @classmethod
    def __fromdefaultdict__(self, cls):
        return cls.__init__(self, cls.KeysetDict)


class SampleMetrics(LockedDictionary):
    KeysetDict = {"Name": None, "TotalReadCount": None,
                  "MergedReadCount": None, "MeanFamSize": None}


class ReviewDirComponents(LockedDictionary):
    KeysetDict = {"bam": None, "vcf": None}


def getReviewDir():
    if("review_dir" in globals()):
        return globals()['review_dir']
    else:
        raise MissingGlobalVariable(
            'review_dir',
            "Could not load review dir path! "
            "Globals: %s" % repr(globals()))
