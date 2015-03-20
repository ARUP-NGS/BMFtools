# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

"""
Contains tools for accurately identifying somatic mutations based on
different experimental 
"""

from MawCluster.BCVCF import IterativeVCFFile
import numpy as np
# PROTOCOLS is a list of lowered strings.
PROTOCOLS = ["ffpe", "amplicon", "cf", "other"]


class ExpProtocol:
    """
    Contains information for an experimental protocol for the purposes of
    informing analysis.
    exp should be a string of protocols separated by commas.
    For example, if an experiment is an FFPE amplicon assay, a valid exp
    argument would be: "ffpe,amplicon"
    """
    def __init__(self, exp="default"):
        if(exp == "default"):
            raise ThisIsMadness("Experimental protocol must be set!")
        if(sum([i not in PROTOCOLS for i in exp.lower().split(",")]) != 0):
            raise ThisIsMadness("Experiment must one of those in PROTOCOL"
                                "S: {}".format(PROTOCOLS))
        self.exp = exp
    
    def getEx(self):
        return self.exp
    
    def setEx(self, newExp):
        self.exp = newExp
        if(sum([i not in PROTOCOLS for i
                in self.exp.lower().split(",")]) != 0):
            raise ThisIsMadness("Experiment must one of those in PROTOCOL"
                                "S: {}".format(PROTOCOLS))


def GetDeaminationBackgroundFrequencies(inVCF):
    
    """
    Returns a list of raw base frequencies for G->A and C->T.
    Only accepts SNVCrawler's VCF output as it requires those INFO fields.
    """
    IVCFObj = IterativeVCFFile(inVCF)
    GAFreqArray = []
    CTFreqArray = []
    for line in IVCFObj:
        if(line.REF == "C"):
            CTFreqArray.append(float(
                line.InfoDict["MAFS"].split(";")[2].split(">")[1]))
        elif(line.REF == "G"):
            GAFreqArray.append(float(
                line.InfoDict["MAFS"].split(";")[0].split(">")[1]))
    raise ThisIsMadness("This isn't finished yet.")
    return

# somehow I need to work this into the analysis.
# apply the CA/GT filter to FFPE samples, apply the primer filter to amplicon