# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True


import pysam
import cython

"""
Contains utilities for working with indels for HTS data.
"""

def FilterByIndelRelevance(inBAM, indelOutputBAM="default", otherOutputBAM="default",
                           minFamSize=3):
    """
    Writes reads potentially relevant to an indel to indelOutputBAM and
    other reads to the otherOutputBAM.
    idRel stands for indel relevant.
    idIrl stands for indel irrelevant.
    """
    if(indelOutputBAM == "default"):
        indelOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indRel", "bam"])
    if(otherOutputBAM == "default"):
        otherOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indIrl", "bam"])
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    indelHandle = pysam.AlignmentFile(indelOutputBAM, "wb", template=inHandle)
    otherHandle = pysam.AlignmentFile(otherOutputBAM, "wb", template=inHandle)
    for entry in inHandle:
        if entry.is_read1:
            read1 = entry
            continue
        else:
            read2 = entry
        assert read1.query_name == read2.query_name
        if(IsIndelRelevant(read1, minFam=minFamSize) or
           IsIndelRelevant(read2, minFam=minFamSize)):
            indelHandle.write(read1)
            indelHandle.write(read2)
        else:
            otherHandle.write(read1)
            otherHandle.write(read2)
    inHandle.close()
    otherHandle.close()
    indelHandle.close()
    return indelOutputBAM, otherOutputBAM


@cython.locals(keepSoft=cython.bint, keepUnmapped=cython.bint,
               minFam=cython.long)
def IsIndelRelevant(read, minFam=1, keepSoft=False,
                    keepUnmapped=False):
    """
    True if considered relevant to indels.
    False otherwise.
    """
    if(minFam != 1):
        if(read.opt("FM") < minFam):
            return False
    if(read.cigarstring is None):
        # This read is simply unmapped. Let's give it a chance!
        if(keepUnmapped):
            return True
        return False
    if("I" in read.cigarstring or "D" in read.cigarstring):
        return True
    if(keepSoft is True and "S" in read.cigarstring):
        return True
    try:
        if(read.opt("SV") != "NF"):
            return True
    except KeyError:
        # No SV tag was found, not much we can do.
        pass
    return False
