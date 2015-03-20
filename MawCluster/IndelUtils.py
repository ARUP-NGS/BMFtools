# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True


import pysam

"""
Contains utilities for working with indels for HTS data.
"""

def FilterByIndelRelevance(inBAM, indelOutputBAM="default", otherOutputBAM="default"):
    """
    Writes reads potentially relevant to an indel to indelOutputBAM and
    other reads to the otherOutputBAM.
    idRel stands for indel relevant.
    idIrl stands for indel irrelevant.
    """
    if(indelOutputBAM == "default"):
        indelOutputBAM = ".".join([inBAM.split(".") + ["idRel", "bam"]])
    if(otherOutputBAM == "default"):
        otherOutputBAM = ".".join([inBAM.split(".") + ["idIrl", "bam"]])
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    indelHandle = pysam.AlignmentFile(indelOutputBAM, "wb")
    otherHandle = pysam.AlignmentFile(otherOutputBAM, "wb")
    for entry in inHandle:
        if entry.is_read1:
            read1 = entry
            continue
        else:
            read2 = entry
        assert read1.query_name == read2.query_name
        if(IsIndelRelevant(read1) or IsIndelRelevant(read2)):
            indelHandle.write(read1)
            indelHandle.write(read2)
        else:
            otherHandle.write(read1)
            otherHandle.write(read2)
    inHandle.close()
    otherHandle.close()
    indelHandle.close()
    return indelOutputBAM, otherOutputBAM


def IsIndelRelevant(read):
    """
    True if considered relevant to indels.
    False otherwise.
    """
    if(read.cigarstring is None):
        # This read is simply unmapped. Let's give it a chance!
        return True
    if("I" in read.cigarstring or "D" in read.cigarstring
       or "S" in read.cigarstring):
        return True
    try:
        if(read.opt("SV") != "NF"):
            return True
    except KeyError:
        # No SV tag was found, not much we can do.
        pass
    return False