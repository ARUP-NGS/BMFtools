cimport cython
cimport pysam.calignmentfile
cimport pysam.calignedsegment
cimport pysam.cfaidx
cimport numpy as np
cimport utilBMF.HTSUtils
from cython cimport bint
from cpython cimport array
from libc.stdint cimport *
from pysam.chtslib cimport *
from numpy cimport ndarray
from utilBMF.HTSUtils cimport PysamToChrDict
from utilBMF.Inliners cimport Num2Nuc
from utilBMF.PysamUtils cimport PysamToChrInline
from utilBMF.cstring cimport struct_str
from pysam.calignedsegment cimport pysam_get_l_qname, bam1_t
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t
ctypedef cython.str cystr
ctypedef TagBamPipeHG37 TagBamPipeHG37_t
ctypedef TagBamPipe TagBamPipe_t
ctypedef BamPipe BamPipe_t
ctypedef struct_str struct_str_t
ctypedef array.array py_array
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t

import cython


cdef cystr cBarcodeTagCOBam(pysam.calignmentfile.AlignmentFile bam,
                            pysam.calignmentfile.AlignmentFile outbam)
cpdef AlignedSegment_t pTagAlignedSegmentHG37(AlignedSegment_t read)
cpdef cystr pBarcodeTagCOBam(cystr bam, cystr outbam=?)

cdef dict cGetCOTagDict(AlignedSegment_t read)

cpdef dict pGetCOTagDict(AlignedSegment_t read)

cdef inline char * get_Z_tag(AlignedSegment_t read, cystr btag):
    cdef uint8_t * data
    data = bam_aux_get(read._delegate, btag)
    return <char*>bam_aux2Z(data)


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline int8_t BarcodeHD(bam1_t * query, bam1_t * cmp,
                             int8_t length) nogil:
    cdef char* str1
    cdef char* str2
    cdef int8_t ret = 0
    cdef uint8_t index
    str1 = <char*>bam_aux2Z(bam_aux_get(query, "BS"))
    str2 = <char*>bam_aux2Z(bam_aux_get(query, "BS"))
    for index from 0 <= index < length:
        if(str1[index] != str2[index]):
            ret += 1
    return ret


cdef double getAF(AlignedSegment_t read)
cdef double getSF(AlignedSegment_t read)
# cpdef AlignedSegment_t TagAlignedSegment(AlignedSegment_t read)
cdef AlignedSegment_t TagAlignedSegment(
        AlignedSegment_t read, dict RefIDDict=?)
cdef AlignedSegment_t TagAlignedSegmentHG37(
        AlignedSegment_t read)

cdef cystr RPStringNonHG37(AlignedSegment_t read, dict RefIDDict=?)


cdef inline cystr RPString(AlignedSegment_t read):
    return (PysamToChrInline(read.reference_id) + ":%s," % read.pos +
            PysamToChrInline(read.rnext) +
            ":%s" % read.mpos)


cdef inline cystr RescueFlag(AlignedSegment_t read):
    return "%s:%s,%s:%s-%s" % (read.pos, read.reference_id,
                               read.mpos, read.rnext,
                               read.is_read1)


cdef inline bint IS_READ1(AlignedSegment_t read):
    return read.is_read1


cdef inline int MPOS(AlignedSegment_t read):
    return read.mpos


cdef inline int POS(AlignedSegment_t read):
    return read.pos


cdef inline int REF_ID(AlignedSegment_t read):
    return read.reference_id


cdef inline int RNEXT(AlignedSegment_t read):
    return read.rnext


cdef class BamPipe:
    """
    Creates a callable function which acts on a BAM stream.

    :param function - callable function which returns an input BAM object.
    :param bin_input - boolean - true if input is BAM
    false for TAM/SAM
    :param bin_output - boolean - true to output in BAM format.
    :param uncompressed_output - boolean - true to output uncompressed
    BAM records.
    """
    cdef public object function
    cdef public bint bin_input, bin_output, uncompressed_output
    cdef public AlignmentFile_t inHandle, outHandle
    cpdef process(self)
    cdef write(self, AlignedSegment_t read)


cdef class TagBamPipeHG37(BamPipe):
    pass


cdef class TagBamPipe:
    cdef public bint bin_input, bin_output, uncompressed_output
    cdef public AlignmentFile_t inHandle, outHandle
    cdef write(self, AlignedSegment_t read)
    cdef public dict RefIDDict
    cpdef process(self)


cdef inline cystr AStostring(AlignedSegment_t read, AlignmentFile_t handle):
    cdef cython.str cigarstring, mate_ref, ref
    if(read.is_unmapped):
        ref = "*"
    else:
        ref = handle.getrname(read.reference_id)
    if(read.mate_is_unmapped):
        mate_ref = "*"
    else:
        mate_ref = handle.getrname(read.rnext) if(
            read.rnext != read.reference_id and read.rnext != -1) else "="
    cigarstring = read.cigarstring if(
        read.cigarstring is not None) else "*"
    return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        read.query_name, read.flag, ref,
        read.pos + 1, read.mapq, cigarstring, mate_ref, read.mpos + 1,
        read.template_length, read.seq, read.qual, read.get_tag_string())
