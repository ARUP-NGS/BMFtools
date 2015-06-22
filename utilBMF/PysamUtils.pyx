from pysam import AlignedSegment
import pysam


cdef AlignedSegment_t CopyAlignedSegment(AlignedSegment_t template):
    cdef py_array tmpQual = template.query_qualities
    retValue = AlignedSegment()
    retValue.query_sequence = template.query_sequence
    retValue.query_qualities = tmpQual
    retValue.query_name = template.query_name + "_copy"
    retValue.tlen = template.tlen
    retValue.flag = template.flag
    retValue.mapq = template.mapq
    retValue.cigar = template.cigar
    retValue.rnext = template.rnext
    retValue.pnext = template.pnext
    retValue.tlen = template.tlen
    retValue.tags = template.tags
    retValue.reference_id = template.reference_id
    return retValue


cdef class PairwiseAlignmentFile:
    def __init__(self, object initializer):
        if(isinstance(initializer, str)):
            self.handle = pysam.AlignmentFile(initializer)
        elif(isinstance(initializer, pysam.calignmentfile.AlignmentFile)):
            self.handle = initializer

    def __next__(self):
        cdef tuple reads
        cdef AlignedSegment_t read1, read2, read
        while(len(reads) < 2):
            read = self.handle.next()
            if(read.is_supplementary or read.is_secondary):
                continue
            if(read.is_read1 is True):
                read1 = read
                continue
            elif(read.is_read2 is True):
                read2 = read
            reads = (read1, read2)
        assert reads[0].qname == reads[1].qname
        return reads

    def __iter__(self):
        return self