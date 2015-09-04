from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf, sprintf, fprintf, stdin, stdout, fopen, fclose
from posix.stdio cimport fdopen, fileno
cimport cython


cdef extern from "zlib.h" nogil:
    ctypedef void * gzFile

    int gzclose(gzFile fp)
    int gzread(gzFile fp, void *buf, unsigned int n)
    char *gzerror(gzFile fp, int *errnum)

    gzFile gzopen( char *path, char *mode)
    gzFile gzdopen (int fd, char *mode)
    char * gzgets(gzFile file, char *buf, int len)
    int gzeof(gzFile file)

cdef extern from "../src/dmp/include/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s

cdef extern from "../src/dmp/dmp_interface.h" nogil:
    struct KingFisher:
        int **nuc_counts
        double **phred_sums
        int length
        int readlen
        char *max_phreds
        char *barcode
        int pass_fail
    ctypedef KingFisher KingFisher_t
    KingFisher_t init_kf(int readlen)
    void destroy_kf(KingFisher_t *kfp)
    int kseq_read(kseq_t *seq)

    ctypedef struct kseq_t:
        kstring_t name
        kstring_t comment
        kstring_t seq
        kstring_t qual

    kseq_t *kseq_init(gzFile)
    int kseq_read(kseq_t *)
    void kseq_destroy(kseq_t *)
    void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen)
    char *barcode_mem_view(kseq_t *seq)
    int infer_barcode_length(char *bs_ptr)
    uint64_t get_binnerl(char *barcode, int length)

cdef class KFWrapper:
    """
    Simple cython wrapper for KingFisher_t
    """
    cdef KingFisher_t *src
