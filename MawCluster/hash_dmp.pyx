cdef class KFWrapper:
    def __cinit__(self, int readlen):
        cdef KingFisher_t tmp
        tmp = init_kf(readlen)
        self.src = &tmp

    def __dealloc__(self):
        destroy_kf(self.src)

cdef hash_dmp_write(cython.str infname, cython.str outfname):
    cdef dict hashmaster = {}
    cdef FILE *in_handle
    cdef FILE *out_handle
    cdef gzFile fp
    cdef kseq_t *seq
    cdef int l
    cdef char *bs_ptr
    cdef int blen, readlen
    cdef bytes infname_buf = <bytes>infname  # Temp variable to convert
    cdef bytes outfname_buf = <bytes>outfname # As above
    cdef KFWrapper tmp_wrapper
    cdef int *nuc_indices = <int *>malloc(2 * sizeof(int))
    if(infname == '-' or infname is None):
        in_handle = stdin
    else:
        in_handle = fopen(infname_buf, "r")
    if(outfname is None):
        out_handle = stdout
    else:
        out_handle = fopen(outfname_buf, "w")
    fp = gzdopen(fileno(in_handle), "r")
    seq = kseq_init(fp)
    l = kseq_read(seq)
    readlen = seq.seq.l
    bs_ptr = barcode_mem_view(seq)
    blen = infer_barcode_length(bs_ptr)
    tmp_wrapper = KFWrapper(readlen)
    pushback_kseq(tmp_wrapper.src, seq, nuc_indices, blen)
    hashmaster[get_binnerl(bs_ptr, blen)] = tmp_wrapper
    fclose(in_handle)
    fclose(out_handle)
    return
