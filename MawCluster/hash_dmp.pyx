# distutils: sources = src/dmp/igamc_cephes.c
# cython: c_string_type=str, c_string_encoding=ascii

cdef class KFWrapper:
    def __init__(self, int readlen):
        cdef KingFisher_t tmp
        tmp = init_kf(readlen)
        self.src = &tmp

    def __dealloc__(self):
        destroy_kf(self.src)
'''
cdef hash_dmp_write(cython.str infname, cython.str outfname):
    fprintf(stderr, "You're using code that doesn't work written by "
                    "someone who isn't excited about it. Abort mission!\n")
    exit(1)
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
    cdef uint64_t bin
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
    cdef cython.str wtf = seq.seq.s
    print("WTF: %s" % wtf)
    readlen = len(wtf)
    fprintf(stderr, "OMGZ THE SEQ IS %s\n", seq.seq.s)
    bs_ptr = barcode_mem_view(seq)
    blen = infer_barcode_length(bs_ptr)
    tmp_wrapper = KFWrapper(readlen)
    pushback_kseq(tmp_wrapper.src, seq, nuc_indices, blen)
    hashmaster[get_binnerl(bs_ptr, blen)] = tmp_wrapper
    while True:
        l = kseq_read(seq)
        if(l < 0):
            break
        bin = get_binnerl(bs_ptr, blen)
        bs_ptr = barcode_mem_view(seq)
        try:
            tmp_wrapper = hashmaster[bin]
            pushback_kseq(tmp_wrapper.src, seq, nuc_indices, blen)
        except KeyError:
            tmp_wrapper = KFWrapper(readlen)
            pushback_kseq(tmp_wrapper.src, seq, nuc_indices, blen)
            hashmaster[bin] = tmp_wrapper
    cdef KFWrapper flavaflav
    for flavaflav in hashmaster.itervalues():
        dmp_process_write(flavaflav.src, out_handle, blen)
    fclose(in_handle)
    fclose(out_handle)
    return


cpdef hash_in_a_tin(cython.str infname, cython.str outfname):
    hash_dmp_write(infname, outfname)
    return
'''
