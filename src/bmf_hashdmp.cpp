#include <assert.h>
#include "bmf_dmp.h"
#include "dlib/io_util.h"
#include "lib/mseq.h"


namespace bmf {
void hash_inmem_inline_core(char *in1, char *in2, char *out1, char *out2,
                            char *homing, int blen, int threshold, int level=0, int mask=0,
                            int max_blen=-1);


void hashdmp_usage() {
    fprintf(stderr,
                    "Molecularly demultiplexes marked temporary fastqs into final unique observation records.\n"
                    "bmftools hashdmp does so in one large hashmap. This may require huge amounts of memory.\n"
                    "For samples which can't fit into memory, use bmftools dmp or sdmp, which subsets the sample,"
                    " reducing memory requirements while keeping the constant-time hashmap collapsing.\n"
                    "Usage: bmftools hashdmp <opts> <input_filename>.\n"
                    "Flags:\n"
                    "-s\tPerform secondary index consolidation rather than Loeb-like inline consolidation.\n"
                    "-o\tOutput filename.\n"
                    "If output file is unset, defaults to stdout. If input filename is not set, defaults to stdin.\n"
            );
}
void inmem_usage() {
    fprintf(stderr,
                    "Molecularly demultiplexes marked temporary fastqs into final unique observation records.\n"
                    "bmftools hashdmp does so in one large hashmap. This may require huge amounts of memory.\n"
                    "For samples which can't fit into memory, use bmftools dmp or sdmp, which subsets the sample,"
                    " reducing memory requirements while keeping the constant-time hashmap collapsing.\n"
                    "Usage: bmftools inmem <opts> -1 <out.r1.fastq> -2 <out.r2.fastq> <r1.fastq> <r2.fastq>.\n"
                    "Flags:\n"
                    "-s:\tHoming sequence -- REQUIRED.\n"
                    "-1:\tPath to output fastq for read 1. REQUIRED.\n"
                    "-2:\tPath to output fastq for read 2. REQUIRED.\n"
                    "-l:\tBarcode length. If using variable-length barcodes, this is the minimum barcode length.\n"
                    "-t:\tHomopolymer failure threshold. A molecular barcode with"
                    " a homopolymer of length >= this limit is flagged as QC fail. Default: 10.\n"
                    "-v:\tMaximum barcode length. (Set only if using variable-length barcodes.)\n"
                    "-m:\tSkip the first <INT> bases from each inline barcode. Default: 0\n"
                    "-L:\tOutput fastq compression level (Default: plain text).\n"
                    "If output file is unset, defaults to stdout. If input filename is not set, defaults to stdin.\n"
            );
}

tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen)
{
    tmpvars_t *ret((tmpvars_t *)calloc(1, sizeof(tmpvars_t)));
    ret->blen = blen;
    ret->readlen = readlen;
    ret->bs_ptr = bs_ptr;
    ret->buffers = (tmpbuffers_t *)malloc(sizeof(tmpbuffers_t));
    ret->buffers->name_buffer[0] = '@';
    ret->buffers->name_buffer[blen] = '\0';
    ret->buffers->cons_seq_buffer[readlen] = '\0';
    return ret;
}



int hashdmp_main(int argc, char *argv[])
{
    if(argc == 1) hashdmp_usage(), exit(EXIT_FAILURE);
    char *outfname = nullptr, *infname = nullptr;
    int c;
    int stranded_analysis = 1;
    int level = -1;
    while ((c = getopt(argc, argv, "l:o:sh?")) >= 0) {
        switch(c) {
            case 'l': level = atoi(optarg)%10; break;
            case 'o': outfname = optarg; break;
            case 's': stranded_analysis = 0; break;
            case '?': case 'h': hashdmp_usage(); return EXIT_SUCCESS;
        }
    }
    if(argc < 2) {
        fprintf(stderr, "[E:%s] Required arguments missing. See usage.\n", __func__);
        hashdmp_usage();
        exit(EXIT_FAILURE);
    }
    if(argc - 1 == optind) infname = argv[optind];
    else LOG_WARNING("Note: no input filename provided. Defaulting to stdin.\n");
    stranded_analysis ? stranded_hash_dmp_core(infname, outfname, level)
                      : hash_dmp_core(infname, outfname, level);
    LOG_INFO("Successfully complete bmftools hashdmp!\n");
    return EXIT_SUCCESS;
}


int hashdmp_inmem_main(int argc, char *argv[])
{
    if(argc == 1) inmem_usage(), exit(EXIT_FAILURE);
    char *outfname1 = const_cast<char *>("-");
    char *outfname2 = const_cast<char *>("-");
    char *homing = nullptr;
    int c;
    int blen = -1;
    int max_blen = -1;
    int mask = 0;
    int threshold = 10;
    int level = 0; // uncompressed
    while ((c = getopt(argc, argv, "1:2:v:l:L:l:m:s:t:h?")) >= 0) {
        switch(c) {
            case '1': outfname1 = optarg; break;
            case '2': outfname2 = optarg; break;
            case 'm': mask = atoi(optarg); break;
            case 'v': max_blen = atoi(optarg); break;
            case 'l': blen = atoi(optarg); break;
            case 'L': level = atoi(optarg) % 10; break;
            case 's': homing = optarg; break;
            case 't': threshold = atoi(optarg); break;
            case '?': case 'h': inmem_usage(); return EXIT_SUCCESS;
        }
    }
    if(argc < 2) {
        fprintf(stderr, "[E:%s] Required arguments missing. See usage.\n", __func__);
        hashdmp_usage();
        exit(EXIT_FAILURE);
    }
    if(argc - 2 != optind) {
        if(argc - 1 != optind)
            LOG_EXIT("Require at least one input fastq.\n");
        LOG_EXIT("raise NotImplementedError(\"Single-end inmem not implemented.\")\n")
    }
    if(blen < 0) LOG_EXIT("Barcode length required.");
    if(!homing) LOG_EXIT("Homing sequence required.\n");
    if(strcmp(outfname1, outfname2) == 0) LOG_EXIT("read 1 and read 2 must be separate files. Abort!\n");

    hash_inmem_inline_core(argv[optind], argv[optind + 1], outfname1, outfname2,
                           homing, blen, threshold, level, mask,
                           max_blen);
    LOG_INFO("Successfully complete bmftools hashdmp!\n");
    return EXIT_SUCCESS;
}

namespace {
    class hashtmp_t {
        kstring_t barcode;
        ~hashtmp_t() {
            if(barcode.s) free(barcode.s);
        }
    };
}

inline int get_blen(char *seq, char *homing, int homing_len, int blen, int max_blen, int mask) {
    for(int i = blen; i <= max_blen; ++i)
        if(memcmp(seq + i, homing, homing_len) == 0)
            return i - mask;
    return -1;
}

void hash_inmem_inline_core(char *in1, char *in2, char *out1, char *out2,
                            char *homing, int blen, int threshold, int level, int mask,
                            int max_blen) {
    if(max_blen < 0) max_blen = blen;
    char mode[4];
#if ZLIB_VER_MAJOR <= 1 && ZLIB_VER_MINOR <= 2 && ZLIB_VER_REVISION < 5
#pragma message("Note: zlib version < 1.2.5 doesn't support transparent file writing. Writing uncompressed temporary gzip files by default.")
// If not set, zlib compresses all our files enormously.
sprintf(mode, level > 0 ? "wb%i": "wb0", level % 10);
#else
sprintf(mode, level > 0 ? "wb%i": "wT", level % 10);
#endif
    if(level > 0) {
        if(strcmp(out1, "-") && strcmp(strrchr(out1, '\0') - 3, ".gz") != 0) {
            LOG_WARNING("Output gzip compressed but filename not terminated with .gz. FYI\n");
        } else if(strcmp(out2, "-") && strcmp(strrchr(out2, '\0') - 3, ".gz") != 0) {
            LOG_WARNING("Output gzip compressed but filename not terminated with .gz. FYI\n");
        }
    } else {
        if(strcmp(out1, "-") && strcmp(strrchr(out1, '\0') - 3, ".gz") == 0) {
            LOG_WARNING("Output filename stats with .gz but output is not compressed. FYI\n");
        } else if(strcmp(out2, "-") && strcmp(strrchr(out2, '\0') - 3, ".gz") == 0) {
            LOG_WARNING("Output filename stats with .gz but output is not compressed. FYI\n");
        }
    }
    FILE *in_handle1(dlib::open_ifp(in1));
    FILE *in_handle2(dlib::open_ifp(in2));
    gzFile out_handle1(dlib::open_gzfile(out1, mode));
    gzFile out_handle2(dlib::open_gzfile(out2, mode));
    const int homing_len = strlen(homing);
    if(!in_handle1) {
        if(dlib::isfile(in1)) {
            LOG_DEBUG("%s a file, but it's empty....\n", in1);
            gzclose(out_handle1);
            return;
        }
        LOG_EXIT("Could not open %s for reading. Abort mission!\n", in1);
    }
    if(!in_handle2) {
         if(dlib::isfile(in2)) {
            LOG_DEBUG("%s a file, but it's empty....\n", in2);
            gzclose(out_handle2);
             return;
         }
         LOG_EXIT("Could not open %s for reading. Abort mission!\n", in2);
     }
    gzFile fp1(gzdopen(fileno(in_handle1), "r"));
    gzFile fp2(gzdopen(fileno(in_handle2), "r"));
    kseq_t *seq1(kseq_init(fp1));
    kseq_t *seq2(kseq_init(fp2));
    kingfisher_hash_t *hash1f = nullptr, *hash2f = nullptr, *hash1r = nullptr, *hash2r = nullptr;
    kingfisher_hash_t *ce1 = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
    kingfisher_hash_t *ce2 = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
    kingfisher_hash_t *tmp_hk1 = ce1, *tmp_hk2 = ce2; // Save the pointer location for later comparison.
    kstring_t barcode = {0, 32, (char *)malloc(32uL * sizeof(char))};
    unsigned blen1, blen2;
    unsigned offset1, offset2;
    char pass;
    size_t barcode_count{0};
    while(LIKELY(kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0)) {
        pass = 1;
        blen1 = get_blen(seq1->seq.s, homing, homing_len, blen, max_blen, mask);
        blen2 = get_blen(seq2->seq.s, homing, homing_len, blen, max_blen, mask);
        if(switch_test(seq1, seq2, mask)) {
            if(blen2 != (unsigned)-1) {
                memcpy(barcode.s, seq2->seq.s + mask, blen2);
            } else {
                pass = 0;
                blen2 = blen - mask;
                memset(barcode.s, 'N', blen2);
            }
            barcode.l = blen2;
            //barcode.s[barcode.l] = '\0';
            if(blen1 != (unsigned)-1) {
                while(blen1 + blen2 >= barcode.m) {
                    kroundup32(barcode.m);
                    barcode.s = (char *)realloc(barcode.s, barcode.m);
                }
                memcpy(barcode.s + barcode.l, seq1->seq.s + mask, blen2);
            } else {
                pass = 0;
                blen1 = blen - mask;
                memset(barcode.s + barcode.l, 'N', blen1);
            }
            barcode.l = barcode.l + blen1;
            barcode.s[barcode.l] = '\0';
            //LOG_DEBUG("Looking for barcode %s.\n", barcode.s);
            HASH_FIND_STR(hash1r, barcode.s, tmp_hk1);
            HASH_FIND_STR(hash2r, barcode.s, tmp_hk2);
            pass &= test_hp(barcode.s, threshold);
            assert(!tmp_hk1 == !tmp_hk2); // Make sure that both have the same keyset.
            offset1 = blen1 + homing_len + mask;
            offset2 = blen2 + homing_len + mask;
            if(!tmp_hk1) {
                if(UNLIKELY(++barcode_count % 1000000 == 0))
                    LOG_INFO("Number of unique barcodes loaded: %lu\n", barcode_count);
                tmp_hk1 = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                tmp_hk2 = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                tmp_hk1->value = init_kfp(seq2->seq.l - offset2);
                tmp_hk2->value = init_kfp(seq1->seq.l - offset1);
                memcpy(tmp_hk1->id, barcode.s, barcode.l);
                memcpy(tmp_hk2->id, barcode.s, barcode.l);
                tmp_hk1->id[barcode.l] = '\0';
                tmp_hk2->id[barcode.l] = '\0';
                memcpy(tmp_hk1->value->barcode + 1, barcode.s, barcode.l);
                tmp_hk1->value->barcode[0] = '@';
                tmp_hk1->value->barcode[barcode.l + 1] = '\0';
                memcpy(tmp_hk2->value->barcode + 1, barcode.s, barcode.l);
                tmp_hk2->value->barcode[0] = '@';
                tmp_hk2->value->barcode[barcode.l + 1] = '\0';
                pushback_inmem(tmp_hk2->value, seq1, offset1, pass);
                pushback_inmem(tmp_hk1->value, seq2, offset2, pass);
                HASH_ADD_STR(hash1r, id, tmp_hk1);
                HASH_ADD_STR(hash2r, id, tmp_hk2);
            } else {
                pushback_inmem(tmp_hk2->value, seq1, offset1, pass);
                pushback_inmem(tmp_hk1->value, seq2, offset2, pass);
            }
        } else {
            if(blen1 != (unsigned)-1) memcpy(barcode.s, seq1->seq.s + mask, blen1);
            else { // Fail!
                pass = 0;
                blen1 = blen - mask;
                memset(barcode.s, 'N', blen1);
            }
            barcode.l = blen1;
            barcode.s[barcode.l] = '\0';
            if(blen2 != (unsigned)-1) {
                while(blen1 + blen2 >= barcode.m) {
                    kroundup32(barcode.m);
                    barcode.s = (char *)realloc(barcode.s, barcode.m);
                }
                memcpy(barcode.s + barcode.l, seq2->seq.s + mask, blen2);
            } else {
                pass = 0;
                blen2 = blen - mask;
                memset(barcode.s + barcode.l, 'N', blen2);
            }
            barcode.l = barcode.l + blen2;
            barcode.s[barcode.l] = '\0';
            pass &= test_hp(barcode.s, threshold);
            offset1 = blen1 + homing_len + mask;
            offset2 = blen2 + homing_len + mask;
            HASH_FIND_STR(hash1f, barcode.s, tmp_hk1);
            if(!tmp_hk1) {
                if(UNLIKELY(++barcode_count % 1000000 == 0))
                    LOG_INFO("Number of unique barcodes loaded: %lu\n", barcode_count);
                // Create
                tmp_hk1 = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                tmp_hk2 = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                tmp_hk1->value = init_kfp(seq1->seq.l - offset1);
                tmp_hk2->value = init_kfp(seq2->seq.l - offset2);
                memcpy(tmp_hk1->id, barcode.s, barcode.l);
                memcpy(tmp_hk2->id, barcode.s, barcode.l);
                tmp_hk1->id[barcode.l] = '\0';
                tmp_hk2->id[barcode.l] = '\0';
                memcpy(tmp_hk1->value->barcode + 1, barcode.s, barcode.l);
                tmp_hk1->value->barcode[barcode.l + 1] = '\0';
                tmp_hk1->value->barcode[0] = '@';
                memcpy(tmp_hk2->value->barcode + 1, barcode.s, barcode.l);
                tmp_hk2->value->barcode[barcode.l + 1] = '\0';
                tmp_hk2->value->barcode[0] = '@';
                HASH_ADD_STR(hash1f, id, tmp_hk1);
                HASH_ADD_STR(hash2f, id, tmp_hk2);
                pushback_inmem(tmp_hk1->value, seq1, offset1, pass);
                pushback_inmem(tmp_hk2->value, seq2, offset2, pass);
            } else {
                HASH_FIND_STR(hash2f, barcode.s, tmp_hk2);
                assert(!!tmp_hk2);
                pushback_inmem(tmp_hk1->value, seq1, offset1, pass);
                pushback_inmem(tmp_hk2->value, seq2, offset2, pass);
            }
        }
    }
    free(barcode.s);
    fclose(in_handle1), in_handle1 = nullptr;
    fclose(in_handle2), in_handle1 = nullptr;
    gzclose(fp1), fp1 = nullptr;
    gzclose(fp2), fp2 = nullptr;
    kseq_destroy(seq1), seq1 = nullptr;
    kseq_destroy(seq2), seq2 = nullptr;
    LOG_DEBUG("Loaded all records into memory.\n");

    kstring_t ks1{0, 0, nullptr};
    kstring_t ks2{0, 0, nullptr};
    tmpbuffers_t tmp;
    kingfisher_hash_t *t2 = nullptr;
    HASH_ITER(hh, hash1f, ce1, tmp_hk1) {
        HASH_FIND_STR(hash1r, ce1->id, t2);
        HASH_FIND_STR(hash2f, ce1->id, ce2);
        if(t2) {
            assert(ce2);
            HASH_FIND_STR(hash2r, ce1->id, tmp_hk2);
            assert(tmp_hk2);
            zstranded_process_write(ce1->value, t2->value, &ks1, &tmp);
            destroy_kf(t2->value);
            HASH_DEL(hash1r, t2);
            zstranded_process_write(ce2->value, tmp_hk2->value, &ks2, &tmp);
            destroy_kf(tmp_hk2->value);
            HASH_DEL(hash2r, tmp_hk2);
        } else {
            dmp_process_write(ce1->value, &ks1, &tmp, 0);
            dmp_process_write(ce2->value, &ks2, &tmp, 0);
        }
        destroy_kf(ce1->value);
        HASH_DEL(hash1f, ce1);
        destroy_kf(ce2->value);
        HASH_DEL(hash2f, ce2);
        gzputs(out_handle1, const_cast<const char *>(ks1.s));
        ks1.l = 0;
        gzputs(out_handle2, const_cast<const char *>(ks2.s));
        ks2.l = 0;
    }
    HASH_ITER(hh, hash1r, ce1, tmp_hk1) {
        HASH_FIND_STR(hash2r, ce1->id, ce2);
        dmp_process_write(ce1->value, &ks1, &tmp, 1);
        HASH_DEL(hash1r, ce1);
        gzputs(out_handle1, const_cast<const char *>(ks1.s));
        ks1.l = 0;
        destroy_kf(ce1->value);
        dmp_process_write(ce2->value, &ks2, &tmp, 1);
        gzputs(out_handle2, const_cast<const char *>(ks2.s));
        ks2.l = 0;
        destroy_kf(ce2->value);
        HASH_DEL(hash2r, ce2);
    }
    gzclose(out_handle1);
    gzclose(out_handle2);
    free(ks1.s);
    free(ks2.s);
}

void hash_dmp_core(char *infname, char *outfname, int level)
{
    char mode[4];
#if ZLIB_VER_MAJOR <= 1 && ZLIB_VER_MINOR <= 2 && ZLIB_VER_REVISION < 5
#pragma message("Note: zlib version < 1.2.5 doesn't support transparent file writing. Writing uncompressed temporary gzip files by default.")
// If not set, zlib compresses all our files enormously.
    sprintf(mode, level > 0 ? "wb%i": "wb0", level % 10);
#else
    sprintf(mode, level > 0 ? "wb%i": "wT", level % 10);
#endif
    LOG_DEBUG("zlib write mode: %s.\n", mode);
    FILE *in_handle(dlib::open_ifp(infname));
    gzFile out_handle(gzopen(outfname, mode));
    if(!in_handle) {
        if(dlib::isfile(infname)) {
            LOG_DEBUG("It's a file, but it's empty....\n");
            gzclose(out_handle);
            return;
        }
        LOG_EXIT("Could not open %s for reading. Abort mission!\n", infname);
    }
    gzFile fp(gzdopen(fileno(in_handle), "r"));
    kseq_t *seq(kseq_init(fp));
    // Initialized kseq
    int l = kseq_read(seq);
    if(l < 0) {
        gzclose(fp);
        fclose(in_handle);
        gzclose(out_handle);
        kseq_destroy(seq);
        return;
    }
    char *bs_ptr(barcode_mem_view(seq));
    const int blen(infer_barcode_length(bs_ptr));
    LOG_DEBUG("Barcode length (inferred): %i.\n", blen);
    tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
    memcpy(tmp->key, bs_ptr, blen);
    tmp->key[blen] = '\0';
    // Start hash table
    kingfisher_hash_t *hash = nullptr;
    kingfisher_hash_t *current_entry = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
    kingfisher_hash_t *tmp_hk = current_entry; // Save the pointer location for later comparison.
    cp_view2buf(bs_ptr + 1, current_entry->id);
    current_entry->value = init_kfp(tmp->readlen);
    HASH_ADD_STR(hash, id, current_entry);
    pushback_kseq(current_entry->value, seq, blen);

    uint64_t count = 0;
    // Add barcodes to the hash table
    while(LIKELY((l = kseq_read(seq)) >= 0)) {
        if(UNLIKELY(++count % 1000000 == 0))
            fprintf(stderr, "[%s::%s] Number of records read: %lu.\n", __func__,
                    strcmp("-", infname) == 0 ? "stdin": infname,count);
        cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
        HASH_FIND_STR(hash, tmp->key, tmp_hk);
        if(!tmp_hk) {
            tmp_hk = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
            tmp_hk->value = init_kfp(tmp->readlen);
            cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hk->id);
            pushback_kseq(tmp_hk->value, seq, blen);
            HASH_ADD_STR(hash, id, tmp_hk);
        } else pushback_kseq(tmp_hk->value, seq, blen);
    }
    LOG_DEBUG("Loaded all records into memory. Writing out to %s!\n", ifn_stream(outfname));
    count = 0;
    kstring_t ks{0, 0, nullptr};
    HASH_ITER(hh, hash, current_entry, tmp_hk) {
        ++count;
        dmp_process_write(current_entry->value, &ks, tmp->buffers, -1);
        gzputs(out_handle, (const char *)ks.s);
        ks.l = 0;
        destroy_kf(current_entry->value);
        HASH_DEL(hash, current_entry);
        free(current_entry);
    }
    // Demultiplex and write out.
#if !NDEBUG
    fprintf(stderr, "[D:%s::%s] Total number of collapsed observations: %lu.\n", __func__, ifn_stream(infname), count);
#endif
    free(ks.s);
    gzclose(fp);
    gzclose(out_handle);
    kseq_destroy(seq);
    tmpvars_destroy(tmp);
}
#if !NDEBUG
KHASH_MAP_INIT_INT(hd, uint64_t)
#endif

void stranded_hash_dmp_core(char *infname, char *outfname, int level)
{
#if !NDEBUG
    khash_t(hd) *hds = kh_init(hd);
#endif
    char mode[4] = "wT"; // Defaults to uncompressed "transparent" gzip output.
    if(level > 0) sprintf(mode, "wb%i", level % 10);
    LOG_DEBUG("Writing stranded hash dmp information with mode: '%s'.\n", mode);
    gzFile out_handle(gzopen(outfname, mode));
    gzFile fp(gzopen((infname && *infname) ? infname: "-", "r"));
    if(!fp) {
        LOG_EXIT("Could not open %s for reading. Abort mission!\n", infname);
    }
    if(!out_handle) {
        LOG_EXIT("Could not open %s for reading. Abort mission!\n", outfname);
    }
    kseq_t *seq(kseq_init(fp));
    // Initialized kseq
    int l(kseq_read(seq));
    if(l < 0) {
        gzclose(out_handle);
        gzclose(fp);
        kseq_destroy(seq);
        return;
    }
    char *bs_ptr = barcode_mem_view(seq);
    int blen = infer_barcode_length(bs_ptr);
    LOG_DEBUG("Barcode length (inferred): %i. First barcode: %s.\n", blen, bs_ptr);
    tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
    memcpy(tmp->key, bs_ptr, blen);
    tmp->key[blen] = '\0';
    // Start hash table
    kingfisher_hash_t *hfor = nullptr, *hrev = nullptr; // Hash forward, hash reverse
    kingfisher_hash_t *crev = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t)); // Current reverse, current forward.
    kingfisher_hash_t *cfor = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
    kingfisher_hash_t *tmp_hkr = crev, *tmp_hkf = cfor;
    uint64_t count = 1, fcount = 0;
    /* Handle first record.
     * We read in a record from the fastq to get the length of the reads
     * and the barcodes.
    */
    if(*bs_ptr == 'F') {
        ++fcount;
        cp_view2buf(bs_ptr + 1, cfor->id);
        cfor->value = init_kfp(tmp->readlen);
        HASH_ADD_STR(hfor, id, cfor);
        pushback_kseq(cfor->value, seq, blen);
    } else {
        cp_view2buf(bs_ptr + 1, crev->id);
        crev->value = init_kfp(tmp->readlen);
        HASH_ADD_STR(hrev, id, crev);
        pushback_kseq(hrev->value, seq, blen);
    }

    // Add reads to the hash
    while(LIKELY((l = kseq_read(seq)) >= 0)) {
#if !NDEBUG
        if(UNLIKELY(++count % 1000000 == 0))
            fprintf(stderr, "[%s::%s] Number of records processed: %lu.\n", __func__,
                    *infname == '-' ? "stdin" : infname, count);
#else
        ++count;
#endif
        if(seq->comment.s[HASH_DMP_OFFSET] == 'F') {
            ++fcount;
            cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
            HASH_FIND_STR(hfor, tmp->key, tmp_hkf);
            if(tmp_hkf) pushback_kseq(tmp_hkf->value, seq, blen);
            else {
                tmp_hkf = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                tmp_hkf->value = init_kfp(tmp->readlen);
                cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkf->id);
                pushback_kseq(tmp_hkf->value, seq, blen);
                HASH_ADD_STR(hfor, id, tmp_hkf);
            }
        } else {
            cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
            HASH_FIND_STR(hrev, tmp->key, tmp_hkr);
            if(tmp_hkr) pushback_kseq(tmp_hkr->value, seq, blen);
            else {
                tmp_hkr = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                tmp_hkr->value = init_kfp(tmp->readlen);
                cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkr->id);
                pushback_kseq(tmp_hkr->value, seq, blen);
                HASH_ADD_STR(hrev, id, tmp_hkr);
            }
        }
    }
#if !NDEBUG
    const uint64_t rcount = count - fcount;
#endif
    LOG_DEBUG("Number of reverse reads: %lu. Number of forward reads: %lu.\n", rcount, fcount);
    LOG_DEBUG("Loaded all records into memory. Writing out to %s!\n", ifn_stream(outfname));
    // Write out all unmatched in forward and handle all barcodes handled from both strands.
    uint64_t duplex = 0, non_duplex = 0, non_duplex_fm = 0;
    kstring_t ks{0, 0, nullptr};
    // Demultiplex and empty the hash.
#if !NDEBUG
    khiter_t ki;
    int hamming_distance, khr;
#endif
    HASH_ITER(hh, hfor, cfor, tmp_hkf) {
        HASH_FIND_STR(hrev, cfor->id, crev);
        if(crev) {
#if !NDEBUG
            hamming_distance = kf_hamming(cfor->value, crev->value);
            if((ki = kh_get(hd, hds, hamming_distance)) == kh_end(hds)) {
                ki = kh_put(hd, hds, hamming_distance, &khr);
                kh_val(hds, ki) = 1;
            } else ++kh_val(hds, ki);
#endif
            ++duplex;
            zstranded_process_write(cfor->value, crev->value, &ks, tmp->buffers); // Found from both strands!
            destroy_kf(cfor->value); destroy_kf(crev->value);
            HASH_DEL(hrev, crev); HASH_DEL(hfor, cfor);
            if(crev) free(crev), crev = nullptr;
            if(cfor) free(cfor), cfor = nullptr;

        } else {
            ++non_duplex;
            if(cfor->value->length > 1) ++non_duplex_fm;
            dmp_process_write(cfor->value, &ks, tmp->buffers, 0); // No reverse strand found. \='{
            destroy_kf(cfor->value);
            gzputs(out_handle, (const char *)ks.s);
            ks.l = 0;
            HASH_DEL(hfor, cfor);
            free(cfor);
        }
    }
#if !NDEBUG
    fprintf(stderr, "#HD\tCount\n");
    for(khiter_t ki = kh_begin(hds); ki != kh_end(hds); ++ki)
        if(kh_exist(hds, ki))
            fprintf(stderr, "%i\t%lu\n", kh_key(hds, ki), kh_val(hds, ki));
    kh_destroy(hd, hds);
#endif
    LOG_DEBUG("Before handling reverse only counts for non_duplex: %lu.\n", non_duplex);
    HASH_ITER(hh, hrev, crev, tmp_hkr) {
        ++non_duplex;
        if(crev->value->length > 1) ++non_duplex_fm;
        dmp_process_write(crev->value, &ks, tmp->buffers, 1); // Only reverse strand found. \='{
        destroy_kf(crev->value);
        gzputs(out_handle, (const char *)ks.s);
        ks.l = 0;
        HASH_DEL(hrev, crev);
        free(crev);
    }
    LOG_DEBUG("Cleaning up.\n");
    LOG_DEBUG("Number of duplex observations: %lu.\t"
              "Number of non-duplex observations: %lu.\t"
              "Non-duplex families: %lu\n",
              duplex, non_duplex, non_duplex_fm);
    free(ks.s);
    gzclose(fp); gzclose(out_handle);
    kseq_destroy(seq);
    tmpvars_destroy(tmp);
}

} /* namespace bmf */
