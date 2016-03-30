#include "bmf_hashdmp.h"

namespace BMF {

    namespace {
        size_t buf_set_size = 50000;
    }

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



    int hash_dmp_main(int argc, char *argv[])
    {
        if(argc == 1) hashdmp_usage(), exit(EXIT_SUCCESS);
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


    void duplex_hash_process(kingfisher_hash_t *hfor, kingfisher_hash_t *cfor, kingfisher_hash_t *tmp_hkf, kingfisher_hash_t *crev, kingfisher_hash_t *hrev, gzFile out_handle, tmpvars_t *tmp)
    {
        kstring_t ks = {0, 0, nullptr};
        size_t buf_record_count = 0;

        HASH_ITER(hh, hfor, cfor, tmp_hkf) {
            if(buf_record_count++ == buf_set_size) {
                buf_record_count = 0;
                gzputs(out_handle, (const char *)ks.s);
                ks.l = 0;
            }
            HASH_FIND_STR(hrev, cfor->id, crev);
            if(crev) {
                zstranded_process_write(cfor->value, crev->value, &ks, tmp->buffers); // Found from both strands!
                destroy_kf(cfor->value), destroy_kf(crev->value);
                HASH_DEL(hrev, crev); HASH_DEL(hfor, cfor);
                free(crev); free(cfor);
            } else {
                dmp_process_write(cfor->value, &ks, tmp->buffers, 0); // No reverse strand found. \='{
                destroy_kf(cfor->value);
                HASH_DEL(hfor, cfor);
                free(cfor);
            }
        }
        if(buf_record_count) gzputs(out_handle, (const char *)ks.s);
        free(ks.s);
    }

    void hash_dmp_core(char *infname, char *outfname, int level)
    {
        char mode[4]("wT");
        if(level > 0) sprintf(mode, "wb%i", level % 10);
        FILE *in_handle(dlib::open_ifp(infname));
        gzFile out_handle(gzopen(outfname, mode));
        if(!in_handle) {
            if(dlib::isfile(infname)) {
                LOG_DEBUG("It's a file, but it's empty....\n");
                gzclose(out_handle);
                return;
            }
            LOG_EXIT("[E:%s] Could not open %s for reading. Abort mission!\n", infname);
            exit(EXIT_FAILURE);
        }
        gzFile fp(gzdopen(fileno(in_handle), "r"));
        kseq_t *seq(kseq_init(fp));
        // Initialized kseq
        int l = kseq_read(seq);
        if(l < 0) {
            fprintf(stderr, "[E:%s]: Could not open fastq file (%s). Abort mission!\n",
                    __func__, strcmp(infname, "-") == 0 ? "stdin": infname);
            exit(EXIT_FAILURE);
        }
        char *bs_ptr(barcode_mem_view(seq));
        const int blen(infer_barcode_length(bs_ptr));
    #if !NDEBUG
        fprintf(stderr, "[D:%s] Barcode length (inferred): %i.\n", __func__, blen);
    #endif
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
        fprintf(stderr, "[%s::%s] Loaded all records into memory. Writing out to file!\n", __func__, ifn_stream(infname));
        count = 0;
        uint64_t buf_record_count(0);
        kstring_t ks{0, 0, nullptr};
        HASH_ITER(hh, hash, current_entry, tmp_hk) {
            if(buf_record_count++ == buf_set_size) {
                count += buf_record_count;
                buf_record_count = 0;
                gzputs(out_handle, (const char *)ks.s);
                ks.l = 0;
            }
            dmp_process_write(current_entry->value, &ks, tmp->buffers, 0); // No reverse strand found. \='{
            destroy_kf(current_entry->value);
            HASH_DEL(hash, current_entry);
            free(current_entry);
        }
        count += buf_record_count;
        buf_record_count = 0;
        gzputs(out_handle, (const char *)ks.s);
        // Demultiplex and write out.
        fprintf(stderr, "[%s::%s] Total number of collapsed observations: %lu.\n", __func__, ifn_stream(infname), count);
        free(ks.s);
        gzclose(fp);
        gzclose(out_handle);
        kseq_destroy(seq);
        tmpvars_destroy(tmp);
    }

    void stranded_hash_dmp_core(char *infname, char *outfname, int level)
    {
        char mode[4] = "wT"; // Defaults to uncompressed "transparent" gzip output.
        if(level > 0) sprintf(mode, "wb%i", level % 10);
        LOG_DEBUG("Writing stranded hash dmp information with mode: '%s'.\n", mode);
        gzFile out_handle(gzopen(outfname, mode));
        gzFile fp(gzopen((infname && *infname) ? infname: "-", "r"));
        if(!fp) {
            fprintf(stderr, "[E:%s] Could not open %s for reading. Abort mission!\n", __func__, infname);
            exit(EXIT_FAILURE);
        }
        kseq_t *seq(kseq_init(fp));
        // Initialized kseq
        int l(kseq_read(seq));
        if(l < 0) {
            fprintf(stderr, "[%s]: Could not open fastq file (%s). Abort mission!\n",
                    __func__, strcmp(infname, "-") == 0 ? "stdin": infname);
            exit(EXIT_FAILURE);
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
            if(UNLIKELY(++count % 1000000 == 0))
                fprintf(stderr, "[%s::%s] Number of records processed: %lu.\n", __func__,
                        *infname == '-' ? "stdin" : infname, count);
            if(seq->comment.s[HASH_DMP_OFFSET] == 'F') {
                ++fcount;
                cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
                HASH_FIND_STR(hfor, tmp->key, tmp_hkf);
                if(!tmp_hkf) {
                    tmp_hkf = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                    tmp_hkf->value = init_kfp(tmp->readlen);
                    cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkf->id);
                    pushback_kseq(tmp_hkf->value, seq, blen);
                    HASH_ADD_STR(hfor, id, tmp_hkf);
                } else pushback_kseq(tmp_hkf->value, seq, blen);
            } else {
                cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
                HASH_FIND_STR(hrev, tmp->key, tmp_hkr);
                if(!tmp_hkr) {
                    tmp_hkr = (kingfisher_hash_t *)malloc(sizeof(kingfisher_hash_t));
                    tmp_hkr->value = init_kfp(tmp->readlen);
                    cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkr->id);
                    pushback_kseq(tmp_hkr->value, seq, blen);
                    HASH_ADD_STR(hrev, id, tmp_hkr);
                } else pushback_kseq(tmp_hkr->value, seq, blen);
            }
        }
        uint64_t rcount = count - fcount;
        LOG_INFO("Number of reverse reads: %lu. Number of forward reads: %lu.\n", rcount, fcount);
#if !NDEBUG
        fprintf(stderr, "[%s::%s] Loaded all records into memory. Writing out to file!\n", __func__, ifn_stream(outfname));
#endif
        // Write out all unmatched in forward and handle all barcodes handled from both strands.
        uint64_t duplex = 0, non_duplex = 0, non_duplex_fm = 0;
        size_t buf_record_count = 0;
        kstring_t ks = {0, 0, nullptr};
        // Demultiplex and empty the hash.
        HASH_ITER(hh, hfor, cfor, tmp_hkf) {
            if(buf_record_count++ == buf_set_size) {
                buf_record_count = 0;
                gzputs(out_handle, (const char *)ks.s);
                ks.l = 0;
            }
            HASH_FIND_STR(hrev, cfor->id, crev);
            if(!crev) {
                ++non_duplex;
                if(cfor->value->length > 1) ++non_duplex_fm;
                dmp_process_write(cfor->value, &ks, tmp->buffers, 0); // No reverse strand found. \='{
                destroy_kf(cfor->value);
                HASH_DEL(hfor, cfor);
                free(cfor);
            } else {
                ++duplex;
                zstranded_process_write(cfor->value, crev->value, &ks, tmp->buffers); // Found from both strands!
                destroy_kf(cfor->value); destroy_kf(crev->value);
                HASH_DEL(hrev, crev); HASH_DEL(hfor, cfor);
                cond_free(crev); cond_free(cfor);
            }
        }
        gzputs(out_handle, (const char *)ks.s);
        LOG_DEBUG("Before handling reverse only counts for non_duplex: %lu.\n", non_duplex);
        buf_record_count = 0;
        HASH_ITER(hh, hrev, crev, tmp_hkr) {
            if(buf_record_count++ == buf_set_size) {
                buf_record_count = 0;
                gzputs(out_handle, (const char *)ks.s);
                ks.l = 0;
            }
            ++non_duplex;
            if(crev->value->length > 1) ++non_duplex_fm;
            dmp_process_write(crev->value, &ks, tmp->buffers, 0); // No reverse strand found. \='{
            destroy_kf(crev->value);
            HASH_DEL(hrev, crev);
            free(crev);
        }
        gzputs(out_handle, (const char *)ks.s);
        LOG_DEBUG("Cleaning up.\n");
        LOG_INFO("Number of duplex observations: %lu. Number of non-duplex observations: %lu. Non-duplex families: %lu\n", duplex, non_duplex, non_duplex_fm);
        free(ks.s);
        gzclose(fp); gzclose(out_handle);
        kseq_destroy(seq);
        tmpvars_destroy(tmp);
    }

} /* namespace BMF */
