#include "bmf_err.h"

namespace BMF {

    namespace {
        uint64_t default_min_obs = 10000uL;
    }

    int err_main_main(int argc, char *argv[]);
    int err_fm_main(int argc, char *argv[]);
    int err_cycle_main(int argc, char *argv[]);
    int err_region_main(int argc, char *argv[]);

    RegionErr::RegionErr(region_set_t set, int i):
            counts({0}),
            name("")
        {
        LOG_DEBUG("Starting to make RegionErr for region_set_t with contig name at pos %p.\n", (void *)set.contig_name);
        kstring_t tmp{0, 0, nullptr};
        LOG_DEBUG("Contig name: %s. Strlen: %lu.\n", set.contig_name, strlen(set.contig_name));
        ksprintf(&tmp, "%s:%i:%i", set.contig_name, get_start(set.intervals[i]), get_stop(set.intervals[i]));
        name = std::string(tmp.s);
        free(tmp.s);
    }

    uint64_t get_max_obs(khash_t(obs) *hash)
    {
        uint64_t ret = 0;
        for(khiter_t k = kh_begin(hash); k != kh_end(hash); ++k)
            if(kh_exist(hash, k))
                if(kh_val(hash, k).obs > ret)
                    ret = kh_val(hash, k).obs;
        return ret;
    }


    uint64_t get_max_err(khash_t(obs) *hash)
    {
        uint64_t ret = 0;
        for(khiter_t k = kh_begin(hash); k != kh_end(hash); ++k)
            if(kh_exist(hash, k))
                if(kh_val(hash, k).err > ret)
                    ret = kh_val(hash, k).err;
        return ret;
    }

    int err_main_usage(int exit_status)
    {
        fprintf(stderr,
                        "Calculates error rates over a variety of variables."
                        "The primary output format consists of quality scores "
                        "Usage: bmftools err main <reference.fasta> <input.csrt.bam>\n"
                        "Flags:\n"
                        "-h/-?\t\tThis helpful help menu!\n"
                        "-o\t\tPath to output file. Set to '-' or 'stdout' to emit to stdout.\n"
                        "-a\t\tSet minimum mapping quality for inclusion.\n"
                        "-S\t\tSet minimum calculated PV tag value for inclusion.\n"
                        "-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
                        "-3:\t\tPath to write the 3d offset array in tabular format.\n"
                        "-f:\t\tPath to write the full measured error rates in tabular format.\n"
                        "-n:\t\tPath to write the cycle/nucleotide call error rates in tabular format.\n"
                        "-c:\t\tPath to write the cycle error rates in tabular format.\n"
                        "-b:\t\tPath to bed file for restricting analysis.\n"
                        "-m:\t\tMinimum family size for inclusion. Default: 0.\n"
                        "-M:\t\tMaximum family size for inclusion. Default: %i.\n"
                        "-d:\t\tFlag to only calculate error rates for duplex reads.\n"
                        "-D:\t\tFlag to only calculate error rates for non-duplex reads.\n"
                        "-p:\t\tSet padding for bed region. Default: %i.\n"
                        "-P:\t\tOnly include proper pairs.\n"
                        "-O:\t\tSet minimum number of observations for imputing quality Default: %lu.\n"
                , INT_MAX, DEFAULT_PADDING, default_min_obs);
        exit(exit_status);
        return exit_status;
    }

    int err_fm_usage(int exit_status)
    {
        fprintf(stderr,
                        "Calculates error rates by family size.\n"
                        "Usage: bmftools err fm -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
                        "Flags:\n"
                        "-h/-?\t\tThis helpful help menu!\n"
                        "-o\t\tPath to output file. Set to '-' or 'stdout' to emit to stdout.\n"
                        "-S\t\tSet minimum calculated PV tag value for inclusion.\n"
                        "-a\t\tSet minimum mapping quality for inclusion.\n"
                        "-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
                        "-b:\t\tPath to bed file for restricting analysis.\n"
                        "-d:\t\tFlag to only calculate error rates for duplex reads.\n"
                        "-p:\t\tSet padding for bed region. Default: %i.\n"
                        "-P:\t\tOnly include proper pairs.\n"
                        "-F:\t\tRequire that the FP tag be present and nonzero.\n"
                        "-f:\t\tRequire that the fraction of family members agreed on a  base be <FLOAT> or greater. Default: 0.0\n"
                , DEFAULT_PADDING);
        exit(exit_status);
        return exit_status; // This never happens.
    }

    int err_region_usage(int exit_status)
    {
        fprintf(stderr,
                        "Calculates error rates by genomic region.\n"
                        "Usage: bmftools err region <reference.fasta> <input.csrt.bam>\n"
                        "Flags:\n"
                        "-b\t\tPath to bed file. REQUIRED."
                        "-h/-?\t\tThis helpful help menu!\n"
                        "-o\t\tPath to output file. Leave unset or set to '-' or 'stdout' to emit to stdout.\n"
                        "-a\t\tSet minimum mapping quality for inclusion.\n"
                        "-p:\t\tSet padding for bed region. Default: %i.\n"
                , DEFAULT_PADDING);
        exit(exit_status);
        return exit_status; // This never happens.
    }


    int err_cycle_usage(int exit_status)
    {
        fprintf(stderr,
                        "Calculates error rate by cycle.\n"
                        "Usage: bmftools err cycle <opts> <reference.fasta> <input.csrt.bam>\n"
                        "Flags:\n"
                        "-h/-?\t\tThis helpful help menu!\n"
                        "-o\t\tPath to output file. Set to '-' or 'stdout' to emit to stdout.\n"
                        "-a\t\tSet minimum mapping quality for inclusion.\n"
                        "-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
                        "-b:\t\tPath to bed file for restricting analysis.\n"
                        "-p:\t\tSet padding for bed region. Default: 50.\n"
                        "-P:\t\tOnly include proper pairs.\n"
                );
        exit(exit_status);
        return exit_status; // This never happens.
    }


    void write_final(FILE *fp, fullerr_t *e)
    {
        for(uint32_t cycle = 0; cycle < e->l; ++cycle) {
            for(uint32_t qn = 0; qn < NQSCORES; ++qn) {
                fprintf(fp, "%i", e->r1->final[0][qn][cycle]);
                for(uint32_t bn = 1; bn < 4; ++bn)
                    fprintf(fp, ":%i", e->r1->final[bn][qn][cycle]);
                if(qn != NQSCORES - 1) fprintf(fp, ",");
            }
            fputc('|', fp);
            for(uint32_t qn = 0; qn < NQSCORES; ++qn) {
                fprintf(fp, "%i", e->r2->final[0][qn][cycle]);
                for(uint32_t bn = 1; bn < 4; ++bn)
                    fprintf(fp, ":%i", e->r2->final[bn][qn][cycle]);
                if(qn != NQSCORES - 1) fprintf(fp, ",");
            }
            fputc('\n', fp);
        }
    }

    void err_cycle_report(FILE *fp, cycle_err_t *ce)
    {
        LOG_DEBUG("Beginning err cycle report.\n");
        fprintf(fp, "#Cycle number\tRead 1 error\tRead 2 error\tRead 1 error counts"
                "\tRead 1 observation counts\tRead 2 error counts\tRead 2 observation counts\n");
        for(int i = 0; i < ce->rlen; ++i) {
            fprintf(fp, "%i\t%0.12f\t%lu\t%lu\t%0.12f\t%lu\t%lu\n",
                    i + 1,
                    (double)ce->r1[i].err / ce->r1[i].obs,
                    ce->r1[i].err, ce->r1[i].obs,
                    (double)ce->r2[i].err / ce->r2[i].obs, ce->r2[i].err, ce->r2[i].obs);
        }
    }

    void err_fm_report(FILE *fp, fmerr_t *f)
    {
        LOG_DEBUG("Beginning err fm report.\n");
        int khr, fm;
        khiter_t k, k1, k2;
        // Make a set of all FMs to print out.
        khash_t(obs_union) *key_union = kh_init(obs_union);
        for(k1 = kh_begin(f->hash1); k1 != kh_end(f->hash1); ++k1)
            if(kh_exist(f->hash1, k1))
                if((k = kh_get(obs_union, key_union, kh_key(f->hash1, k1))) == kh_end(key_union))
                    k = kh_put(obs_union, key_union, kh_key(f->hash1, k1), &khr);
        for(k2 = kh_begin(f->hash2); k2 != kh_end(f->hash2); ++k2)
            if(kh_exist(f->hash2, k2))
                if((k = kh_get(obs_union, key_union, kh_key(f->hash2, k2))) == kh_end(key_union))
                    k = kh_put(obs_union, key_union, kh_key(f->hash2, k2), &khr);

        // Write  header
        fprintf(fp, "##PARAMETERS\n##refcontig:\"%s\"\n##bed:\"%s\"\n"
                "##minMQ:%i\n##Duplex Required: %s.\n##Duplex Refused: %s.\n", f->refcontig ? f->refcontig: "N/A",
                f->bedpath? f->bedpath: "N/A", f->minMQ,
                f->flag & REQUIRE_DUPLEX ? "True": "False",
                f->flag & REFUSE_DUPLEX ? "True": "False");
        fprintf(fp, "##STATS\n##nread:%lu\n##nskipped:%lu\n", f->nread, f->nskipped);
        fprintf(fp, "#FM\tRead 1 Error\tRead 2 Error\tRead 1 Errors\tRead 1 Counts\tRead 2 Errors\tRead 2 Counts\n");
        int *tmp = (int *)malloc(key_union->n_occupied * sizeof(int));
        for(k = kh_begin(key_union), khr = 0; k != kh_end(key_union); ++k)
            if(kh_exist(key_union, k))
                tmp[khr++] = kh_key(key_union, k);
        std::sort(tmp, tmp + key_union->n_occupied);
        for(unsigned i = 0; i < key_union->n_occupied; ++i) {
            fm = tmp[i];
            fprintf(fp, "%i\t", fm);

            if((k1 = kh_get(obs, f->hash1, fm)) == kh_end(f->hash1))
                fprintf(fp, "-nan\t");
            else
                fprintf(fp, "%0.12f\t", (double)kh_val(f->hash1, k1).err / kh_val(f->hash1, k1).obs);

            if((k2 = kh_get(obs, f->hash2, fm)) == kh_end(f->hash2))
                fprintf(fp, "-nan\t");
            else
                fprintf(fp, "%0.12f\t", (double)kh_val(f->hash2, k2).err / kh_val(f->hash2, k2).obs);

            if(k1 != kh_end(f->hash1))
                fprintf(fp, "%lu\t%lu\t",
                    kh_val(f->hash1, k1).err, kh_val(f->hash1, k1).obs);
            else fputs("0\t0\t", fp);
            if(k2 != kh_end(f->hash2))
                fprintf(fp, "%lu\t%lu\n",
                    kh_val(f->hash2, k2).err, kh_val(f->hash2, k2).obs);
            else fputs("0\t0\n", fp);
        }
        free(tmp);
        kh_destroy(obs_union, key_union);
    }

    void err_report(FILE *fp, fullerr_t *f)
    {
        LOG_DEBUG("Beginning error main report.\n");
        fprintf(fp, "{\n{\"total_read\": %lu},\n{\"total_skipped\": %lu},\n", f->nread, f->nskipped);
        uint64_t n1_obs = 0, n1_err = 0, n1_ins = 0;
        uint64_t n2_obs = 0, n2_err = 0, n2_ins = 0;
        // n_ins is number with insufficient observations to report.
        for(int i = 0; i < 4; ++i) {
            for(unsigned j = 0; j < NQSCORES; ++j) {
                for(unsigned k = 0; k < f->l; ++k) {
                    n1_obs += f->r1->obs[i][j][k]; n1_err += f->r1->err[i][j][k];
                    n2_obs += f->r2->obs[i][j][k]; n2_err += f->r2->err[i][j][k];
                    if(f->r1->obs[i][j][k] < f->min_obs) ++n1_ins;
                    if(f->r2->obs[i][j][k] < f->min_obs) ++n2_ins;
                }
            }
        }
        uint64_t n_cases = NQSCORES * 4 * f->l;
        fprintf(stderr, "{\"read1\": {\"total_error\": %f},\n{\"total_obs\": %lu},\n{\"total_err\": %lu}"
                ",\n{\"number_insufficient\": %lu},\n{\"n_cases\": %lu}},",
                (double)n1_err / n1_obs, n1_obs, n1_err, n1_ins, n_cases);
        fprintf(stderr, "{\"read2\": {\"total_error\": %f},\n{\"total_obs\": %lu},\n{\"total_err\": %lu}"
                ",\n{\"number_insufficient\": %lu},\n{\"n_cases\": %lu}},",
                (double)n2_err / n2_obs, n2_obs, n2_err, n2_ins, n_cases);
        fprintf(fp, "}");
    }

    void readerr_destroy(readerr_t *e){
        for(int i = 0; i < 4; ++i) {
            for(unsigned j = 0; j < NQSCORES; ++j) {
                cond_free(e->obs[i][j]);
                cond_free(e->err[i][j]);
                cond_free(e->final[i][j]);
            }
            cond_free(e->obs[i]);
            cond_free(e->err[i]);
            cond_free(e->qobs[i]);
            cond_free(e->qerr[i]);
            cond_free(e->final[i]);
            cond_free(e->qpvsum[i]);
            cond_free(e->qdiffs[i]);
        }
        cond_free(e->obs);
        cond_free(e->err);
        cond_free(e->qerr);
        cond_free(e->qobs);
        cond_free(e->final);
        cond_free(e->qpvsum);
        cond_free(e->qdiffs);
        cond_free(e);
    }


    void err_fm_core(char *fname, faidx_t *fai, fmerr_t *f, htsFormat *open_fmt)
    {
        samFile *fp = sam_open_format(fname, "r", open_fmt);
        bam_hdr_t *hdr = sam_hdr_read(fp);
        if (!hdr) LOG_EXIT("Failed to read input header from bam %s. Abort!\n", fname);
        bam1_t *b = bam_init1();
        int32_t is_rev, ind, s, i, fc, rc, r, khr, DR, FP, FM, reflen, length, pos, tid_to_study = -1, last_tid = -1;
        char *ref = nullptr; // Will hold the sequence for a  chromosome
        khash_t(obs) *hash;
        uint8_t *seq, *drdata, *fpdata;
        uint32_t *cigar, *pv_array, *fa_array;
        khiter_t k;
        if(f->refcontig) {
            for(int i = 0; i < hdr->n_targets; ++i) {
                if(!strcmp(hdr->target_name[i], f->refcontig)) {
                    tid_to_study = i; break;
                }
            }
            if(tid_to_study < 0)
                LOG_EXIT("Contig %s not found in bam header. Abort mission!\n", f->refcontig);
        }
        while(LIKELY((r = sam_read1(fp, hdr, b)) != -1)) {
            if(++f->nread % 1000000 == 0) LOG_INFO("Records read: %lu.\n", f->nread);
            pv_array = (uint32_t *)dlib::array_tag(b, "PV");
            fa_array = (uint32_t *)dlib::array_tag(b, "FA");
            FM = bam_itag(b, "FM");
            drdata = bam_aux_get(b, "DR");
            fpdata = bam_aux_get(b, "FP");
            DR = dlib::int_tag_zero(drdata);
            FP = dlib::int_tag_zero(fpdata);
            // Pass reads without FP tag.
            if(b->core.flag & (BAM_FSECONDARY | BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP)) {
                ++f->nskipped;
                continue;
            }
            if(b->core.qual < f->minMQ) {
                ++f->nskipped;
                continue;
            }
            if(f->refcontig && tid_to_study != b->core.tid) {
                continue;
                ++f->nskipped;
            }
            if((f->flag & REQUIRE_PROPER) && (!(b->core.flag & BAM_FPROPER_PAIR))) {
                ++f->nskipped;
                continue;
            }
            if(f->bed && dlib::bed_test(b, f->bed) == 0) {
                ++f->nskipped;
                continue;
            }
            if((f->flag & REQUIRE_DUPLEX) && !DR) {
                ++f->nskipped;
                continue;
            }
            if((f->flag & REFUSE_DUPLEX) && DR) {
                ++f->nskipped;
                continue;
            }
            if((f->flag & REQUIRE_FP_PASS) && FP == 0) {
                ++f->nskipped;
                continue;
            }
            seq = (uint8_t *)bam_get_seq(b);
            cigar = bam_get_cigar(b);
            hash = (b->core.flag & BAM_FREAD1) ? f->hash1: f->hash2;
            if(b->core.tid != last_tid) {
                cond_free(ref);
                LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
                if((ref = fai_fetch(fai, hdr->target_name[b->core.tid], &reflen)) == nullptr)
                    LOG_EXIT("Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);
                LOG_DEBUG("Fetched.\n");
                last_tid = b->core.tid;
            }
            pos = b->core.pos;
            is_rev = b->core.flag & BAM_FREVERSE;
            //LOG_DEBUG("Max err count in hashmap: %lu.\n", get_max_err(hash));
            //LOG_DEBUG("Max obs count in hashmap: %lu.\n", get_max_obs(hash));
            if((k = kh_get(obs, hash, FM)) == kh_end(hash)) {
                k = kh_put(obs, hash, FM, &khr);
                memset(&kh_val(hash, k), 0, sizeof(obserr_t));
            }
            for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
                length = bam_cigar_oplen(cigar[i]);
                switch(bam_cigar_op(cigar[i])) {
                case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                    if(is_rev) {
                        for(ind = 0; ind < length; ++ind) {
                            if(pv_array[b->core.l_qseq - 1 - ind - rc] < f->minPV)
                                continue;
                            if(((double)fa_array[b->core.l_qseq - 1 - ind - rc] / FM) < f->minFR)
                                continue;
                            s = bam_seqi(seq, ind + rc);
                            //fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
                            if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                            ++kh_val(hash, k).obs;
                            if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
                                ++kh_val(hash, k).err;
                        }
                    } else {
                        for(ind = 0; ind < length; ++ind) {
                            if(pv_array[ind + rc] < f->minPV)
                                continue;
                            if(((double)fa_array[ind + rc] / FM) < f->minFR)
                                continue;
                            s = bam_seqi(seq, ind + rc);
                            //fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
                            if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                            ++kh_val(hash, k).obs;
                            if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
                                ++kh_val(hash, k).err;
                        }
                    }
                    rc += length; fc += length;
                    break;
                case BAM_CSOFT_CLIP: case BAM_CHARD_CLIP: case BAM_CINS:
                    rc += length;
                    break;
                case BAM_CREF_SKIP: case BAM_CDEL:
                    fc += length;
                    break;
                }
            }
        }
        LOG_INFO("Total records read: %lu. Total records skipped: %lu.\n", f->nread, f->nskipped);
        cond_free(ref);
        bam_destroy1(b);
        bam_hdr_destroy(hdr), sam_close(fp);
    }


    cycle_err_t *err_cycle_core(char *fname, faidx_t *fai, htsFormat *open_fmt,
                                char *bedpath, char *refcontig, int padding, unsigned minMQ, int flag)
    {
        bam1_t *b = bam_init1();
        samFile *fp = sam_open_format(fname, "r", open_fmt);
        bam_hdr_t *hdr = sam_hdr_read(fp);
        if (!hdr)
            LOG_EXIT("Failed to read input header from bam %s. Abort!\n", fname);
        const int32_t rlen = b->core.l_qseq;
        cycle_err_t *ce = cycle_init(bedpath, hdr, refcontig, padding, minMQ, rlen, flag);
        int32_t is_rev, ind, s, i, fc, rc, r, reflen, length, cycle, pos, tid_to_study = -1, last_tid = -1;
        uint8_t *seq, *fpdata;
        uint32_t *cigar;
        obserr_t *arr;
        char *ref = nullptr; // Will hold the sequence for a  chromosome
        if(ce->refcontig) {
            for(int i = 0; i < hdr->n_targets; ++i) {
                if(!strcmp(hdr->target_name[i], ce->refcontig)) {
                    tid_to_study = i; break;
                }
            }
            if(tid_to_study < 0)
                LOG_EXIT("Contig %s not found in bam header. Abort mission!\n", ce->refcontig);
        }
        while(LIKELY((r = sam_read1(fp, hdr, b)) != -1)) {
            if(++ce->nread % 1000000 == 0) LOG_INFO("Records read: %lu.\n", ce->nread);
            // Filter
            if((b->core.flag & 1796) /* unmapped, secondary, qc fail, duplicate*/||
                b->core.qual < ce->minMQ ||
                (ce->refcontig && tid_to_study != b->core.tid) ||
                ((ce->flag & REQUIRE_FP_PASS) && ((fpdata = bam_aux_get(b, "FP")) != nullptr) && bam_aux2i(fpdata) == 0) ||
                ((ce->flag & REQUIRE_PROPER) && (!(b->core.flag & BAM_FPROPER_PAIR))) ||
                (ce->bed && dlib::bed_test(b, ce->bed) == 0) /* Outside of region */) {
                ++ce->nskipped;
                /* LOG_DEBUG("Skipped record with name %s.\n", bam_get_qname(b)); */
                continue;
            }
            seq = (uint8_t *)bam_get_seq(b);
            cigar = bam_get_cigar(b);
            if(b->core.tid != last_tid) {
                cond_free(ref);
                LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
                ref = fai_fetch(fai, hdr->target_name[b->core.tid], &reflen);
                if(!ref) {
                    LOG_EXIT("Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);
                }
                last_tid = b->core.tid;
            }
            pos = b->core.pos;
            is_rev = b->core.flag & BAM_FREVERSE;
            arr = (b->core.flag & BAM_FREAD1) ? ce->r1: ce->r2;
            for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
                length = bam_cigar_oplen(cigar[i]);
                switch(bam_cigar_op(cigar[i])) {
                case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                    if(is_rev) {
                        for(ind = 0; ind < length; ++ind) {
                            s = bam_seqi(seq, ind + rc);
                            if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                            cycle = rlen - 1 - ind - rc;
                            ++arr[cycle].obs;
                            if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s) ++arr[cycle].err;
                        }
                    } else {
                        for(ind = 0; ind < length; ++ind) {
                            s = bam_seqi(seq, ind + rc);
                            if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                            cycle = ind + rc;
                            ++arr[cycle].obs;
                            if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s) ++arr[cycle].err;
                        }
                    }
                    rc += length; fc += length;
                    break;
                case BAM_CSOFT_CLIP: case BAM_CHARD_CLIP: case BAM_CINS:
                    rc += length;
                    break;
                case BAM_CREF_SKIP: case BAM_CDEL:
                    fc += length;
                    break;
                }
            }
        }
        LOG_INFO("Total records read: %lu. Total records skipped: %lu.\n", ce->nread, ce->nskipped);
        cond_free(ref);
        bam_destroy1(b);
        bam_hdr_destroy(hdr), sam_close(fp);
        return ce;
    }

    void err_main_core(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt)
    {
        if(!f->r1) f->r1 = readerr_init(f->l);
        if(!f->r2) f->r2 = readerr_init(f->l);
        samFile *fp = sam_open_format(fname, "r", open_fmt);
        bam_hdr_t *hdr = sam_hdr_read(fp);
        if (!hdr)
            LOG_EXIT("Failed to read input header from bam %s. Abort!\n", fname);
        int32_t i, s, c, len, pos, FM, RV, rc, fc, last_tid = -1, tid_to_study = -1, cycle, is_rev;
        unsigned ind;
        bam1_t *b = bam_init1();
        char *ref = nullptr; // Will hold the sequence for a  chromosome
        if(f->refcontig) {
            for(i = 0; i < hdr->n_targets; ++i) {
                if(!strcmp(hdr->target_name[i], f->refcontig)) {
                    tid_to_study = i; break;
                }
            }
            if(tid_to_study < 0) LOG_EXIT("Contig %s not found in bam header. Abort mission!\n", f->refcontig);
        }
        uint8_t *fdata, *rdata, *pdata, *seq, *qual;
        uint32_t *cigar, *pv_array, length;
        readerr_t *r;
        while(LIKELY((c = sam_read1(fp, hdr, b)) != -1)) {
            fdata = bam_aux_get(b, "FM");
            rdata = bam_aux_get(b, "RV");
            pdata = bam_aux_get(b, "FP");
            FM = dlib::int_tag_zero(fdata);
            RV = dlib::int_tag_zero(rdata);
            // Filters... WOOF
            if((b->core.flag & 1796) || b->core.qual < f->minMQ || (f->refcontig && tid_to_study != b->core.tid) ||
                (f->bed && dlib::bed_test(b, f->bed) == 0) || // Outside of region
                (FM < f->minFM) || (FM > f->maxFM) || // minFM
                ((f->flag & REQUIRE_PROPER) && (!(b->core.flag & BAM_FPROPER_PAIR))) || // skip improper pairs
                ((f->flag & REQUIRE_DUPLEX) ? (RV == FM || RV == 0): ((f->flag & REFUSE_DUPLEX) && (RV != FM && RV != 0))) || // Requires
                ((f->flag & REQUIRE_FP_PASS) && pdata && bam_aux2i(pdata) == 0) /* Fails barcode QC */) {
                    ++f->nskipped;
                    //LOG_DEBUG("Skipped record with name %s.\n", bam_get_qname(b));
                    continue;
            }
            seq = (uint8_t *)bam_get_seq(b);
            qual = (uint8_t *)bam_get_qual(b);
            cigar = bam_get_cigar(b);

            if(++f->nread % 1000000 == 0) LOG_INFO("Records read: %lu.\n", f->nread);
            if(b->core.tid != last_tid) {
                last_tid = b->core.tid;
                cond_free(ref);
                LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
                ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
                if(ref == nullptr) LOG_EXIT("[Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);
            }
            r = (b->core.flag & BAM_FREAD1) ? f->r1: f->r2;
            pos = b->core.pos;
            is_rev = (b->core.flag & BAM_FREVERSE);
            if(f->minPV) {
                pv_array = (uint32_t *)dlib::array_tag(b, "PV");
                for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
                    length = bam_cigar_oplen(cigar[i]);
                    switch(bam_cigar_op(cigar[i])) {
                    case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                        if(is_rev) {
                            for(ind = 0; ind < length; ++ind) {
                                cycle = b->core.l_qseq - 1 - ind - rc;
                                if(pv_array[cycle] < f->minPV)
                                    continue;
                                s = bam_seqi(seq, ind + rc);
                                //fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
                                if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                                ++r->obs[bamseq2i[s]][qual[ind + rc] - 2][cycle];
                                if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
                                    ++r->err[bamseq2i[s]][qual[ind + rc] - 2][cycle];
                            }
                        } else {
                            for(ind = 0; ind < length; ++ind) {
                                cycle = ind + rc;
                                if(pv_array[cycle] < f->minPV)
                                    continue;
                                s = bam_seqi(seq, cycle);
                                //fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
                                if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                                ++r->obs[bamseq2i[s]][qual[cycle] - 2][cycle];
                                if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
                                    ++r->err[bamseq2i[s]][qual[cycle] - 2][cycle];
                            }
                        }
                        rc += length; fc += length;
                        break;
                    case BAM_CSOFT_CLIP: case BAM_CHARD_CLIP: case BAM_CINS:
                        rc += length;
                        break;
                    case BAM_CREF_SKIP: case BAM_CDEL:
                        fc += length;
                        break;
                    }
                }
            } else {
                for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
                    length = bam_cigar_oplen(cigar[i]);
                    switch(bam_cigar_op(cigar[i])) {
                    case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                        if(is_rev) {
                            for(ind = 0; ind < length; ++ind) {
                                s = bam_seqi(seq, ind + rc);
                                if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                                cycle = b->core.l_qseq - 1 - ind - rc;
                                ++r->obs[bamseq2i[s]][qual[ind + rc] - 2][cycle];
                                if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
                                    ++r->err[bamseq2i[s]][qual[ind + rc] - 2][cycle];
                            }
                        } else {
                            for(ind = 0; ind < length; ++ind) {
                                cycle = ind + rc;
                                s = bam_seqi(seq, cycle);
                                if(s == dlib::htseq::HTS_N || ref[pos + fc + ind] == 'N') continue;
                                ++r->obs[bamseq2i[s]][qual[cycle] - 2][cycle];
                                if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
                                    ++r->err[bamseq2i[s]][qual[cycle] - 2][cycle];
                            }
                        }
                        rc += length; fc += length;
                        break;
                    case BAM_CSOFT_CLIP: case BAM_CHARD_CLIP: case BAM_CINS:
                        rc += length;
                        break;
                    case BAM_CREF_SKIP: case BAM_CDEL:
                        fc += length;
                        break;
                    }
                }
            }
        }
        cond_free(ref);
        bam_destroy1(b);
        bam_hdr_destroy(hdr), sam_close(fp);
    }


    void write_full_rates(FILE *fp, fullerr_t *f)
    {
        uint64_t l;
        unsigned i, j;
        for(l = 0; l < f->l; ++l) {
            for(j = 0; j < NQSCORES; ++j) {
                for(i = 0; i < 4u; ++i) {
                    if(f->r1->obs[i][j][l])
                        fprintf(fp, i ? ":%0.12f": "%0.12f", (double)f->r1->err[i][j][l] / f->r1->obs[i][j][l]);
                    else fputs(i ? ":-1337": "-1337", fp);
                }
                if(j != NQSCORES - 1) fputc(',', fp);
            }
            fputc('|', fp);
            for(j = 0; j < NQSCORES; ++j) {
                for(i = 0; i < 4u; ++i) {
                    if(f->r2->obs[i][j][l])
                        fprintf(fp, i ? ":%0.12f": "%0.12f", (double)f->r2->err[i][j][l] / f->r2->obs[i][j][l]);
                    else fputs(i ? ":-1337": "-1337", fp);
                }
                if(j != NQSCORES - 1) fputc(',', fp);
            }
            fputc('\n', fp);
        }
    }



    void write_base_rates(FILE *fp, fullerr_t *f)
    {
        fputs("#Cycle\tR1A\tR1C\tR1G\tR1T\tR2A\tR2C\tR2G\tR2T\n", fp);
        for(uint64_t l = 0; l < f->l; ++l) {
            int i;
            fprintf(fp, "%lu\t", l + 1);
            for(i = 0; i < 4; ++i) {
                LOG_DEBUG("obs: %lu. err: %lu.\n", f->r1->qerr[i][l], f->r1->qobs[i][l]);
                fprintf(fp, i ? "\t%0.12f": "%0.12f", (double)f->r1->qerr[i][l] / f->r1->qobs[i][l]);
            }
            fputc('|', fp);
            for(i = 0; i < 4; ++i)
                fprintf(fp, i ? "\t%0.12f": "%0.12f", (double)f->r2->qerr[i][l] / f->r2->qobs[i][l]);
            fputc('\n', fp);
        }
    }


    void write_global_rates(FILE *fp, fullerr_t *f)
    {
        fprintf(fp, "##Parameters: minFM %i. maxFM %i.", f->minFM, f->maxFM);
        fputs("Duplex required: ", fp);
        fputs((f->flag & REQUIRE_DUPLEX) ? "True": "False", fp);
        fputc('\n', fp);
        uint64_t sum1 = 0, sum2 = 0, counts1 = 0, counts2 = 0;
        for(uint64_t l = 0; l < f->l; ++l) {
            for(int i = 0; i < 4; ++i) {
                sum1 += f->r1->qerr[i][l];
                counts1 += f->r1->qobs[i][l];
                sum2 += f->r2->qerr[i][l];
                counts2 += f->r2->qobs[i][l];
            }
        }
        fprintf(fp, "#Global Error Rates\t%0.12f\t%0.12f\n", (double)sum1 / counts1, (double)sum2 / counts2);
        fprintf(fp, "#Global Sum/Count\t%lu/%lu\t%lu/%lu\n", sum1, counts1, sum2, counts2);
    }

    void write_cycle_rates(FILE *fp, fullerr_t *f)
    {
        fputs("#Cycle\tRead 1 Error Rate\tRead 2 Error Rate\tRead 1 Error Count\t"
                "Read 1 Obs Count\tRead 2 Error Count\tRead 2 Obs Count\n", fp);
        for(uint64_t l = 0; l < f->l; ++l) {
            fprintf(fp, "%lu\t", l + 1);
            uint64_t sum1 = 0, sum2 = 0, counts1 = 0, counts2 = 0;
            for(int i = 0; i < 4; ++i) {
                sum1 += f->r1->qerr[i][l];
                counts1 += f->r1->qobs[i][l];
                sum2 += f->r2->qerr[i][l];
                counts2 += f->r2->qobs[i][l];
            }
            fprintf(fp, "%0.12f\t%0.12f\t%lu\t%lu\t%lu\t%lu\n", (double)sum1 / counts1, (double)sum2 / counts2,
                    sum1, counts1, sum2, counts2);
        }
    }

    void impute_scores(fullerr_t *f)
    {
        for(unsigned i = 0; i < 4u; ++i) {
            for(uint64_t l = 0; l < f->l; ++l) {
                // Handle qscores of 2
                f->r1->final[i][0][l] = f->r1->obs[i][0][l] >= f->min_obs ? pv2ph((double)f->r1->err[i][0][l] / f->r1->obs[i][0][l])
                                                                          : 2;
                f->r2->final[i][0][l] = f->r2->obs[i][0][l] >= f->min_obs ? pv2ph((double)f->r2->err[i][0][l] / f->r2->obs[i][0][l])
                                                                          : 2;
                for(unsigned j = 1; j < NQSCORES; ++j) {
                    f->r1->final[i][j][l] = f->r1->qdiffs[i][l] + j + 2;
                    if(f->r1->final[i][j][l] < 2) f->r1->final[i][j][l] = 2;
                    f->r2->final[i][j][l] = f->r2->qdiffs[i][l] + j + 2;
                    if(f->r2->final[i][j][l] < 2) f->r2->final[i][j][l] = 2;
                }
            }
        }
    }

    void fill_qvals(fullerr_t *f)
    {
        int i;
        uint64_t l;
        for(i = 0; i < 4; ++i) {
            for(l = 0; l < f->l; ++l) {
                for(unsigned j = 1; j < NQSCORES; ++j) { // Skip qualities of 2
                    f->r1->qpvsum[i][l] += std::pow(10., (double)(-0.1 * (j + 2))) * f->r1->obs[i][j][l];
                    f->r1->qobs[i][l] += f->r1->obs[i][j][l];
                    f->r1->qerr[i][l] += f->r1->err[i][j][l];
                    f->r2->qpvsum[i][l] += std::pow(10., (double)(-0.1 * (j + 2))) * f->r2->obs[i][j][l];
                    f->r2->qobs[i][l] += f->r2->obs[i][j][l];
                    f->r2->qerr[i][l] += f->r2->err[i][j][l];
                }
            }
        }
        for(i = 0; i < 4; ++i) {
            for(l = 0; l < f->l; ++l) {
                f->r1->qpvsum[i][l] /= f->r1->qobs[i][l]; // Get average ILMN-reported quality
                f->r2->qpvsum[i][l] /= f->r2->qobs[i][l]; // Divide by observations of cycle/base call
                //LOG_DEBUG("qpvsum: %f. Mean p value: %f.  pv2ph mean p value: %i.\n", f->r1->qpvsum[i][l] * f->r1->qobs[i][l], f->r1->qpvsum[i][l], pv2ph(f->r1->qpvsum[i][l]));
                f->r1->qdiffs[i][l] = (f->r1->qobs[i][l] >= f->min_obs) ? pv2ph((double)f->r1->qerr[i][l] / f->r1->qobs[i][l]) - pv2ph(f->r1->qpvsum[i][l])
                                                                     : 0;
                f->r2->qdiffs[i][l] = (f->r2->qobs[i][l] >= f->min_obs) ? pv2ph((double)f->r2->qerr[i][l] / f->r2->qobs[i][l]) - pv2ph(f->r2->qpvsum[i][l])
                                                                     : 0;
            }
        }
    }

    void fill_sufficient_obs(fullerr_t *f)
    {
        for(int i = 0; i < 4; ++i) {
            for(unsigned j = 0; j < NQSCORES; ++j) {
                for(uint64_t l = 0; l < f->l; ++l) {
                    if(f->r1->obs[i][j][l] >= f->min_obs)
                        f->r1->final[i][j][l] = pv2ph((double)f->r1->err[i][j][l] / f->r1->obs[i][j][l]);
                    if(f->r2->obs[i][j][l] >= f->min_obs)
                        f->r2->final[i][j][l] = pv2ph((double)f->r2->err[i][j][l] / f->r2->obs[i][j][l]);
                }
            }
        }
    }


    void write_counts(fullerr_t *f, FILE *cp, FILE *ep)
    {
        FILE *dictwrite = fopen("dict.txt", "w");
        fprintf(dictwrite, "{\n\t");
        unsigned i, j, l;
        for(l = 0; l < f->l; ++l) {
            for(j = 0; j < NQSCORES; ++j) {
                for(i = 0; i < 4u; ++i) {
                    fprintf(dictwrite, "'r1,%c,%i,%u,obs': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r1->obs[i][j][l]);
                    fprintf(dictwrite, "'r2,%c,%i,%u,obs': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r2->obs[i][j][l]);
                    fprintf(dictwrite, "'r1,%c,%i,%u,err': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r1->err[i][j][l]);
                    if(i == 3 && j == NQSCORES - 1 && l == f->l - 1)
                        fprintf(dictwrite, "'r2,%c,%i,%u,err': %lu\n}", NUM2NUC_STR[i], j + 2, l + 1, f->r2->err[i][j][l]);
                    else
                        fprintf(dictwrite, "'r2,%c,%i,%u,err': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r2->err[i][j][l]);
                    fprintf(cp, i ? ":%lu": "%lu", f->r1->obs[i][j][l]);
                    fprintf(ep, i ? ":%lu": "%lu", f->r1->err[i][j][l]);
                }
                if(j != NQSCORES - 1) {
                    fprintf(ep, ","); fprintf(cp, ",");
                }
            }
            fprintf(ep, "|"); fprintf(cp, "|");
            for(j = 0; j < NQSCORES; ++j) {
                for(i = 0; i < 4; ++i) {
                    fprintf(cp, i ? ":%lu": "%lu", f->r2->obs[i][j][l]);
                    fprintf(ep, i ? ":%lu": "%lu", f->r2->err[i][j][l]);
                }
                if(j != NQSCORES - 1) {
                    fprintf(ep, ","); fprintf(cp, ",");
                }
            }
            fprintf(ep, "\n"); fprintf(cp, "\n");
        }
        fclose(dictwrite);
    }

    void write_3d_offsets(FILE *fp, fullerr_t *f)
    {
        fprintf(fp, "#Cycle\tR1A\tR1C\tR1G\tR1T\tR2A\tR2C\tR2G\tR2T\n");
        for(uint64_t l = 0; l < f->l; ++l) {
            fprintf(fp, "%lu\t", l + 1);
            int i;
            for(i = 0; i < 4; ++i) fprintf(fp, i ? "\t%i": "%i", f->r1->qdiffs[i][l]);
            fputc('|', fp);
            for(i = 0; i < 4; ++i) fprintf(fp, i ? "\t%i": "%i", f->r2->qdiffs[i][l]);
            fputc('\n', fp);
        }
        return;
    }

    readerr_t *readerr_init(size_t l) {
        readerr_t *ret = (readerr_t *)calloc(1, sizeof(readerr_t));
        arr3d_init(ret->obs, l, uint64_t);
        arr3d_init(ret->err, l, uint64_t);
        arr3d_init(ret->final, l, int);
        arr2d_init(ret->qdiffs, l, int);
        arr2d_init(ret->qpvsum, l, double);
        arr2d_init(ret->qobs, l, uint64_t);
        arr2d_init(ret->qerr, l, uint64_t);
        ret->l = l;
        return ret;
    }

    fullerr_t *fullerr_init(size_t l, char *bedpath, bam_hdr_t *hdr,
                            int padding, int minFM, int maxFM, int flag,
                            int minMQ, uint32_t minPV, uint64_t min_obs) {
        fullerr_t *ret = (fullerr_t *)calloc(1, sizeof(fullerr_t));
        ret->l = l;
        ret->r1 = readerr_init(l);
        ret->r2 = readerr_init(l);
        if(bedpath)
            ret->bed = dlib::parse_bed_hash(bedpath, hdr, padding);
        ret->minFM = minFM;
        ret->maxFM = maxFM;
        ret->flag = flag;
        ret->minMQ = minMQ;
        ret->minPV = minPV;
        ret->min_obs = min_obs;
        return ret;
    }

    void fullerr_destroy(fullerr_t *e) {
        if(e->r1) readerr_destroy(e->r1), e->r1 = nullptr;
        if(e->r2) readerr_destroy(e->r2), e->r2 = nullptr;
        cond_free(e->refcontig);
        if(e->bed) {
            kh_destroy(bed, e->bed);
        }
        free(e);
    }

    fmerr_t *fm_init(char *bedpath, bam_hdr_t *hdr, const char *refcontig, int padding, int flag, int minMQ, uint32_t minPV, double minFR) {
        fmerr_t *ret = (fmerr_t *)calloc(1, sizeof(fmerr_t));
        if(bedpath && *bedpath) {
            ret->bed = dlib::parse_bed_hash(bedpath, hdr, padding);
            ret->bedpath = strdup(bedpath);
        }
        if(refcontig && *refcontig) {
            ret->refcontig = strdup(refcontig);
        }
        ret->hash1 = kh_init(obs);
        ret->hash2 = kh_init(obs);
        ret->flag = flag;
        ret->minMQ = minMQ;
        ret->minPV = minPV;
        ret->minFR = minFR;
        return ret;
    }

    void fm_destroy(fmerr_t *fm) {
        if(fm->bed) kh_destroy(bed, fm->bed);
        kh_destroy(obs, fm->hash1);
        kh_destroy(obs, fm->hash2);
        cond_free(fm->refcontig);
        cond_free(fm->bedpath);
        free(fm);
    }

    cycle_err_t *cycle_init(char *bedpath, bam_hdr_t *hdr, char *refcontig, int padding, int minMQ, int rlen, int flag)
    {
        cycle_err_t *ret = (cycle_err_t *)calloc(1, sizeof(cycle_err_t));
        if(bedpath && *bedpath) {
            ret->bed = dlib::parse_bed_hash(bedpath, hdr, padding);
            ret->bedpath = strdup(bedpath);
        }
        if(refcontig && *refcontig) {
            ret->refcontig = strdup(refcontig);
        }
        ret->rlen = rlen;
        ret->r1 = (obserr_t *)calloc(rlen, sizeof(obserr_t));
        ret->r2 = (obserr_t *)calloc(rlen, sizeof(obserr_t));
        ret->minMQ = minMQ;
        ret->flag = flag;
        return ret;
    }

    void cycle_destroy(cycle_err_t *c)
    {
        if(c->bed) kh_destroy(bed, c->bed), c->bed = nullptr;
        cond_free(c->bedpath);
        cond_free(c->refcontig);
        cond_free(c->r1); cond_free(c->r2);
        free(c);
    }


    int err_usage(int exit_status)
    {
        fprintf(stderr,
                        "bmftools err\nSubcommands:\n"
                        "\tmain:\n"
                        "\t\tCalculates error rates by family size, cycle, base call, and quality score.\n"
                        "\t\tProvides a variety of ways of slicing the data.\n"
                        "\tfm:\n"
                        "\t\tCalculates error rates by family size.\n"
                        "\tcycle:\n"
                        "\t\tCalculates error rates by cycle.\n"
                        "\tregion:\n"
                        "\t\tCalculates error rates by bed region.\n"
                );
        exit(exit_status);
        return exit_status; // This never happens
    }

    int err_main(int argc, char *argv[])
    {
        if(argc < 2) return err_usage(EXIT_FAILURE);
        if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
            return err_usage(EXIT_SUCCESS);
        if(strcmp(argv[1], "main") == 0)
            return err_main_main(argc - 1, argv + 1);
        if(strcmp(argv[1], "fm") == 0)
            return err_fm_main(argc - 1, argv + 1);
        if(strcmp(argv[1], "cycle") == 0)
            return err_cycle_main(argc - 1, argv + 1);
        if(strcmp(argv[1], "region") == 0)
            return err_region_main(argc - 1, argv + 1);
        LOG_EXIT("Unrecognized subcommand '%s'. Abort!\n", argv[1]);
        return EXIT_FAILURE;
    }


    int err_main_main(int argc, char *argv[])
    {
        htsFormat open_fmt = {sequence_data, bam, {1, 3}, gzip, 0, nullptr};
        samFile *fp = nullptr;
        bam_hdr_t *header = nullptr;
        int c, minMQ = 0;
        std::string outpath("");
        if(argc < 2) return err_main_usage(EXIT_FAILURE);

        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) err_main_usage(EXIT_SUCCESS);


        FILE *ofp = nullptr, *d3 = nullptr, *df = nullptr, *dbc = nullptr, *dc = nullptr, *global_fp = nullptr;
        char refcontig[200] = "";
        char *bedpath = nullptr;
        int padding = -1;
        int minFM = 0;
        int maxFM = INT_MAX;
        int flag = 0;
        uint32_t minPV = 0;
        uint64_t min_obs = default_min_obs;
        while ((c = getopt(argc, argv, "a:p:b:r:c:n:f:3:o:g:m:M:S:O:h?FdDP")) >= 0) {
            switch (c) {
            case 'a': minMQ = atoi(optarg); break;
            case 'd': flag |= REQUIRE_DUPLEX; break;
            case 'D': flag |= REFUSE_DUPLEX; break;
            case 'P': flag |= REQUIRE_PROPER; break;
            case 'F': flag |= REQUIRE_FP_PASS; break;
            case 'm': minFM = atoi(optarg); break;
            case 'M': maxFM = atoi(optarg); break;
            case 'f': df = dlib::open_ofp(optarg); break;
            case 'o': outpath = optarg; break;
            case 'O': min_obs = strtoull(optarg, nullptr, 10); break;
            case '3': d3 = dlib::open_ofp(optarg); break;
            case 'c': dc = dlib::open_ofp(optarg); break;
            case 'n': dbc = dlib::open_ofp(optarg); break;
            case 'r': strcpy(refcontig, optarg); break;
            case 'b': bedpath = strdup(optarg); break;
            case 'p': padding = atoi(optarg); break;
            case 'g': global_fp = dlib::open_ofp(optarg); break;
            case 'S': minPV = strtoul(optarg, nullptr, 0); break;
            case '?': case 'h': return err_main_usage(EXIT_SUCCESS);
            }
        }

        if(padding < 0 && bedpath && *bedpath)
            LOG_INFO((char *)"Padding not set. Setting to default value %i.\n", DEFAULT_PADDING);

        if (argc != optind+2)
            return err_main_usage(EXIT_FAILURE);

        faidx_t *fai = fai_load(argv[optind]);

        if ((fp = sam_open_format(argv[optind + 1], "r", &open_fmt)) == nullptr)
            LOG_EXIT("Cannot open input file \"%s\"", argv[optind]);
        if ((header = sam_hdr_read(fp)) == nullptr)
            LOG_EXIT("Failed to read header for \"%s\"", argv[optind]);

        if(minPV) dlib::check_bam_tag_exit(argv[optind + 1], "PV");
        if(minFM || maxFM != INT_MAX) dlib::check_bam_tag_exit(argv[optind + 1], "FM");

        // Get read length from the first read
        bam1_t *b = bam_init1();
        c = sam_read1(fp, header, b);
        fullerr_t *f = fullerr_init((size_t)b->core.l_qseq, bedpath, header,
                                     padding, minFM, maxFM, flag, minMQ, minPV, min_obs);
        sam_close(fp);
        fp = nullptr;
        bam_destroy1(b);
        if(*refcontig) f->refcontig = strdup(refcontig);
        bam_hdr_destroy(header), header = nullptr;
        err_main_core(argv[optind + 1], fai, f, &open_fmt);
        fai_destroy(fai);
        fill_qvals(f);
        impute_scores(f);
        //fill_sufficient_obs(f); Try avoiding the fill sufficients and only use observations.
        if(outpath.size()) {
            ofp = fopen(outpath.c_str(), "w");
            write_final(ofp, f);
            fclose(ofp);
        }

        if(d3) {
            write_3d_offsets(d3, f);
            fclose(d3), d3 = nullptr;
        }
        if(df) {
            write_full_rates(df, f);
            fclose(df), df = nullptr;
        }
        if(dbc) {
            write_base_rates(dbc, f);
            fclose(dbc), dbc = nullptr;
        }
        if(dc) {
            write_cycle_rates(dc, f);
            fclose(dc), dc = nullptr;
        }
        if(!global_fp) {
            LOG_INFO("No global rate outfile provided. Defaulting to stdout.\n");
            global_fp = stdout;
        }
        write_global_rates(global_fp, f); fclose(global_fp);
        fullerr_destroy(f);
        LOG_INFO("Successfully completed bmftools err main!\n");
        return EXIT_SUCCESS;
    }

    int err_cycle_main(int argc, char *argv[])
    {
        htsFormat open_fmt = {sequence_data, bam, {1, 3}, gzip, 0, nullptr};

        if(argc < 2) return err_cycle_usage(EXIT_FAILURE);

        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
            return err_cycle_usage(EXIT_SUCCESS);

        FILE *ofp = nullptr;
        int padding = -1, minMQ = 0, flag = 0, c;
        char *bedpath = nullptr, *outpath = nullptr, *refcontig = nullptr;
        while ((c = getopt(argc, argv, "p:b:r:o:a:h?P")) >= 0) {
            switch (c) {
            case 'a': minMQ = atoi(optarg); break;
            case 'o': outpath = optarg; break;
            case 'r': refcontig = strdup(optarg); break;
            case 'b': bedpath = strdup(optarg); break;
            case 'p': padding = atoi(optarg); break;
            case 'P': flag |= REQUIRE_PROPER; break;
            case '?': case 'h': return err_cycle_usage(EXIT_SUCCESS);
            }
        }

        if(padding < 0 && bedpath) {
            LOG_INFO("Padding not set. Setting to default value %i.\n", DEFAULT_PADDING);
        }

        if(!outpath) {
            LOG_WARNING("Output path not set. Defaulting to stdout.\n");
        }
        ofp = dlib::open_ofp(outpath);

        if (argc != optind+2)
            return err_cycle_usage(EXIT_FAILURE);

        faidx_t *fai = fai_load(argv[optind]);

        cycle_err_t *ce = err_cycle_core(argv[optind + 1], fai, &open_fmt, bedpath, refcontig, padding, minMQ, flag);
        err_cycle_report(ofp, ce); fclose(ofp);
        cycle_destroy(ce);
        fai_destroy(fai);
        cond_free(refcontig);
        cond_free(bedpath);
        LOG_INFO("Successfully completed bmftools err cycle!\n");
        return EXIT_SUCCESS;
    }


    int err_fm_main(int argc, char *argv[])
    {
        htsFormat open_fmt;
        memset(&open_fmt, 0, sizeof(htsFormat));
        open_fmt.category = sequence_data;
        open_fmt.format = bam;
        open_fmt.version.major = 1;
        open_fmt.version.minor = 3;
        samFile *fp = nullptr;
        bam_hdr_t *header = nullptr;
        char *outpath = nullptr;

        if(argc < 2) return err_fm_usage(EXIT_FAILURE);

        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) return err_fm_usage(EXIT_SUCCESS);



        FILE *ofp = nullptr;
        std::string refcontig("");
        char *bedpath = nullptr;
        int flag = 0, padding = -1, minMQ = 0, c;
        uint32_t minPV = 0;
        double minFR = 0.;
        while ((c = getopt(argc, argv, "S:p:b:r:o:a:f:Fh?dP")) >= 0) {
            switch (c) {
            case 'a': minMQ = atoi(optarg); break;
            case 'd': flag |= REQUIRE_DUPLEX; break;
            case 'o': outpath = optarg; break;
            case 'r': refcontig = optarg; break;
            case 'b': bedpath = strdup(optarg); break;
            case 'p': padding = atoi(optarg); break;
            case 'P': flag |= REQUIRE_PROPER; break;
            case 'S': minPV = strtoul(optarg, nullptr, 0); break;
            case 'F': flag |= REQUIRE_FP_PASS; break;
            case 'f':
                minFR = atof(optarg);
                if(minFR < 0.0 || minFR > 1.0) LOG_EXIT("minFR must be between 0 and 1. Given: %f.\n", minFR);
                break;
            case '?': case 'h': return err_fm_usage(EXIT_SUCCESS);
            }
        }

        if(bedpath && padding < 0) {
            LOG_INFO("Padding not set. Setting to default value %i.\n", DEFAULT_PADDING);
        }

        if(!outpath) {
            LOG_WARNING("Output path unset. Defaulting to stdout.\n");
        }
        ofp = dlib::open_ofp(outpath);

        if (argc != optind+2)
            return err_fm_usage(EXIT_FAILURE);

        faidx_t *fai = fai_load(argv[optind]);

        if ((fp = sam_open_format(argv[optind + 1], "r", &open_fmt)) == nullptr) {
            LOG_EXIT("Cannot open input file \"%s\"", argv[optind]);
        }
        if ((header = sam_hdr_read(fp)) == nullptr) {
            LOG_EXIT("Failed to read header for \"%s\"", argv[optind]);
        }
        for(auto tag: {"FM", "FP", "RV"})
            dlib::check_bam_tag_exit(argv[optind + 1], tag);
        if(flag & (REQUIRE_DUPLEX | REFUSE_DUPLEX))
            dlib::check_bam_tag_exit(argv[optind + 1], "DR");

        fmerr_t *f = fm_init(bedpath, header, refcontig.c_str(), padding, flag, minMQ, minPV, minFR);
        // Get read length from the first
        bam_hdr_destroy(header); header = nullptr;
        err_fm_core(argv[optind + 1], fai, f, &open_fmt);
        err_fm_report(ofp, f); fclose(ofp);
        fai_destroy(fai);
        fm_destroy(f);
        LOG_INFO("Successfully completed bmftools err fm!\n");
        return EXIT_SUCCESS;
    }

    static int read_bam(RegionExpedition *navy, bam1_t *b)
    {
        int ret;
        uint8_t *fmdata, *fpdata;
        for(;;)
        {

            if((ret = sam_itr_next(navy->fp, navy->iter, b)) < 0) break;
            if((b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) ||
                    (int)b->core.qual < navy->minMQ) continue;
            fmdata = bam_aux_get(b, "FM");
            fpdata = bam_aux_get(b, "FP");
            if ((fmdata && bam_aux2i(fmdata) < navy->minFM) ||
                (navy->requireFP && fpdata && bam_aux2i(fpdata) == 0))
                    continue;
            break;
        }
        return ret;
    }

    inline void region_loop(RegionErr& counter, char *ref, bam1_t *b)
    {
        int i, rc, fc, length, ind, s;
        uint32_t *const cigar = bam_get_cigar(b);
        uint8_t *seq = bam_get_seq(b);
        for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
            length = bam_cigar_oplen(cigar[i]);
            switch(bam_cigar_op(cigar[i])) {
            case BAM_CMATCH: case BAM_CEQUAL: case BAM_CDIFF:
                for(ind = 0; ind < length; ++ind) {
                    s = bam_seqi(seq, ind + rc);
                    if(s == dlib::htseq::HTS_N || ref[b->core.pos + fc + ind] == 'N') continue;
                    counter.inc_obs();
                    if(seq_nt16_table[(int8_t)ref[b->core.pos + fc + ind]] != s) counter.inc_err();
                }
                rc += length; fc += length;
                break;
            case BAM_CSOFT_CLIP: case BAM_CHARD_CLIP: case BAM_CINS:
                rc += length;
                break;
            case BAM_CREF_SKIP: case BAM_CDEL:
                fc += length;
                break;
            }
        }
        // Add actual error code!
    }

    void err_region_core(RegionExpedition *Holloway) {
        // Make region_counts classes. These can now be filled from the bam.
        char *ref = nullptr;
        int len;
        bam1_t *b = bam_init1();
        std::vector<khiter_t> sorted_keys(dlib::make_sorted_keys(Holloway->bed));
        for(khiter_t& k: sorted_keys) {
            if(!kh_exist(Holloway->bed, k)) continue;
            if(ref) free(ref);
            LOG_DEBUG("Fetching ref sequence for contig id %i.\n", kh_key(Holloway->bed, k));
            ref = fai_fetch(Holloway->fai, Holloway->hdr->target_name[kh_key(Holloway->bed, k)], &len);
            LOG_DEBUG("Fetched! Length: %i\n", len);
            for(unsigned i = 0; i < kh_val(Holloway->bed, k).n; ++i) {
                Holloway->region_counts.emplace_back(kh_val(Holloway->bed, k), i);
                LOG_DEBUG("Add new RegionErr. Size of vector: %lu.\n", Holloway->region_counts.size());
                const int start = get_start(kh_val(Holloway->bed, k).intervals[i]);
                const int stop = get_stop(kh_val(Holloway->bed, k).intervals[i]);
                LOG_DEBUG("Iterating through bam region with start %i and stop %i.\n", start, stop);
                if(Holloway->iter) hts_itr_destroy(Holloway->iter);
                Holloway->iter = sam_itr_queryi(Holloway->bam_index, kh_key(Holloway->bed, k),
                                                start, stop);
                while(read_bam(Holloway, b) >= 0) {
                    assert((unsigned)b->core.tid == kh_key(Holloway->bed, k));
                    if(bam_getend(b) <= start) {
                        LOG_DEBUG("Pos (%i) less than region start (%i).\n", b->core.pos, start);
                        continue;
                    } else if(b->core.pos >= stop) {
                        LOG_DEBUG("Pos (%i) g/e region stop (%i).\n", b->core.pos, stop);
                        break;
                    }
                    region_loop(Holloway->region_counts[Holloway->region_counts.size() - 1], ref, b);
                }
            }
        }
        bam_destroy1(b);
    }

    void write_region_rates(FILE *fp, RegionExpedition& Holloway) {
        fprintf(fp, "#Region name\t%%Error Rate\t#Errors\t#Obs\n");
        for(RegionErr& re: Holloway.region_counts)
            re.write_report(fp);
    }

    int err_region_main(int argc, char *argv[])
    {
        if(argc < 2) return err_region_usage(EXIT_FAILURE);

        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
            return err_region_usage(EXIT_SUCCESS);

        FILE *ofp = nullptr;
        int padding = -1, minMQ = 0, minFM = 0, c, requireFP = 0;
        char *bedpath = nullptr, *outpath = nullptr;
        faidx_t *fai;
        while ((c = getopt(argc, argv, "p:b:r:o:a:h?q")) >= 0) {
            switch (c) {
            case 'q': requireFP = 1; break;
            case 'a': minMQ = atoi(optarg); break;
            case 'f': minFM = atoi(optarg); break;
            case 'o': outpath = strdup(optarg); break;
            case 'b': bedpath = strdup(optarg); break;
            case 'p': padding = atoi(optarg); break;
            case '?': case 'h': return err_region_usage(EXIT_SUCCESS);
            }
        }

        if(padding < 0 && bedpath) {
            LOG_INFO("Padding not set. Setting to default value %i.\n", DEFAULT_PADDING);
        }

        if(!outpath) {
            outpath = strdup("-");
            LOG_WARNING("Output path not set. Defaulting to stdout.\n");
        }
        if(!bedpath) {
            LOG_EXIT("Bed file required for bmftools err region.\n");
        }
        ofp = dlib::open_ofp(outpath);

        if (argc != optind+2)
            return err_region_usage(EXIT_FAILURE);

        fai = fai_load(argv[optind]);
        if(!fai) {
            LOG_EXIT("Could not load fasta index for %s. Abort!\n", argv[optind]);
        }

        RegionExpedition Holloway(argv[optind + 1], bedpath, fai, minMQ, padding, minFM, requireFP);
        err_region_core(&Holloway);
        write_region_rates(ofp, Holloway), fclose(ofp);
        fai_destroy(fai);
        cond_free(bedpath);
        cond_free(outpath);
        LOG_INFO("Successfully completed bmftools err region!\n");
        return EXIT_SUCCESS;
    }

}
