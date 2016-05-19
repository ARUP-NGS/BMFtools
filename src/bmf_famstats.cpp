#include <getopt.h>
#include <algorithm>
#include "dlib/bam_util.h"

namespace BMF {

    const char *tags_to_check[] = {"FP", "FM", "FA"};

    int famstats_frac_usage(int exit_status) {
        fprintf(stderr,
                        "Calculates the fraction of raw reads with family size >= parameter.\n"
                        "Usage: bmftools famstats frac <opts> <minFM> <in.bam>\n"
                        "-n: Set notification interval. Default: 1000000.\n"
                        "-h, -?: Return usage.\n"
                );
        exit(exit_status);
        return exit_status; // This never happens
    }


    KHASH_MAP_INIT_INT64(fm, uint64_t)


    struct famstats_t {
        uint64_t n_pass;
        uint64_t fp_fail_count;
        uint64_t n_fm_fail;
        uint64_t n_flag_fail;
        uint64_t allfm_sum;
        uint64_t allfm_counts;
        uint64_t allrc_sum;
        uint64_t realfm_sum;
        uint64_t realfm_counts;
        uint64_t realrc_sum;
        uint64_t dr_sum;
        uint64_t dr_counts;
        uint64_t dr_rc_sum;
        double dr_rc_frac_sum;
        khash_t(fm) *fm;
        khash_t(fm) *np;
        khash_t(fm) *rc;
        khiter_t ki;
        uint8_t *data;
    };

    struct famstats_fm_settings_t {
        uint64_t notification_interval;
        uint32_t minmq:8;
        uint32_t skip_fp_fail:1;
        uint32_t minFM:16;
    };

    struct fm_t{
        uint64_t n; // Number of times observed
        uint64_t fm; // Number of times observed
    };


    int get_nbins(khash_t(fm) *table)
    {
        int ret = 0;
        for(khiter_t k = kh_begin(table); k != kh_end(table); ++k)
            if(kh_exist(table, k))
                ++ret;
        return ret;
    }


    static void print_hashstats(famstats_t *stats, FILE *fp)
    {
        std::vector<fm_t> fms(stats->fm->n_occupied);
        unsigned i;
        khiter_t ki;
        fprintf(fp, "#Family size\tNumber of families\n");
        for(i = 0, ki = kh_begin(stats->fm); ki != kh_end(stats->fm); ++ki)
            if(kh_exist(stats->fm, ki))
                fms[i++] = {kh_val(stats->fm, ki), kh_key(stats->fm, ki)};
        std::sort(fms.begin(), fms.end(), [](const fm_t a, const fm_t b){
            return a.fm < b.fm;
        });
        for(i = 0; i < stats->fm->n_occupied; ++i)
            fprintf(fp, "%lu\t%lu\n", fms[i].fm, fms[i].n);

        fms.resize(stats->rc->n_occupied);
        for(i = 0, ki = kh_begin(stats->rc); ki != kh_end(stats->rc); ++ki)
            if(kh_exist(stats->rc, ki))
                fms[i++] = {kh_val(stats->rc, ki), kh_key(stats->rc, ki)};
        std::sort(fms.begin(), fms.end(), [](const fm_t a, const fm_t b){
            return a.fm < b.fm;
        });
        if(fms[0].fm != (uint64_t)-1) {
            fprintf(fp, "#RV'd in family\tNumber of families\n");
            for(i = 0; i < stats->rc->n_occupied; ++i)
                fprintf(fp, "%lu\t%lu\n", fms[i].fm, fms[i].n);
        }
        // Handle stats->np
        fms.resize(stats->np->n_occupied);
        size_t n_rsq_fams = 0;
        for(i = 0, ki = kh_begin(stats->np); ki != kh_end(stats->np); ++ki) {
            if(kh_exist(stats->np, ki)) {
                fms[i++] = {kh_val(stats->np, ki), kh_key(stats->np, ki)};
                n_rsq_fams += kh_val(stats->np, ki);
            }
        }
        std::sort(fms.begin(), fms.end(), [](const fm_t a, const fm_t b){
            return a.fm < b.fm;
        });
        fprintf(fp, "#Number of families that were rescued: %lu\n", n_rsq_fams);
        fputs("#Number of pre-rescue reads in rescued\tNumber of families\n", fp);
        for(i = 0; i < stats->np->n_occupied; ++i)
            fprintf(fp, "%lu\t%lu\n", fms[i].fm, fms[i].n);
    }


    static void print_stats(famstats_t *stats, FILE *fp, famstats_fm_settings_t *settings)
    {
        fprintf(fp, "#Number passing filters: %lu\n", stats->n_pass);
        fprintf(fp, "#Number failing filters: %lu\n", stats->fp_fail_count + stats->n_fm_fail + stats->n_flag_fail);
        if(settings->skip_fp_fail)
            fprintf(fp, "#Number failing FP filters: %lu\n", stats->fp_fail_count);
        else
            fprintf(fp, "#Count for FP failed reads, still included in total counts: %lu\n", stats->fp_fail_count);
        fprintf(fp, "#Number failing FM filters: %lu\n", stats->n_fm_fail);
        fprintf(fp, "#Number failing flag filters (secondary, supplementary): %lu\n", stats->n_flag_fail);
        fprintf(fp, "#Summed FM (total founding reads): %lu\n", stats->allfm_sum);
        fprintf(fp, "#Summed FM (total founding reads), (FM > 1): %lu\n", stats->realfm_sum);
        fprintf(fp, "#Summed RV (total reverse-complemented reads): %lu\n", stats->allrc_sum);
        fprintf(fp, "#Summed RV (total reverse-complemented reads), (FM > 1): %lu\n", stats->realrc_sum);
        fprintf(fp, "#RV fraction for all read families: %f\n", (double)stats->allrc_sum / (double)stats->allfm_sum);
        fprintf(fp, "#RV fraction for real read families: %f\n", (double)stats->realrc_sum / (double)stats->realfm_sum);
        fprintf(fp, "#Mean Family Size (all)\t%f\n", (double)stats->allfm_sum / (double)stats->allfm_counts);
        fprintf(fp, "#Mean Family Size (real)\t%f\n", (double)stats->realfm_sum / (double)stats->realfm_counts);
        if(stats->allrc_sum) {
            fprintf(fp, "#Duplex fraction of unique observations\t%0.12f\n", (double)stats->dr_counts / stats->n_pass);
            fprintf(fp, "#Fraction of raw reads in duplex families\t%0.12f\n", (double)stats->dr_sum / stats->allfm_sum);
            fprintf(fp, "#Mean fraction of reverse reads within each duplex family\t%0.12f\n", stats->dr_rc_frac_sum / stats->dr_counts);
            fprintf(fp, "#Mean fraction of reverse reads within all duplex families\t%0.12f\n", (double)stats->dr_rc_sum / stats->dr_sum);
        }
        print_hashstats(stats, fp);
    }



    static inline void famstats_fm_loop(famstats_t *s, bam1_t *b, famstats_fm_settings_t *settings)
    {
        uint8_t *data;
        if(b->core.flag & BAM_FREAD2) return; // Silently skip all read 2s since they have the same FM values.
        if((b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ||
                b->core.qual < settings->minmq) {
            ++s->n_flag_fail;
            return;
        }
        const int FM = ((data = bam_aux_get(b, "FM")) != nullptr ? bam_aux2i(data) : 0);
        const int NP = ((data = bam_aux_get(b, "NP")) != nullptr ? bam_aux2i(data) : -1);
        int RV = ((data = bam_aux_get(b, "RV")) != nullptr ? bam_aux2i(data) : -1);
        if(UNLIKELY(FM == 0)) LOG_EXIT("Missing required FM tag. Abort!\n");
        if(FM < settings->minFM) {
            ++s->n_fm_fail;
            return;
        }
        if(bam_itag(b, "FP") == 0) {
            ++s->fp_fail_count;
            if(settings->skip_fp_fail) return;
        }
        ++s->n_pass;

        if(FM > 1) {
            ++s->realfm_counts;
            s->realfm_sum += FM;
            s->realrc_sum += RV < 0 ? 0 : RV;
        }
        ++s->allfm_counts;
        s->allfm_sum += FM;
        s->allrc_sum += RV < 0 ? 0 : RV;

        int khr;
        // Have we seen this family size before?
        if((s->ki = kh_get(fm, s->fm, FM)) == kh_end(s->fm))
            // If not, put it into the hash table with a count of 1.
            s->ki = kh_put(fm, s->fm, FM, &khr), kh_val(s->fm, s->ki) = 1;
        else ++kh_val(s->fm, s->ki); // Otherwise increment counts
        // Same, but for RV
        if((s->ki = kh_get(fm, s->rc, RV)) == kh_end(s->rc))
            s->ki = kh_put(fm, s->rc, RV, &khr), kh_val(s->rc, s->ki) = 1;
        else ++kh_val(s->rc, s->ki);
        // Same, but for NP
        if(NP > 0) {
            if((s->ki = kh_get(fm, s->np, NP)) == kh_end(s->np))
                s->ki = kh_put(fm, s->np, NP, &khr), kh_val(s->np, s->ki) = 1;
            else ++kh_val(s->np, s->ki);
        }

        // If the Duplex Read tag is present, increment duplex read counts
        uint8_t *dr_data = bam_aux_get(b, "DR");
        if(dr_data && bam_aux2i(dr_data)) {
            if(RV < 0) RV = 0;
            s->dr_sum += FM;
            ++s->dr_counts;
            s->dr_rc_sum += RV;
            s->dr_rc_frac_sum += (double)RV / FM;
        }
    }


    famstats_t *famstats_fm_core(dlib::BamHandle& handle, famstats_fm_settings_t *settings)
    {
        uint64_t count = 0;
        famstats_t *s = (famstats_t*)calloc(1, sizeof(famstats_t));
        int ret;
        s->fm = kh_init(fm);
        s->rc = kh_init(fm);
        s->np = kh_init(fm);
        s->data = nullptr;
        while (LIKELY((ret = handle.next()) >= 0)) {
            famstats_fm_loop(s, handle.rec, settings);
            if(UNLIKELY(++count % settings->notification_interval == 0))
                LOG_INFO("Number of records processed: %lu.\n", count);
        }
        if (ret != -1) LOG_WARNING("Truncated file? Continue anyway.\n");
        return s;
    }

    static int famstats_usage_exit(int exit_status)
    {
        fprintf(stderr,
                        "Calculates various utilities regarding family size on a given bam.\n"
                        "Usage: bmftools famstats <subcommand> <subcommand opts>\n"
                        "Subcommands: \nfm\tFamily Size stats\n"
                        "frac\tFraction of raw reads in family sizes >= minFM parameter.\n"
                );
        exit(exit_status);
        return exit_status;
    }

    static int famstats_fm_usage(int exit_status)
    {
        fprintf(stderr,
                        "Produces a histogram of family sizes and reverse counts.\n"
                        "Usage: bmftools famstats fm <opts> <in.bam>\n"
                        "Flags:\n"
                        "-m Set minimum mapping quality. Default: 0.\n"
                        "-f Set minimum family size. Default: 0.\n"
                );
        exit(exit_status);
        return exit_status;
    }

    int famstats_fm_main(int argc, char *argv[])
    {
        famstats_t *s;
        int c;


        famstats_fm_settings_t settings{0};
        settings.notification_interval = 1000000uL;

        while ((c = getopt(argc, argv, "m:f:n:Fh?")) >= 0) {
            switch (c) {
            case 'm':
                settings.minmq = atoi(optarg); break;
                break;
            case 'f':
                settings.minFM = atoi(optarg); break;
                break;
            case 'F':
                settings.skip_fp_fail = 1; break;
            case 'n': settings.notification_interval = strtoull(optarg, nullptr, 0); break;
            case '?': case 'h':
                return famstats_fm_usage(EXIT_SUCCESS);
            }
        }

        if (argc != optind+1) {
            return famstats_fm_usage((argc == optind) ? EXIT_SUCCESS: EXIT_FAILURE);
        }

        LOG_INFO("Running famstats fm with minmq %i and minFM %i.\n", settings.minmq, settings.minFM);
        for(const char *tag: tags_to_check)
            dlib::check_bam_tag_exit(argv[optind], tag);

        dlib::BamHandle handle(argv[optind]);
        s = famstats_fm_core(handle, &settings);
        print_stats(s, stdout, &settings);
        kh_destroy(fm, s->fm);
        kh_destroy(fm, s->np);
        kh_destroy(fm, s->rc);
        free(s);
        LOG_INFO("Successfully complete bmftools famstats fm.\n");
        return EXIT_SUCCESS;
    }

    int famstats_frac_main(int argc, char *argv[])
    {
        int c;
        uint64_t notification_interval = 1000000;

        if(argc < 2) famstats_frac_usage(EXIT_FAILURE);
        if(strcmp(argv[1], "--help") == 0) famstats_frac_usage(EXIT_SUCCESS);

        while ((c = getopt(argc, argv, "n:m:h?")) >= 0) {
            switch (c) {
            case 'n':
                notification_interval = strtoull(optarg, nullptr, 0); break;
            case '?': case 'h':
                return famstats_frac_usage(EXIT_SUCCESS);
            }
        }

        if (argc != optind+2) {
            if (argc == optind) famstats_frac_usage(EXIT_SUCCESS);
            else famstats_frac_usage(EXIT_FAILURE);
        }

        uint32_t minFM = strtoul(argv[optind], nullptr, 10);
        LOG_INFO("MinFM %i.\n", minFM);
        for(const char *tag: tags_to_check) dlib::check_bam_tag_exit(argv[optind+1], tag);
        dlib::BamHandle handle(argv[optind + 1]);
        uint64_t fm_above = 0, total_fm = 0, count = 0;
        // Check to see if the required tags are present before starting.
        int FM;
        int ret;
        while (LIKELY((ret = handle.next()) >= 0)) {
            // Filter reads
            if((handle.rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FREAD2)) ||
                    bam_itag(handle.rec, "FP") == 0)
                continue;
            FM = bam_itag(handle.rec, "FM");
            total_fm += FM;
            if((unsigned)FM >= minFM) fm_above += FM;
            if(UNLIKELY(!(++count % notification_interval)))
                LOG_INFO("Number of records processed: %lu.\n", count);
        }
        if (ret != -1) LOG_WARNING("Truncated file? Continue anyway.\n");
        fprintf(stdout, "#Fraction of raw reads with >= minFM %u:\t%f\n",
                minFM, (double)fm_above / total_fm);
        LOG_INFO("Successfully complete bmftools famstats frac.\n");
        return EXIT_SUCCESS;
    }

    int famstats_main(int argc, char *argv[])
    {
        if(argc < 2)
            return famstats_usage_exit(EXIT_FAILURE);
        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
            return famstats_usage_exit(EXIT_SUCCESS);
        if(strcmp(argv[1], "fm") == 0)
            return famstats_fm_main(argc - 1, argv + 1);
        if(strcmp(argv[1], "frac") == 0)
            return famstats_frac_main(argc - 1, argv + 1);
        fprintf(stderr, "[E:%s] Unrecognized subcommand '%s'. See usage.\n", __func__, argv[1]);
        return famstats_usage_exit(EXIT_FAILURE);
    }

}
