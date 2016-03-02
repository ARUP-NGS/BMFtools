#include "bmf_famstats.h"

const char *tags_to_check[] = {"FP", "RV", "FM", "FA"};
int RVWarn = 1;

int famstats_frac_usage_exit(FILE *fp, int code) {
    fprintf(fp, "bmftools famstats frac <opts> <in.bam>\n"
            "Opts:\n-m minFM to accept. REQUIRED.\n"
            "-h, -?: Return usage.\n");
    exit(code);
    return code; // This never happens
}

typedef struct {
    uint64_t n; // Number of times observed
    uint64_t fm; // Number of times observed
} fm_t;


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
    std::vector<fm_t> fms = std::vector<fm_t>(stats->fm->n_occupied);
    unsigned i;
    fprintf(fp, "#Family size\tNumber of families\n");
    for(i = 0, stats->ki = kh_begin(stats->fm); stats->ki != kh_end(stats->fm); ++stats->ki)
        if(kh_exist(stats->fm, stats->ki))
            fms[i++] = {kh_val(stats->fm, stats->ki), kh_key(stats->fm, stats->ki)};
    std::sort(fms.begin(), fms.end(), [](const fm_t a, const fm_t b){
        return a.fm < b.fm;
    });
    for(i = 0; i < stats->fm->n_occupied; ++i)
        fprintf(fp, "%lu\t%lu\n", fms[i].fm, fms[i].n);
    fprintf(fp, "#RV'd in family\tNumber of families\n");

    fms.resize(stats->rc->n_occupied);
    for(i = 0, stats->ki = kh_begin(stats->rc); stats->ki != kh_end(stats->rc); ++stats->ki)
        if(kh_exist(stats->rc, stats->ki))
            fms[i++] = {kh_val(stats->rc, stats->ki), kh_key(stats->rc, stats->ki)};
    std::sort(fms.begin(), fms.end(), [](const fm_t a, const fm_t b){
        return a.fm < b.fm;
    });
    for(i = 0; i < stats->rc->n_occupied; ++i)
        fprintf(fp, "%lu\t%lu\n", fms[i].fm, fms[i].n);
}


int famstats_target_usage_exit(FILE *fp, int success)
{
    fprintf(fp, "Usage: bmftools famstats target <opts> <in.bam>\nOpts:\n-b Path to bed file.\n"
            "-p padding. Number of bases around bed regions to pad. Default: %i.\n"
            "-h, -?: Return usage.\n", DEFAULT_PADDING);
    exit(success);
    return success;
}


int famstats_target_main(int argc, char *argv[])
{
    int c;
    khash_t(bed) *bed = NULL;
    char *bedpath = NULL;
    uint32_t padding = (uint32_t)-1;
    uint64_t notification_interval = 1000000;

    if(argc < 4) return famstats_target_usage_exit(stderr, EXIT_SUCCESS);

    if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
        return famstats_target_usage_exit(stderr, EXIT_SUCCESS);

    while ((c = getopt(argc, argv, "b:p:n:h?")) >= 0) {
        switch (c) {
        case 'b':
            bedpath = strdup(optarg);
            break;
        case 'p':
            padding = strtoul(optarg, NULL, 0);
            break;
        case 'n':
            notification_interval = strtoull(optarg, NULL, 0);
            break;
        case '?': case 'h':
            return famstats_target_usage_exit(stderr, EXIT_SUCCESS);
        }
    }


    if(padding == (uint32_t)-1) {
        padding = DEFAULT_PADDING;
        LOG_INFO("Padding not set. Set to default value (%u).\n", DEFAULT_PADDING);
    }

    if (argc != optind+1)
        return famstats_target_usage_exit(stderr, EXIT_FAILURE);

    if(!bedpath) {
        fprintf(stderr, "[E:%s] Bed path required for famstats target. See usage.\n", __func__);
        return famstats_target_usage_exit(stderr, EXIT_FAILURE);
    }

    dlib::BamHandle handle(argv[optind]);
    // If bedfile provided, use it. If not, calculate by contig. Why? Why not?
    if(!bedpath) {
        LOG_WARNING("Can't calculate on-target without a bed file. Abort!\n");
    }
    bed = parse_bed_hash(bedpath, handle.header, padding);
    uint64_t fm_target = 0, total_fm = 0, count = 0, n_flag_skipped = 0, n_fp_skipped = 0;
    uint8_t *fpdata = NULL;
    while(LIKELY(handle.next() >= 0)) {
        if((handle.rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FSUPPLEMENTARY))) {
            ++n_flag_skipped;
            continue;
        } else if((fpdata = bam_aux_get(handle.rec, "FP")) != NULL && !bam_aux2i(fpdata)) {
            ++n_fp_skipped;
            continue;
        }
        const int FM = bam_aux2i(bam_aux_get(handle.rec, "FM"));
        total_fm += FM;
        if(bed_test(handle.rec, bed))
            fm_target += FM;
        if(UNLIKELY(++count % notification_interval == 0))
            LOG_INFO("Number of records processed: %lu.\n", count);
    }
    LOG_INFO("#Number of records read: %lu. Number skipped (flag): %lu. Number skipped (FP): %lu.\n", count, n_flag_skipped, n_fp_skipped);
    bed_destroy_hash(bed);
    cond_free(bedpath);
    fprintf(stdout, "Fraction of raw reads on target: %f. \nTotal raw reads: %lu. Raw reads on target: %lu.\n",
            (double)fm_target / total_fm, total_fm, fm_target);
    return EXIT_SUCCESS;
}


static void print_stats(famstats_t *stats, FILE *fp)
{
    fprintf(fp, "#Number passing filters: %lu.\n", stats->n_pass);
    fprintf(fp, "#Number failing filters: %lu.\n", stats->n_fp_fail + stats->n_fm_fail + stats->n_flag_fail);
    fprintf(fp, "#Number failing FP filters: %lu.\n", stats->n_fp_fail);
    fprintf(fp, "#Number failing FM filters: %lu.\n", stats->n_fm_fail);
    fprintf(fp, "#Number failing flag filters: %lu.\n", stats->n_flag_fail);
    fprintf(fp, "#Summed FM (total founding reads): %lu.\n", stats->allfm_sum);
    fprintf(fp, "#Summed FM (total founding reads), (FM > 1): %lu.\n", stats->realfm_sum);
    fprintf(fp, "#Summed RV (total reverse-complemented reads): %lu.\n", stats->allrc_sum);
    fprintf(fp, "#Summed RV (total reverse-complemented reads), (FM > 1): %lu.\n", stats->realrc_sum);
    fprintf(fp, "#RV fraction for all read families: %lf.\n", (double)stats->allrc_sum / (double)stats->allfm_sum);
    fprintf(fp, "#RV fraction for real read families: %lf.\n", (double)stats->realrc_sum / (double)stats->realfm_sum);
    fprintf(fp, "#Mean Family Size (all)\t%lf\n", (double)stats->allfm_sum / (double)stats->allfm_counts);
    fprintf(fp, "#Mean Family Size (real)\t%lf\n", (double)stats->realfm_sum / (double)stats->realfm_counts);
    if(stats->dr_counts) {
        fprintf(fp, "#Duplex fraction of unique observations\t%0.12lf\n", (double)stats->dr_counts / stats->n_pass);
        fprintf(fp, "#Fraction of raw reads in duplex families\t%0.12lf\n", (double)stats->dr_sum / stats->allfm_sum);
        fprintf(fp, "#Mean fraction of reverse reads within each duplex family\t%0.12lf\n", stats->dr_rc_frac_sum / stats->dr_rc_sum);
        fprintf(fp, "#Mean fraction of reverse reads within all duplex families\t%0.12lf\n", (double)stats->dr_rc_sum / stats->dr_sum);
    }
    print_hashstats(stats, fp);
}



static inline void famstats_fm_loop(famstats_t *s, bam1_t *b, famstat_settings_t *settings)
{
    //Since R1 and R2
    if((b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FREAD2)) ||
            b->core.qual < settings->minMQ) {
        ++s->n_flag_fail;
        return;
    }
    const int FM = bam_itag(b, "FM");
    const int RV = bam_itag(b, "RV");
    if(FM < settings->minFM) {
        ++s->n_fm_fail;
        return;
    }
    if(!bam_itag(b, "FP")) {
        ++s->n_fp_fail;
        return;
    }
    ++s->n_pass;

    if(FM > 1) {
        ++s->realfm_counts;
        s->realfm_sum += FM;
        s->realrc_sum += RV;
    }
    ++s->allfm_counts;
    s->allfm_sum += FM;
    s->allrc_sum += RV;

    // Have we seen this family size before?
    if((s->ki = kh_get(fm, s->fm, FM)) == kh_end(s->fm))
        // If not, put it into the hash table with a count of 1.
        s->ki = kh_put(fm, s->fm, FM, &s->khr), kh_val(s->fm, s->ki) = 1;
    else ++kh_val(s->fm, s->ki); // Otherwise increment counts
    // Same, but for RV
    if((s->ki = kh_get(rc, s->rc, RV)) == kh_end(s->rc))
        s->ki = kh_put(rc, s->rc, RV, &s->khr), kh_val(s->rc, s->ki) = 1;
    else ++kh_val(s->rc, s->ki);

    // If the Duplex Read tag is present, incrememnt duplex read counts
    uint8_t *dr_data = bam_aux_get(b, "DR");
    if(dr_data && bam_aux2i(dr_data)) {
        s->dr_sum += FM;
        ++s->dr_counts;
        s->dr_rc_sum += RV;
        s->dr_rc_frac_sum += (double)RV / FM;
    }
}


famstats_t *famstats_fm_core(dlib::BamHandle& handle, famstat_settings_t *settings)
{
    uint64_t count = 0;
    famstats_t *s = (famstats_t*)calloc(1, sizeof(famstats_t));
    int ret;
    s->fm = kh_init(fm);
    s->rc = kh_init(rc);
    s->data = NULL;
    while (LIKELY((ret = handle.next()) >= 0)) {
        famstats_fm_loop(s, handle.rec, settings);
        if(UNLIKELY(++count % settings->notification_interval == 0))
            LOG_INFO("Number of records processed: %lu.\n", count);
    }
    if (ret != -1) LOG_WARNING("Truncated file? Continue anyway.\n");
    return s;
}

static int famstats_usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: bmftools famstats\n");
    fprintf(fp, "Subcommands: \nfm\tFamily Size stats\n"
            "frac\tFraction of raw reads in family sizes >= minFM parameter.\n"
            "target\tFraction of raw reads on target.\n");
    exit(exit_status);
    return exit_status;
}

static int famstats_fm_usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: bmftools famstats fm <opts> <in.bam>\n");
    fprintf(fp, "-m Set minimum mapping quality. Default: 0.\n");
    fprintf(fp, "-f Set minimum family size. Default: 0.\n");
    exit(exit_status);
    return exit_status;
}

int famstats_fm_main(int argc, char *argv[])
{
    famstats_t *s;
    int c;


    famstat_settings_t *settings = (famstat_settings_t *)calloc(1, sizeof(famstat_settings_t));
    settings->notification_interval = 1000000uL;

    while ((c = getopt(argc, argv, "m:f:n:h?")) >= 0) {
        switch (c) {
        case 'm':
            settings->minMQ = atoi(optarg); break;
            break;
        case 'f':
            settings->minFM = atoi(optarg); break;
            break;
        case 'n': settings->notification_interval = strtoull(optarg, NULL, 0); break;
        case '?': case 'h':
            return famstats_fm_usage_exit(stderr, EXIT_SUCCESS);
        }
    }

    if (argc != optind+1) {
        if (argc == optind) {
            famstats_fm_usage_exit(stdout, EXIT_SUCCESS);
        } else famstats_fm_usage_exit(stderr, EXIT_FAILURE);
    }

    LOG_INFO("Running main with minMQ %i and minFM %i.\n", settings->minMQ, settings->minFM);
    for(const char *tag: tags_to_check)
        check_bam_tag_exit(argv[optind], tag);

    dlib::BamHandle handle(argv[optind]);
    s = famstats_fm_core(handle, settings);
    print_stats(s, stdout);
    kh_destroy(fm, s->fm);
    kh_destroy(rc, s->rc);
    free(s);
    free(settings);
    return EXIT_SUCCESS;
}

int famstats_frac_main(int argc, char *argv[])
{
    int c;
    uint32_t minFM = 0;
    uint64_t notification_interval = 1000000;

    if(argc < 3) famstats_frac_usage_exit(stderr, EXIT_FAILURE);
    if(strcmp(argv[1], "--help") == 0) famstats_frac_usage_exit(stderr, EXIT_SUCCESS);

    while ((c = getopt(argc, argv, "n:m:h?")) >= 0) {
        switch (c) {
        case 'm':
            minFM = (uint32_t)atoi(optarg); break;
            break;
        case 'n':
            notification_interval = strtoull(optarg, NULL, 0); break;
        case '?': case 'h':
            return famstats_frac_usage_exit(stderr, EXIT_SUCCESS);
        }
    }

    if(!minFM) {
        LOG_EXIT("minFM not set. famstats_frac_main meaningless without it. Result: 1.0.\n");
    }
    LOG_INFO("Running frac main minFM %i.\n", minFM);

    if (argc != optind+1) {
        if (argc == optind) famstats_frac_usage_exit(stdout, EXIT_SUCCESS);
        else famstats_frac_usage_exit(stderr, EXIT_FAILURE);
    }
    for(const char *tag: tags_to_check) check_bam_tag_exit(argv[optind], tag);
    dlib::BamHandle handle(argv[optind]);
    uint64_t fm_above = 0, total_fm = 0, count = 0;
    // Check to see if the required tags are present before starting.
    int FM;
    int ret;
    while (LIKELY((ret = handle.next()) >= 0)) {
        // Filter reads
        if((handle.rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FREAD2)) ||
                bam_itag(handle.rec, "FP") == 0)
            continue;
        FM = bam_itag(handle.rec, "FM");
        total_fm += FM;
        if((unsigned)FM >= minFM) fm_above += FM;
        if(UNLIKELY(!(++count % notification_interval)))
            LOG_INFO("Number of records processed: %lu.\n", count);
    }
    if (ret != -1) LOG_WARNING("Truncated file? Continue anyway.\n");
    fprintf(stdout, "#Fraction of raw reads with >= minFM %i: %f.\n",
            minFM, (double)fm_above / total_fm);
    return EXIT_SUCCESS;
}

int famstats_main(int argc, char *argv[])
{
    if(argc < 2)
        return famstats_usage_exit(stderr, EXIT_FAILURE);
    if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
        return famstats_usage_exit(stderr, EXIT_SUCCESS);
    if(strcmp(argv[1], "fm") == 0)
        return famstats_fm_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "frac") == 0)
        return famstats_frac_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "target") == 0)
        return famstats_target_main(argc - 1, argv + 1);
    fprintf(stderr, "[E:%s] Unrecognized subcommand '%s'. See usage.\n", __func__, argv[1]);
    return famstats_usage_exit(stderr, EXIT_FAILURE);
}
