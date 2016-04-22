#include "dlib/compiler_util.h"
#include "dlib/bam_util.h"
#include "dlib/bed_util.h"
#include <getopt.h>

namespace BMF {

    struct target_counts_t {
        uint64_t count;
        uint64_t n_skipped;
        uint64_t target;
        uint64_t rfm_count;
        uint64_t rfm_target;
    };

    int target_usage(int retcode)
    {
        fprintf(stderr,
                        "Calculates the fraction of on-target reads."
                        "Usage: bmftools target <opts> <in.bam> \n"
                        "Required arguments:\n"
                        "-b\tPath to bed.\n"
                        "Optional arguments:\n"
                        "-m\tSet minimum mapping quality for inclusion.\n"
                        "-p\tSet padding - number of bases around target region to consider as on-target. Default: 0.\n"
                        "-n\tSet notification interval - number of reads between logging statements. Default: 1000000.\n"
                );
        exit(retcode);
        return retcode; // This never happens.
    }

    target_counts_t target_core(char *bedpath, char *bampath, uint32_t padding, uint32_t minMQ, uint64_t notification_interval)
    {
        dlib::BamHandle handle(bampath);
        khash_t(bed) *bed(dlib::parse_bed_hash(bedpath, handle.header, padding));
        target_counts_t counts({0});
        uint8_t *data;
        int test;
        while (LIKELY(handle.next() >= 0)) {
            if((handle.rec->core.qual < minMQ) || (handle.rec->core.flag & (3844))) { // 3844 is unmapped, secondary, supplementary, qcfail, duplicate
                ++counts.n_skipped;
                continue;
            }
            test = dlib::bed_test(handle.rec, bed);
            counts.target += test;
            if((data = bam_aux_get(handle.rec, "FM")) != nullptr && bam_aux2i(data) > 1) {
                counts.rfm_target += test;
                ++counts.rfm_count;
            }
            if(UNLIKELY(++counts.count % notification_interval == 0)) {
                LOG_INFO("Number of records processed: %lu.\n", counts.count);
            }
        }
        dlib::bed_destroy_hash(bed);
        return counts;
    }


    int target_main(int argc, char *argv[])
    {
        char *bedpath = nullptr;
        uint32_t padding = (uint32_t)-1, minMQ = 0;
        uint64_t notification_interval = 1000000;


        if(argc < 4) return target_usage(EXIT_SUCCESS);

        if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) return target_usage(EXIT_SUCCESS);

        int c;
        while ((c = getopt(argc, argv, "m:b:p:n:h?")) >= 0) {
            switch (c) {
            case 'm':
                minMQ = strtoul(optarg, nullptr, 0); break;
            case 'b':
                bedpath = strdup(optarg); break;
            case 'p':
                padding = strtoul(optarg, nullptr, 0); break;
            case 'n':
                notification_interval = strtoull(optarg, nullptr, 0); break;
            case '?': case 'h':
                return target_usage(EXIT_SUCCESS);
            }
        }


        if(padding == (uint32_t)-1) {
            padding = DEFAULT_PADDING;
            LOG_INFO("Padding not set. Setting to default (%u).\n", padding);
        }

        if (argc != optind+1)
            return target_usage((argc == optind) ?  EXIT_SUCCESS: EXIT_FAILURE);

        target_counts_t counts = target_core(bedpath, argv[optind], padding, minMQ, notification_interval);

        if(!bedpath) {
            fprintf(stderr, "[E:%s] Bed path required for bmftools target. See usage.\n", __func__);
            return target_usage(EXIT_FAILURE);
        }

        // Open [smt]am file.

        if(bedpath) free(bedpath);

        LOG_INFO("Number of reads skipped: %lu.\n", counts.n_skipped);
        LOG_INFO("Number of reads total: %lu.\n", counts.count);
        LOG_INFO("Number of reads on target: %lu.\n", counts.target);
        LOG_INFO("Number of real FM reads total: %lu.\n", counts.rfm_count);
        LOG_INFO("Number of real FM reads on target: %lu.\n", counts.rfm_target);
        fprintf(stdout, "Fraction of reads on target with padding of %u bases and %i minMQ: %0.12f.\n",
                padding, minMQ, (double)counts.target / counts.count);
        if(counts.rfm_count) {
            fprintf(stdout, "Fraction of families of size >= 2 on target with padding of %u bases and %i minMQ: %0.12f.\n",
                    padding, minMQ, (double)counts.rfm_target / counts.rfm_count);
        }
        LOG_INFO("Successfully complete bmftools target!\n");
        return EXIT_SUCCESS;
    }

}
