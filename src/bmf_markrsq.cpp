/*
 * OUTLINE
 * 1. Mark/prepare for sorting.
 *   1. Handle supplementary/secondary
 *     1. How?
 *       1. Make a stack for reads with a given read name
 *         1. Add tags for SU/MU/LM to all reads in the set such that they have the same keys
 *         2. Add a tag to read 1 and read 2 for each (if any) of its supplemental alignments to know what to wait for.
 *           1. Add a tag to the supplementals for their actual position/contig, flag as unmapped, and move to their primaries.
 *         3. Output to stdout
 * 2. Pass to samtools sort and sort by coordinate.
 * 3. Load in a buffer of reads
 *   1. Fill a large stack of buffered reads.
 *   2. Build a hashmap for r1/r2 combinations.
 *   3. Once all reads needed for an alignment signature set have been loaded, collapse them, putting the supplementals in a separate table.
 *     1. If a read name set is not collapsed and there are supplementals, unset the unmapped flag and change the coordinates back to what they should be.
 *     2. Otherwise, ignore the supplementals because they'll be realigned.
 */
#include "bmf_markrsq.h"
#include <getopt.h>

namespace BMF {

    struct markrsq_conf_t {
        struct sort_conf_t {
            std::string pipe_name;
            std::string sortmem; // Memory per thread
            std::string tmp_prefix;
            uint32_t threads:8;
            uint32_t level:4; // Output compression level
            sort_conf_t() :
                sortmem("500M"),
                tmp_prefix(""),
                threads(1),
                level(0)
            {
            }
        } s;
        struct rsq_conf_t {
            std::string pipe_name;
            uint32_t mismatch_limit:4; // Maximum number of barcode mismatches for a rescue
            uint32_t use_unclipped_start:1; // Use unclipped start for rescue rather than pos
            rsq_conf_t() :
                mismatch_limit(2),
                use_unclipped_start(1)
            {
            }
        } r;
        markrsq_conf_t():
            s(),
            r()
        {
        }
    };

    /*
     * opens named pipes for mark -> sort + sort -> rsq
     */
    FILE *pipe_call(markrsq_conf_t *conf) {
        if(conf->s.pipe_name.empty()) {
            conf->s.pipe_name.reserve(21uL);
            dlib::rand_string(const_cast<char *>(conf->s.pipe_name.data()), 20uL);
            conf->s.pipe_name.resize(20uL);
        }
        if(conf->r.pipe_name.empty()) {
            conf->r.pipe_name.reserve(21uL);
            dlib::rand_string(const_cast<char *>(conf->r.pipe_name.data()), 20uL);
            conf->r.pipe_name.resize(20uL);
        }
        kstring_t ks{0};
        ksprintf(&ks, "mkfifo %s", conf->s.pipe_name.c_str());
        int retcode;
        if((retcode = pclose(popen(ks.s, "w"))) != 0)
            LOG_EXIT("mkfifo call for sort failed. Abort!\n");
        ksprintf(&ks, "mkfifo %s", conf->r.pipe_name.c_str());
        if((retcode = pclose(popen(ks.s, "w"))) != 0)
            LOG_EXIT("mkfifo call for rsq failed. Abort!\n");
        ksprintf(&ks, "samtools sort -O bam -m%s -@%i -l%i - < %s > %s",
                 conf->s.sortmem.c_str(),
                 (int)conf->s.threads,
                 (int)conf->s.level,
                 conf->s.pipe_name.c_str(),
                 conf->r.pipe_name.c_str()
                 );
        std::string command(ks.s), free(ks.s);
        LOG_DEBUG("Pipe command: %s.\n", command.c_str());
        return popen(command.c_str(), "w");
    }

    void markrsq_usage(int retcode)
    {
        fprintf(stderr,
                        "Unwritten usage. Eh.\n"
                );
        exit(retcode);
    }

    /*
     * Take in a name-sorted bam, marks reads appropriately, and builds a
     * hashmap of information for read pairs for the rescue step.
     * Writes
     */
    int mark_core(markrsq_conf_t* conf) {
        return 0;
    }

int markrsq_main(int argc, char *argv[]) {
    int c;
    if(argc < 2) markrsq_usage(EXIT_FAILURE);
    struct markrsq_conf_t conf;
    const struct option lopts[] = {
        {0, 0, 0, 0}
    };
    while ((c = getopt_long(argc, argv, "R:D:q:r:2:S:d:a:s:m:p:f:b:v:o:O:c:BP?hV", lopts, nullptr)) >= 0) {
        switch (c) {
            case 'T': /* */ break;
            case 'h': case '?': markrsq_usage(EXIT_SUCCESS);
        }
    }
    if(optind >= argc - 1) LOG_EXIT("Insufficient arguments. Input bam required!\n");
    int mark_ret = mark_core(&conf);
    if(mark_ret) {

    }
    LOG_INFO("Successfully complete bmftools stack!\n");
    return 0;
}

}
