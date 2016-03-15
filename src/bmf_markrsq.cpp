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
        struct mark_conf_t {
            samFile *fp;
            bam_hdr_t *hdr;
            samFile *ofp;
            mark_conf_t() :
                fp(nullptr),
                hdr(nullptr),
                ofp(nullptr)
            {
            }
            ~mark_conf_t() {
                if(fp) sam_close(fp);
                if(hdr) bam_hdr_destroy(hdr);
                if(ofp) sam_close(ofp);
            }
        } m;
        struct sort_conf_t {
            std::string pipe_name;
            std::string sortmem; // Memory per thread
            std::string tmp_prefix;
            uint32_t threads:8;
            sort_conf_t() :
                pipe_name(""),
                sortmem("500M"),
                tmp_prefix(""),
                threads(1)
            {
            }
        } s;
        struct rsq_conf_t {
            std::string pipe_name;
            samFile *fp;
            samFile *ofp;
            uint32_t mismatch_limit:4; // Maximum number of barcode mismatches for a rescue
            uint32_t use_unclipped_start:1; // Use unclipped start for rescue rather than pos
            rsq_conf_t() :
                pipe_name(""),
                fp(nullptr),
                ofp(nullptr),
                mismatch_limit(2),
                use_unclipped_start(1)
            {
            }
        } r;
        FILE *pipe_call;
        std::string final_outname;
        int level;
        markrsq_conf_t():
            m(),
            s(),
            r(),
            pipe_call(nullptr),
            final_outname(""),
            level(6)
        {
        }
        void open_pipes(char *infname);
    };

    /*
     * opens named pipes for mark -> sort + sort -> rsq
     */
    void markrsq_conf_t::open_pipes(char *infname) {
        if(s.pipe_name.empty()) {
            s.pipe_name.reserve(21uL);
            dlib::rand_string(const_cast<char *>(s.pipe_name.data()), 20uL);
            s.pipe_name.resize(20uL);
            s.pipe_name += infname; // Salt in case multiple instances are running in the same folder.
        }
        if(r.pipe_name.empty()) {
            r.pipe_name.reserve(21uL);
            dlib::rand_string(const_cast<char *>(r.pipe_name.data()), 20uL);
            r.pipe_name.resize(20uL);
            r.pipe_name += infname; // Salt in case multiple instances are running in the same folder.
        }
        if(s.tmp_prefix.empty()) {
            s.tmp_prefix.reserve(21uL);
            dlib::rand_string(const_cast<char *>(s.tmp_prefix.data()), 20uL);
            s.tmp_prefix.resize(20uL);
            s.tmp_prefix += infname; // Salt in case multiple instances are running in the same folder.
        }
        if((m.fp = sam_open(infname, "r")) == nullptr)
            LOG_EXIT("Could not open input sam %s.\n", infname);
        if((m.hdr = sam_hdr_read(m.fp)) == nullptr)
            LOG_EXIT("Could not read input sam %s's header.\n", infname);
        if(mkfifo(s.pipe_name.c_str(), 0666))
            LOG_EXIT("Could not open pipe %s.\n", s.pipe_name.c_str());
        if((m.ofp = sam_open(s.pipe_name.c_str(), "wb0")) == nullptr)
            LOG_EXIT("Could not open temporary pipe with name %s from htslib.\n", s.pipe_name.c_str());
        sam_hdr_write(m.ofp, m.hdr);
        if(mkfifo(r.pipe_name.c_str(), 0666))
            LOG_EXIT("Could not open pipe %s.\n", r.pipe_name.c_str());
        kstring_t ks = {0, 0, nullptr};
        ksprintf(&ks, "samtools sort -T%s -Obam -m%s -@%i -l0 -o%s %s",
                 s.tmp_prefix.c_str(),
                 s.sortmem.c_str(),
                 (int)s.threads,
                 r.pipe_name.c_str(),
                 s.pipe_name.c_str()
                 );
        std::string command(ks.s), free(ks.s);
        LOG_DEBUG("Pipe command: %s.\n", command.c_str());
        pipe_call = popen(command.c_str(), "w");
        r.fp = sam_open(r.pipe_name.c_str(), "r");
        r.ofp = sam_open(final_outname.c_str(), ("wb"s + std::to_string(level)).c_str());
        sam_hdr_write(r.ofp, m.hdr);
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
    conf.open_pipes(argv[optind]);;
    int mark_ret = mark_core(&conf);
    if(mark_ret) {

    }
    LOG_INFO("Successfully complete bmftools stack!\n");
    return 0;
}

}
