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
            size_t stack_start;
            mark_conf_t() :
                fp(nullptr),
                hdr(nullptr),
                ofp(nullptr),
                stack_start(10)
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
            std::unordered_map<std::string, size_t> stack_sizes;
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
                use_unclipped_start(0)
            {
            }
            ~rsq_conf_t() {
                if(fp) sam_close(fp), fp = nullptr;
                sam_close(ofp), ofp = nullptr;
            }
        } r;
        FILE *pipe_call;
        int level;
        markrsq_conf_t():
            m(),
            s(),
            r(),
            pipe_call(nullptr),
            level(6)
        {
        }
        void open_pipes(char *infname, char *outfname);
        void mark_core();
        void rsq_core();
        void process(char *infname, char *outfname) {
            open_pipes(infname, outfname);
            mark_core();
            rsq_core();
        }
    };

    /*
     * opens named pipes for mark -> sort + sort -> rsq
     */
    void markrsq_conf_t::open_pipes(char *infname, char *outfname) {
        if(s.pipe_name.empty()) {
            char buffer[32];
            dlib::rand_string(buffer, 20uL);
            s.pipe_name = buffer;
            s.pipe_name += infname; // Salt in case multiple instances are running in the same folder.
        }
        if(r.pipe_name.empty()) {
            char buffer[32];
            dlib::rand_string(buffer, 20uL);
            r.pipe_name = buffer;
            r.pipe_name += infname; // Salt in case multiple instances are running in the same folder.
        }
        if(s.tmp_prefix.empty()) {
            char buffer[32];
            dlib::rand_string(buffer, 20uL);
            s.tmp_prefix = buffer;
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
        r.ofp = sam_open(outfname, ("wb"s + std::to_string(level)).c_str());
        sam_hdr_write(r.ofp, m.hdr);
    }

    /*
     * stack: stack of bam records which all share read names.
     * return: "tid1:ucs1:tid2:ucs2:is_rev1:is_rev2:NumberSupplementary"
     */
    void bam1_stack::make_hashmap_key(std::string &ret) {
        std::vector<int> supplemental_coordinates1;
        std::vector<int> supplemental_coordinates2;
#if !NDEBUG
        for(unsigned i = 0; i < n; ++i) {
            // All reads in stack should have the same name.
            assert(strcmp(bam_get_qname(recs + i), bam_get_qname(recs + i)) == 0);
        }
#endif
        int r1len = -1, r2len;
        int ucs1(-1), ucs2(-1);
        int tid1(-1), tid2(-1);
        int is_rev1(-1), is_rev2(-1);
        for(unsigned i = 0; i < n; ++i) {
            if(((recs + i)->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) == 0) {// Primary alignment
                if((recs + i)->core.flag & BAM_FREAD1) {
                    r1len = (recs + i)->core.l_qseq;
                    ucs1 = dlib::get_unclipped_start((recs + i));
                    tid1 = (recs + i)->core.tid;
                    is_rev1 = (recs + i)->core.flag & BAM_FREVERSE;
                    is_rev2 = (recs + i)->core.flag & BAM_FMREVERSE;
                } else {
                    r2len = (recs + i)->core.l_qseq;
                    ucs2 = dlib::get_unclipped_start((recs + i));
                    tid2 = (recs + i)->core.tid;
                }
            }
            else {
                if((recs + i)->core.flag & BAM_FREAD1) {
                    supplemental_coordinates1.reserve(supplemental_coordinates1.capacity() + 2);
                    supplemental_coordinates1.push_back((recs + i)->core.tid);
                    supplemental_coordinates1.push_back((recs + i)->core.pos);
                } else {
                    supplemental_coordinates2.reserve(supplemental_coordinates2.capacity() + 2);
                    supplemental_coordinates2.push_back((recs + i)->core.tid);
                    supplemental_coordinates2.push_back((recs + i)->core.pos);
                }
            }
        }
        if(r1len < 0) {
            LOG_EXIT("Could not find a primary read 1 alignment. Is the dataset paired?\n");
        }
        dlib::safe_stringprintf(ret, "%i:%i:%i:%i:%i:%i%lu",
                                tid1,
                                ucs1,
                                tid2,
                                ucs2,
                                is_rev1,
                                is_rev2,
                                n - 2
                                );
        for(unsigned i = 0; i < n; ++i) {
            bam_aux_append((recs + i), "as", 'Z', ret.size() + 1, (uint8_t *)ret.data()); // Add the key as a tag
            bam_aux_append((recs + i), "S1", 'B', sizeof(int) * supplemental_coordinates1.size(), reinterpret_cast<uint8_t *>(supplemental_coordinates1.data()));
            bam_aux_append((recs + i), "S2", 'B', sizeof(int) * supplemental_coordinates2.size(), reinterpret_cast<uint8_t *>(supplemental_coordinates2.data()));
            bam_aux_append((recs + 1), "L1", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&r1len));
            bam_aux_append((recs + 1), "L2", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&r2len));
            if((recs + i)->core.flag & BAM_FREAD1) {
                bam_aux_append((recs + i), "SU", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&ucs1));
                bam_aux_append((recs + i), "MU", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&ucs2));
                // Add a tag to tell us where the supplementary alignments are?
            } else {
                bam_aux_append((recs + i), "SU", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&ucs2));
                bam_aux_append((recs + i), "MU", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&ucs1));
            }
        }
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
    void markrsq_conf_t::mark_core() {
        if(!m.fp || !m.hdr) LOG_EXIT("Open pipes before running mark_core!\n");

        bam1_stack stack{(m.stack_start)};
        int c;
        bam1_t *b = bam_init1();
        std::string name_ptr;
        std::string hash_key;
        hash_key.reserve(128uL);
        while((c = sam_read1(m.fp, m.hdr, b) >= 0)) {
            if(name_ptr.empty()) {
                name_ptr = bam_get_qname(b);
                stack.add_rec(b);
                continue;
            } else {
                if(strcmp(bam_get_qname(b), name_ptr.c_str()) == 0) {
                    stack.add_rec(b);
                } else {
                    // Add to key hash
                    stack.make_hashmap_key(hash_key);
                    auto found = r.stack_sizes.find(hash_key);
                    if(found == r.stack_sizes.end()) {
                        ++found->second;
                    } else {
                        r.stack_sizes[hash_key] = 1;
                    }
                }
            }
        }
        bam_destroy1(b);

        // end
        sam_close(m.fp), m.fp = nullptr;
        sam_close(m.ofp), m.ofp = nullptr;
    }

    void markrsq_conf_t::rsq_core() {
        // Clean up after.
        int ret;
        if((ret = pclose(pipe_call)) != 0) {
            LOG_EXIT("Pipe call failed: %i.\n", ret);
        }
        m.~mark_conf_t();
        s.~sort_conf_t();
        r.~rsq_conf_t();
    }

int markrsq_main(int argc, char *argv[]) {
    int c;
    if(argc < 2) markrsq_usage(EXIT_FAILURE);
    struct markrsq_conf_t conf;
    const struct option lopts[] = {
        {0, 0, 0, 0}
    };
    char *outfname(nullptr);
    while ((c = getopt_long(argc, argv, "m:o:t:T:l:u?h", lopts, nullptr)) >= 0) {
        switch (c) {
            case 'L': conf.r.mismatch_limit = atoi(optarg); break;
            case 't': conf.s.threads = atoi(optarg); break;
            case 'm': conf.s.sortmem = optarg; break;
            case 'T': conf.s.tmp_prefix = optarg; break;
            case 'o': outfname = optarg; break;
            case 'u': conf.r.use_unclipped_start = 1; break;
            case 'l': conf.level = atoi(optarg) % 10; break;
            case 'h': case '?': markrsq_usage(EXIT_SUCCESS);
        }
    }
    if(optind >= argc - 1) LOG_EXIT("Insufficient arguments. Input bam required!\n");
    conf.process(argv[optind], outfname ? outfname: const_cast<char *>("-"));
    LOG_INFO("Successfully complete bmftools stack!\n");
    return 0;
}

}
