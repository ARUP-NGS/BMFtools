#include "dlib/bam_util.h"

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *);

int main(int argc, char **argv) {
    dlib::BamHandle in = dlib::BamHandle("bed_test.bam");
    dlib::ParsedBed bed = dlib::ParsedBed("bed_test.bed", in.header);
    bam1_t *b = bam_init1();
    size_t diffs = 0;
    void *lh3bed = bed_read("bed_test.bed");
    samFile *so = sam_open("disagreed.bam", "wb9");
    sam_hdr_write(so, in.header);
    size_t disagrees = 0, agrees = 0;
    int dbr = 0, lh3r = 0;
    while(in.read(b) != -1) {
        if(b->core.flag & (BAM_FUNMAP)) continue;
        if((dbr = bed.bam1_test(b)) != (lh3r = bed_overlap(lh3bed, in.header->target_name[b->core.tid], b->core.pos, bam_endpos(b)))) {
            LOG_EXIT("dbr: %i. lh3r: %i. Contig: %s. Position: %i. endpos; %i\n", dbr, lh3r, in.header->target_name[b->core.tid], b->core.pos, bam_endpos(b));
            if(++disagrees % 100 == 0) LOG_DEBUG("disagrees: %lu.\n", disagrees);
            sam_write1(so, in.header, b);
        } else {
            if(++agrees % 500000 == 0) LOG_DEBUG("agrees: %lu.\n", agrees);
        }
    }
    sam_close(so);
    bam_destroy1(b);
    bed_destroy(lh3bed);
    return EXIT_SUCCESS;
}
