#include "bmf_rsq.h"
#include <cstring>
#include <getopt.h>
#include "dlib/cstr_util.h"
#include "include/igamc_cephes.h" /// for igamc
#include <algorithm>

namespace bmf {

static const char SO_STR[]{"positional_rescue"};

static const int sp(1);

struct rsq_aux_t {
    FILE *fqh;
    samFile *in;
    samFile *out;
    uint32_t mmlim:6;
    uint32_t is_se:1;
    uint32_t write_supp:1; // Write reads with supplementary alignments
    uint32_t infer:1; // Use inference instead of barcodes.
    uint32_t trust_unmasked:1;
    bam_hdr_t *hdr; // BAM header
    std::unordered_map<std::string, std::string> realign_pairs;
};

inline void bam2ffq(bam1_t *b, FILE *fp, const int is_supp=0);
inline void add_dummy_tags(bam1_t *b);

void update_bam1(bam1_t *p, bam1_t *b);
void update_bam1_unmasked(bam1_t *p, bam1_t *b);

template<int (*fn)(bam1_t *, bam1_t *)>
struct Stack {
    uint16_t mmlim:8;
    uint16_t trust_unmasked:1;
    uint16_t infer:1;
    unsigned n; // Number used
    unsigned m; // Maximum allocated
    bam1_t *a; // Array
    bam1_t **stack; // Pointers to reads.

    Stack(rsq_aux_t *settings, unsigned _m=0):
            mmlim(settings->mmlim),
            trust_unmasked(settings->trust_unmasked),
            infer(settings->infer),
            n(0),
            m(_m),
            a((bam1_t *)calloc(m, sizeof(bam1_t))),
            stack((bam1_t **)malloc(m * sizeof(bam1_t *)))
    {
        for(unsigned i(0); i < m; ++i) stack[i] = a + i;
    }
    ~Stack() {
        LOG_DEBUG("m: %u.\n", m);
        for(unsigned i(0); i < m; ++i) {
            LOG_DEBUG("Cleaning up.\n")
            if(a[i].data) {
                LOG_DEBUG("NONZERO %p\n", (void *)a[i].data);
                free(a[i].data);
            }
        }
        free(stack);
        free(a);
    }
    void add(const bam1_t *b) {
        if(n + 1 >= m) {
            m <<= 1;
            LOG_DEBUG("Max increased to %lu.\n", m);
            a = (bam1_t *)realloc(a, sizeof(bam1_t) * m); //
            stack = (bam1_t **)realloc(stack, sizeof(bam1_t *) * m); //
            memset(a + n, 0, (m - n) * sizeof(bam1_t)); // Zero-initialize later records.
            for(unsigned i(n); i < m; ++i) stack[i] = a + i;
            LOG_DEBUG("Finished adding.\n");
        }
        bam_copy1(a + n++, b);
    }
    void clear() {
        for(unsigned i(0); i < n; ++i) free((a + i)->data);
        memset(a, 0, n * sizeof(bam1_t));
        n = 0;
    }
    void write_stack_pe(rsq_aux_t *settings);
    void write_stack_se(rsq_aux_t *settings);
    void flatten();
    inline void flatten_infer();
    void pe_core(rsq_aux_t *settings);
    void pe_core_infer(rsq_aux_t *settings);
    void se_core(rsq_aux_t *settings);
    void se_core_infer(rsq_aux_t *settings);
};

template<int (*fn)(bam1_t *, bam1_t *)>
void Stack<fn>::se_core_infer(rsq_aux_t *settings) {
    // This selects the proper function to use for deciding if reads belong in the same stack.
    // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
    if(strcmp(dlib::get_SO(settings->hdr).c_str(), SO_STR))
        LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                 dlib::get_SO(settings->hdr).c_str(), SO_STR);
    bam1_t *b(bam_init1());
    uint64_t count(0);
    while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
        if(UNLIKELY(++count % 1000000 == 0)) LOG_INFO("Records read: %lu.\n", count);
        add_dummy_tags(b);
        if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if(b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) {
            sam_write1(settings->out, settings->hdr, b); continue;
        }
        //LOG_DEBUG("Read a read!\n");
        if(fn(b, a) == 0) write_stack_se(settings); // Flattens and clears stack.
        add(b);
    }
    write_stack_se(settings);
    bam_destroy1(b);
    // Handle any unpaired reads, though there shouldn't be any in real datasets.
    LOG_DEBUG("Number of orphan reads: %lu.\n", settings->realign_pairs.size());
    if(settings->realign_pairs.size()) {
        LOG_WARNING("There shouldn't be orphan reads in real datasets. Number found: %lu\n", settings->realign_pairs.size());
#if !NDEBUG
        for(auto& pair: settings->realign_pairs)
            puts(pair.second.c_str());
#endif
    }
}

template<int (*fn)(bam1_t *, bam1_t *)>
void Stack<fn>::se_core(rsq_aux_t *settings) {
    // This selects the proper function to use for deciding if reads belong in the same stack.
    // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
    if(infer) return se_core_infer(settings);
    if(strcmp(dlib::get_SO(settings->hdr).c_str(), SO_STR))
        LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                 dlib::get_SO(settings->hdr).c_str(), SO_STR);
    bam1_t *b(bam_init1());
    uint64_t count(0);
    while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
        if(UNLIKELY(++count % 1000000 == 0)) LOG_INFO("Records read: %lu.\n", count);
        if(b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) {
            sam_write1(settings->out, settings->hdr, b);
            continue;
        }
        if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        //LOG_DEBUG("Read a read!\n");
        if(fn(b, a) == 0) write_stack_se(settings); // Flattens and clears stack.
        add(b);
    }
    write_stack_se(settings);
    bam_destroy1(b);
    // Handle any unpaired reads, though there shouldn't be any in real datasets.
    LOG_DEBUG("Number of orphan reads: %lu.\n", settings->realign_pairs.size());
    if(settings->realign_pairs.size()) {
        LOG_WARNING("There shouldn't be orphan reads in real datasets. Number found: %lu\n", settings->realign_pairs.size());
#if !NDEBUG
        for(auto pair: settings->realign_pairs)
            puts(pair.second.c_str());
#endif
    }
}

template<int (*fn)(bam1_t *, bam1_t *)>
void Stack<fn>::pe_core_infer(rsq_aux_t *settings)
{
    // This selects the proper function to use for deciding if reads belong in the same stack.
    // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
    if(strcmp(dlib::get_SO(settings->hdr).c_str(), SO_STR))
        LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                 dlib::get_SO(settings->hdr).c_str(), SO_STR);
    bam1_t *b(bam_init1());
    uint64_t count(0);
    while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
        if(UNLIKELY(++count % 1000000 == 0)) LOG_INFO("Records read: %lu.\n", count);
        add_dummy_tags(b);
        if(b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) {
            sam_write1(settings->out, settings->hdr, b);
            continue;
        }
        if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
            continue;
        if(n == 0 || fn(b, a) == 0)
            write_stack_pe(settings); // Flattens and clears stack.
        add(b);
    }
    write_stack_pe(settings);
    bam_destroy1(b);
    // Handle any unpaired reads, though there shouldn't be any in real datasets.
    LOG_DEBUG("Number of orphan reads: %lu.\n", settings->realign_pairs.size());
    if(settings->realign_pairs.size()) {
        LOG_WARNING("There shouldn't be orphan reads in real datasets. Number found: %lu\n", settings->realign_pairs.size());
#if !NDEBUG
        for(auto pair: settings->realign_pairs)
            puts(pair.second.c_str());
#endif
    }
}

template<int (*fn)(bam1_t *, bam1_t *)>
void Stack<fn>::pe_core(rsq_aux_t *settings)
{
    // This selects the proper function to use for deciding if reads belong in the same stack.
    // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
    if(strcmp(dlib::get_SO(settings->hdr).c_str(), SO_STR))
        LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                 dlib::get_SO(settings->hdr).c_str(), SO_STR);
    bam1_t *b(bam_init1());
    uint64_t count(0);
    while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
        if(UNLIKELY(++count % 1000000 == 0)) LOG_INFO("Records read: %lu.\n", count);
        if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if(b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) {
            sam_write1(settings->out, settings->hdr, b);
            continue;
        }
        //LOG_DEBUG("Read a read!\n");
        if(fn(b, a) == 0) write_stack_pe(settings); // Flattens and clears stack.
        add(b);
    }
    write_stack_pe(settings);
    bam_destroy1(b);
    // Handle any unpaired reads, though there shouldn't be any in real datasets.
    LOG_DEBUG("Number of orphan reads: %lu.\n", settings->realign_pairs.size());
    if(settings->realign_pairs.size()) {
        LOG_WARNING("There shouldn't be orphan reads in real datasets. Number found: %lu\n", settings->realign_pairs.size());
#if !NDEBUG
        for(auto pair: settings->realign_pairs)
            puts(pair.second.c_str());
#endif
    }
}

template<int (*fn)(bam1_t *, bam1_t *)>
inline void Stack<fn>::flatten_infer()
{
    unsigned i, j;
    for(i = 0; i < n; ++i) stack[i] = a + i;
    std::sort(stack, stack + n, [](bam1_t *a, bam1_t *b) {
            return a ? (b ? 0: 1): b ? strcmp(bam_get_qname(a), bam_get_qname(b)): 0;
    });
    for(i = 0; i < n; ++i) {
        for(j = i + 1; j < n; ++j) {
            //assert(key == ucs_sort_core_key(a[j]));
            //assert(ucs == dlib::get_unclipped_start(a[j]));
            //assert(tid == a[j]->core.tid);
            //assert(mucs == bam_itag(a[j], "MU"));
            if(stack[i]->core.l_qseq != stack[j]->core.l_qseq ||
               stack[i]->core.l_qname != stack[j]->core.l_qname)
                continue;
            //LOG_DEBUG("Flattening %s into %s.\n", bam_get_qname(a[i]), bam_get_qname(a[j]));
            if(trust_unmasked) update_bam1_unmasked(a + j, a + i);
            else update_bam1(a + j, a + i);
            (a + i)->data = nullptr;
            break;
            // "break" in case there are multiple within hamming distance.
            // Otherwise, I'll end up having memory mistakes.
            // Besides, that read set will get merged into the later read in the set.
        }
    }
}

template<int (*fn)(bam1_t *, bam1_t *)>
void Stack<fn>::flatten()
{
    if(infer) return flatten_infer();
    unsigned i, j;
    for(i = 0; i < n; ++i) stack[i] = a + i;
    std::sort(stack, stack + n, [](const bam1_t *a, const bam1_t *b) {
            return a ? (b ? 0: 1): b ? strcmp(bam_get_qname(a), bam_get_qname(b)): 0;
    });
    for(i = 0; i < n; ++i) {
        for(j = i + 1; j < n; ++j) {
            //assert(key == ucs_sort_core_key(stack->a[j]));
            //assert(ucs == dlib::get_unclipped_start(stack->a[j]));
            //assert(tid == stack->a[j]->core.tid);
            //assert(mucs == bam_itag(stack->a[j], "MU"));
            if(stack[i]->core.l_qseq != stack[j]->core.l_qseq)
                continue;
            if(stack[i]->core.l_qname != stack[j]->core.l_qname)
                continue;
            if(dlib::stringhd(bam_get_qname(stack[i]), bam_get_qname(stack[j])) <= mmlim) {
                assert(dlib::stringhd(bam_get_qname(stack[i]), bam_get_qname(stack[j])) <= mmlim);
                //LOG_DEBUG("Flattening %s into %s.\n", bam_get_qname(stack[i]), bam_get_qname(stack[j]));
                if(trust_unmasked) update_bam1_unmasked(stack[j], stack[i]);
                else update_bam1(stack[j], stack[i]);
                free(stack[i]->data);
                stack[i]->data = nullptr;
                break;
                // "break" in case there are multiple within hamming distance.
                // Otherwise, I'll end up having memory mistakes.
                // Besides, that read set will get merged into the later read in the set.
            }
        }
    }
}

template<int (*fn)(bam1_t *, bam1_t *)>
void Stack<fn>::write_stack_se(rsq_aux_t *settings)
{
    flatten();
    uint8_t *data;
    for(unsigned i(0); i < n; ++i) {
        if((a + i)->data) {
            if((data = bam_aux_get((a + i), "NC")))
                bam2ffq((a + i), settings->fqh);
            else
                sam_write1(settings->out, settings->hdr, (a + i));
        }
    }
    clear();
}

template<int (*fn)(bam1_t *, bam1_t *)>
void Stack<fn>::write_stack_pe(rsq_aux_t *settings)
{
    flatten();
    if(settings->is_se) return write_stack_se(settings);
    //size_t n = 0;
    //LOG_DEBUG("Starting to write stack\n");
    uint8_t *data;
    std::string qname;
    for(unsigned i(0); i < n; ++i) {
        if(a[i].data) {
            if((data = bam_aux_get(a + i, "NC"))) {
                //LOG_DEBUG("Trying to write.\n");
                qname = bam_get_qname(a + i);
                if(settings->realign_pairs.find(qname) == settings->realign_pairs.end()) {
                    settings->realign_pairs.emplace(qname, dlib::bam2cppstr(a + i));
                } else {
                    // Write read 1 out first.
                    if((a + i)->core.flag & BAM_FREAD2) {
                        fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                        bam2ffq((a + i), settings->fqh);
                    } else {
                        bam2ffq((a + i), settings->fqh);
                        fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                    }
                    // Clear entry, as there can only be two.
                    settings->realign_pairs.erase(qname);
                }
            } else if(settings->write_supp & (bam_aux_get((a + i), "SA") || bam_aux_get((a + i), "ms"))) {
                //LOG_DEBUG("Trying to write write supp or stuff.\n");
                // Has an SA or ms tag, meaning that the read or its mate had a supplementary alignment
                qname = bam_get_qname((a + i));
                bam_aux_append(a + i, "SP", 'i', sizeof(int), const_cast<uint8_t *>(reinterpret_cast<const uint8_t*>(&sp)));
                if(settings->realign_pairs.find(qname) == settings->realign_pairs.end()) {
                    settings->realign_pairs[qname] = dlib::bam2cppstr((a + i));
                } else {
                    // Make sure the read names/barcodes match.
                    //assert(memcmp(settings->realign_pairs[qname].c_str() + 1, qname.c_str(), qname.size() - 1) == 0);
                    // Write read 1 out first.
                    if((a + i)->core.flag & BAM_FREAD2) {
                        fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                        //bam2ffq((a + i), settings->fqh, qname);
                        bam2ffq((a + i), settings->fqh, 1);
                    } else {
                        //bam2ffq((a + i), settings->fqh, qname);
                        bam2ffq((a + i), settings->fqh, 1);
                        fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                    }
                    // Clear entry, as there can only be two.
                    settings->realign_pairs.erase(qname);
                }
            } else {
                for(const char *tag: {"MU", "ms", "LM"})
                    if((data = bam_aux_get((a + i), tag)))
                        bam_aux_del((a + i), data);
                sam_write1(settings->out, settings->hdr, (a + i));
            }
        }
    }
    clear();
}

inline void bam2ffq(bam1_t *b, FILE *fp, const int is_supp)
{
    int i;
    uint8_t *rvdata;
    kstring_t ks{0, 120uL, (char *)malloc(120uL)};
    kputc('@', &ks);
    kputsn(bam_get_qname(b), b->core.l_qname - 1, &ks);
    kputsnl(" PV:B:I", &ks);
    auto fa((uint32_t *)dlib::array_tag(b, "FA"));
    auto pv((uint32_t *)dlib::array_tag(b, "PV"));
    for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
    kputsnl("\tFA:B:I", &ks);
    for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
    ksprintf(&ks, "\tFM:i:%i\tFP:i:%i", bam_itag(b, "FM"), bam_itag(b, "FP"));
    write_if_found(rvdata, b, "RV", ks);
    write_if_found(rvdata, b, "NC", ks);
    write_if_found(rvdata, b, "DR", ks);
    write_if_found(rvdata, b, "NP", ks);
    if(is_supp) kputsnl("\tSP:i:1", &ks);
    kputc('\n', &ks);
    uint8_t *seq(bam_get_seq(b));
    char *seqbuf((char *)malloc(b->core.l_qseq + 1));
    for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[bam_seqi(seq, i)];
    seqbuf[i] = '\0';
    if (b->core.flag & BAM_FREVERSE) { // reverse complement
        for(i = 0; i < b->core.l_qseq>>1; ++i) {
            const int8_t t(nuc_cmpl(seqbuf[b->core.l_qseq - i - 1]));
            seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seqbuf[i]);
            seqbuf[i] = t;
        }
        if(b->core.l_qseq&1) seqbuf[i] = nuc_cmpl(seqbuf[i]);
    }
    seqbuf[b->core.l_qseq] = '\0';
    assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
    kputsn(seqbuf, b->core.l_qseq, &ks);
    kputsnl("\n+\n", &ks);
    uint8_t *qual(bam_get_qual(b));
    for(i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = 33 + qual[i];
    if (b->core.flag & BAM_FREVERSE) { // reverse
        for (i = 0; i < b->core.l_qseq>>1; ++i) {
            const int8_t t(seqbuf[b->core.l_qseq - 1 - i]);
            seqbuf[b->core.l_qseq - 1 - i] = seqbuf[i];
            seqbuf[i] = t;
        }
    }
    kputsn(seqbuf, b->core.l_qseq, &ks), free(seqbuf);
    kputc('\n', &ks);
    fputs(ks.s, fp), free(ks.s);
}


inline int switch_names(char *n1, char *n2) {
    for(;*n1;++n1, ++n2) if(*n1 != *n2) return *n1 < *n2;
    return 0; // If identical, don't switch. Should never happen.
}


void update_bam1_unmasked(bam1_t *p, bam1_t *b)
{
    uint8_t *bdata(bam_aux_get(b, "FM"));
    uint8_t *pdata(bam_aux_get(p, "FM"));
    if(UNLIKELY(!bdata || !pdata)) {
        fprintf(stderr, "Required FM tag not found. Abort mission!\n");
        exit(EXIT_FAILURE);
    }
    int bFM(bam_aux2i(bdata));
    int pFM(bam_aux2i(pdata));
    int pTMP(0);
    if(switch_names(bam_get_qname(p), bam_get_qname(b))) {
        memcpy(bam_get_qname(p), bam_get_qname(b), b->core.l_qname);
        assert(strlen(bam_get_qname(p)) == strlen(bam_get_qname(b)));
    }
    pFM += bFM;
    bam_aux_del(p, pdata);
    bam_aux_append(p, "FM", 'i', sizeof(int), (uint8_t *)&pFM);
    if((pdata = bam_aux_get(p, "RV")) != nullptr) {
        pTMP = bam_aux2i(pdata) + bam_itag(b, "RV");
        bam_aux_del(p, pdata);
        bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pTMP);
    }
    // Handle NC (Number Changed) tag
    pdata = bam_aux_get(p, "NC");
    bdata = bam_aux_get(b, "NC");
    int n_changed(dlib::int_tag_zero(pdata) + dlib::int_tag_zero(bdata));
    const int was_merged(((!!pdata) << 1) | (!!bdata));
    if(pdata) bam_aux_del(p, pdata);
    // If the collapsed observation is now duplex but wasn't before, this updates the DR tag.
    if(pTMP != pFM && pTMP && (pdata = bam_aux_get(p, "DR")) && bam_aux2i(pdata) == 0) {
        pTMP = 1;
        bam_aux_del(p, pdata);
        bam_aux_append(p, "DR", 'i', sizeof(int), (uint8_t *)&pTMP);
    }
    if((pdata = bam_aux_get(p, "NP"))) {
        bdata = bam_aux_get(b, "NP");
        pTMP = bam_aux2i(pdata) + (bdata ? bam_aux2i(bdata) : 1);
        bam_aux_del(p, pdata);
    } else {
        bdata = bam_aux_get(b, "NP");
        pTMP = (bdata ? bam_aux2i(bdata) : 1) + 1;
    }
    bam_aux_append(p, "NP", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&pTMP));
    uint32_t *bPV((uint32_t *)dlib::array_tag(b, "PV")); // Length of this should be b->l_qseq
    uint32_t *pPV((uint32_t *)dlib::array_tag(p, "PV"));
    uint32_t *bFA((uint32_t *)dlib::array_tag(b, "FA"));
    uint32_t *pFA((uint32_t *)dlib::array_tag(p, "FA"));
    uint8_t *bSeq(bam_get_seq(b));
    uint8_t *pSeq(bam_get_seq(p));
    uint8_t *bQual(bam_get_qual(b));
    uint8_t *pQual(bam_get_qual(p));
    const int qlen(p->core.l_qseq);
    int8_t ps, bs;

    if(p->core.flag & (BAM_FREVERSE)) {
        int qleni1;
        for(int i(0); i < qlen; ++i) {
            qleni1 = qlen - i - 1;
            ps = bam_seqi(pSeq, qleni1);
            bs = bam_seqi(bSeq, qleni1);
            if(ps == bs) {
                if((pPV[i] = agreed_pvalues(pPV[i], bPV[i])) < 3) {
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[qleni1] = 2;
                } else {
                    pFA[i] += bFA[i];
                    if(bQual[qleni1] > pQual[qleni1])
                        pQual[qleni1] = bQual[qleni1];
                }
            } else if(ps == dlib::htseq::HTS_N) {
                if(was_merged & 2) {
                    n_base(pSeq, qleni1);
                    pFA[i] = 0;
                    pPV[i] = 0;
                    continue;
                }
                bam_set_base(pSeq, bSeq, qleni1);
                pFA[i] = bFA[i];
                if(bPV[i] < 3) {
                    pPV[i] = 0;
                    pFA[i] = 0;
                    pQual[qleni1] = 2;
                } else {
                    pPV[i] = bPV[i];
                    pFA[i] = bFA[i];
                    pQual[qleni1] = bQual[qleni1];
                    ++n_changed; // Note: goes from N to a useable nucleotide.
                }
            } else if(bs == dlib::htseq::HTS_N) {
                if(was_merged & 1) {
                    n_base(pSeq, qleni1);
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[i] = 2;
                }
                continue;
            } else {
                n_base(pSeq, qleni1);
                pFA[i] = 0;
                pPV[i] = 0;
                pQual[i] = 2;
            }
        }
    } else {
        for(int i(0); i < qlen; ++i) {
            ps = bam_seqi(pSeq, i);
            bs = bam_seqi(bSeq, i);
            if(ps == bs) {
                if((pPV[i] = agreed_pvalues(pPV[i], bPV[i])) > 2) {
                    pFA[i] += bFA[i];
                    if(bQual[i] > pQual[i]) pQual[i] = bQual[i];
                } else {
                    n_base(pSeq, i);
                    pPV[i] = 0;
                    pFA[i] = 0;
                }
            } else if(ps == dlib::htseq::HTS_N) {
                if(was_merged & 2) {
                    n_base(pSeq, i);
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[i] = 2;
                    continue;
                }
                if(bPV[i] > 2) {
                    bam_set_base(pSeq, bSeq, i);
                    pFA[i] = bFA[i];
                    pPV[i] = bPV[i];
                    pQual[i] = bQual[i];
                    ++n_changed; // Note: goes from N to a useable nucleotide.
                } else {
                    n_base(pSeq, i);
                    pPV[i] = 0;
                    pFA[i] = 0;
                    pQual[i] = 2;
                }
                continue;
            } else if(bs == dlib::htseq::HTS_N) {
                if(was_merged & 1) {
                    n_base(pSeq, i);
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[i] = 2;
                    continue;
                }
                if(pPV[i] < 3) {
                    n_base(pSeq, i);
                    pPV[i] = 0;
                    pFA[i] = 0;
                    pQual[i] = 2;
                }
                continue;
            } else {
                n_base(pSeq, i);
                pFA[i] = 0;
                pPV[i] = 0;
                pQual[i] = 2;
            }
        }
    }
    bam_aux_append(p, "NC", 'i', sizeof(int), (uint8_t *)&n_changed);
}

void update_bam1(bam1_t *p, bam1_t *b)
{
    uint8_t *bdata(bam_aux_get(b, "FM"));
    uint8_t *pdata(bam_aux_get(p, "FM"));
    if(UNLIKELY(!bdata || !pdata)) {
        fprintf(stderr, "Required FM tag not found. Abort mission!\n");
        exit(EXIT_FAILURE);
    }
    int bFM(bam_aux2i(bdata));
    int pFM(bam_aux2i(pdata));
    int pTMP(0);
    if(switch_names(bam_get_qname(p), bam_get_qname(b))) {
        memcpy(bam_get_qname(p), bam_get_qname(b), b->core.l_qname);
        assert(strlen(bam_get_qname(p)) == strlen(bam_get_qname(b)));
    }
    pFM += bFM;
    bam_aux_del(p, pdata);
    bam_aux_append(p, "FM", 'i', sizeof(int), (uint8_t *)&pFM);
    if((pdata = bam_aux_get(p, "RV")) != nullptr) {
        pTMP = bam_aux2i(pdata) + bam_itag(b, "RV");
        bam_aux_del(p, pdata);
        bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pTMP);
    }
    // Handle NC (Number Changed) tag
    pdata = bam_aux_get(p, "NC");
    bdata = bam_aux_get(b, "NC");
    int n_changed(dlib::int_tag_zero(pdata) + dlib::int_tag_zero(bdata));
    if(pdata) bam_aux_del(p, pdata);
    // If the collapsed observation is now duplex but wasn't before, this updates the DR tag.
    if(pTMP != pFM && pTMP && (pdata = bam_aux_get(p, "DR")) && bam_aux2i(pdata) == 0) {
        pTMP = 1;
        bam_aux_del(p, pdata);
        bam_aux_append(p, "DR", 'i', sizeof(int), (uint8_t *)&pTMP);
    }
    if((pdata = bam_aux_get(p, "NP"))) {
        bdata = bam_aux_get(b, "NP");
        pTMP = bam_aux2i(pdata) + (bdata ? bam_aux2i(bdata) : 1);
        bam_aux_del(p, pdata);
    } else {
        bdata = bam_aux_get(b, "NP");
        pTMP = (bdata ? bam_aux2i(bdata) : 1) + 1;
    }
    bam_aux_append(p, "NP", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&pTMP));
    uint32_t *bPV((uint32_t *)dlib::array_tag(b, "PV")); // Length of this should be b->l_qseq
    uint32_t *pPV((uint32_t *)dlib::array_tag(p, "PV"));
    uint32_t *bFA((uint32_t *)dlib::array_tag(b, "FA"));
    uint32_t *pFA((uint32_t *)dlib::array_tag(p, "FA"));
    uint8_t *bSeq(bam_get_seq(b));
    uint8_t *pSeq(bam_get_seq(p));
    uint8_t *bQual(bam_get_qual(b));
    uint8_t *pQual(bam_get_qual(p));
    const int qlen(p->core.l_qseq);
    int8_t ps, bs;

    if(p->core.flag & (BAM_FREVERSE)) {
        int qleni1;
        for(int i(0); i < qlen; ++i) {
            qleni1 = qlen - i - 1;
            ps = bam_seqi(pSeq, qleni1);
            bs = bam_seqi(bSeq, qleni1);
            if(ps == bs) {
                if((pPV[i] = agreed_pvalues(pPV[i], bPV[i])) < 3) {
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[qleni1] = 2;
                } else {
                    pFA[i] += bFA[i];
                    if(bQual[qleni1] > pQual[qleni1])
                        pQual[qleni1] = bQual[qleni1];
                }
            } else {
                n_base(pSeq, qleni1);
                pFA[i] = 0;
                pPV[i] = 0;
                pQual[qleni1] = 2;
            }
        }
    } else {
        for(int i(0); i < qlen; ++i) {
            ps = bam_seqi(pSeq, i);
            bs = bam_seqi(bSeq, i);
            if(ps == bs) {
                if((pPV[i] = agreed_pvalues(pPV[i], bPV[i])) > 2) {
                    pFA[i] += bFA[i];
                    if(bQual[i] > pQual[i]) pQual[i] = bQual[i];
                } else {
                    n_base(pSeq, i);
                    pPV[i] = 0;
                    pFA[i] = 0;
                }
            } else {
                n_base(pSeq, i);
                pFA[i] = 0;
                pPV[i] = 0;
                pQual[i] = 2;
            }
        }
    }
    bam_aux_append(p, "NC", 'i', sizeof(int), (uint8_t *)&n_changed);
}


inline int string_linear(char *a, char *b, int mmlim)
{
    int hd(0);
    while(*a) if(*a++ != *b++) if(++hd > mmlim) return 0;
    return 1;
}

static const std::vector<uint32_t> ONES(300uL, 1);

inline void add_dummy_tags(bam1_t *b)
{
    static const int one(1);
    int i;
    std::vector<uint32_t> pvbuf;
    pvbuf.reserve(b->core.l_qseq);
    // Set the read to be a singleton
    bam_aux_append(b, "FM", 'i', sizeof(int), const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(&one)));
    // Pass the read
    bam_aux_append(b, "FP", 'i', sizeof(int), const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(&one)));
    uint8_t *qual(bam_get_qual(b));
    if(b->core.flag & BAM_FREVERSE)
        for(i = b->core.l_qseq; i;)
            pvbuf.push_back(static_cast<uint32_t>(qual[--i]));
    else
        for(i = 0; i < b->core.l_qseq; ++i)
            pvbuf.push_back(static_cast<uint32_t>(qual[i]));
    dlib::bam_aux_array_append(b, "FA", 'I', sizeof(uint32_t), b->core.l_qseq, const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(ONES.data())));
    dlib::bam_aux_array_append(b, "PV", 'I', sizeof(uint32_t), b->core.l_qseq, const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(pvbuf.data())));
}


void bam_rsq_bookends(rsq_aux_t *settings)
{
	if(settings->is_se) {
        Stack<same_stack_pos> stack(settings, 1 << 8);
        stack.pe_core(settings);
	} else {
        Stack<same_stack_pos_se> stack(settings, 1 << 8);
        stack.se_core(settings);
    }
}


int rsq_usage(int retcode)
{
    fprintf(stderr,
                    "Positional rescue. \n"
                    "Reads with the same start position are compared.\n"
                    "If their barcodes are sufficiently similar, they are treated as having originated"
                    "from the same original template molecule.\n"
                    "Usage:  bmftools rsq <input.srt.bam> <output.bam>\n\n"
                    "Flags:\n"
                    "-f      Path for the fastq for reads that need to be realigned. REQUIRED.\n"
                    "-s      Flag to write reads with supplementary alignments . Default: False.\n"
                    "-S      Flag to indicate that this rescue is for single-end data.\n"
                    "-t      Mismatch limit. Default: 2\n"
                    "-l      Set bam compression level. Valid: 0-9. (0 == uncompresed)\n"
                    "-m      Trust unmasked bases if reads being collapsed disagree but one is unmasked. Default: mask anyways.\n"
                    "-i      Flag to ignore barcodes and infer solely by positional information. [THIS IS BROKEN]\n"
                    "This flag adds artificial auxiliary tags to treat unbarcoded reads as if they were singletons.\n"
            );
    return retcode;
}


int rsq_main(int argc, char *argv[])
{
    int c;
    char wmode[4]{"wb"};

    rsq_aux_t settings{0};
    settings.mmlim = 2;

    char *fqname(nullptr);

    if(argc < 3) return rsq_usage(EXIT_FAILURE);

    while ((c = getopt(argc, argv, "l:f:t:miSHsh?")) >= 0) {
        switch (c) {
        case 's': settings.write_supp = 1; break;
        case 'S': settings.is_se = 1; break;
        case 'm': settings.trust_unmasked = 1; break;
        case 't': settings.mmlim = atoi(optarg); break;
        case 'f': fqname = optarg; break;
        case 'l': wmode[2] = atoi(optarg)%10 + '0';break;
        case 'i': settings.infer = 1; break;
        case '?': case 'h': case 'H': return rsq_usage(EXIT_SUCCESS);
        }
    }
    if (optind + 2 > argc)
        return rsq_usage(EXIT_FAILURE);

    if(!fqname) {
        fprintf(stderr, "Fastq path for rescued reads required. Abort!\n");
        return rsq_usage(EXIT_FAILURE);
    }

    settings.fqh = fopen(fqname, "w");

    if(!settings.fqh)
        LOG_EXIT("Failed to open output fastq for writing. Abort!\n");


    if(!settings.infer)
        for(const char *tag: {"FM", "FA", "PV", "FP"})
            dlib::check_bam_tag_exit(argv[optind], tag);
    settings.in = sam_open(argv[optind], "r");
    settings.hdr = sam_hdr_read(settings.in);

    if (settings.hdr == nullptr || settings.hdr->n_targets == 0)
        LOG_EXIT("input SAM does not have header. Abort!\n");

    dlib::add_pg_line(settings.hdr, argc, argv, "bmftools rsq", BMF_VERSION, "bmftools",
            "Uses positional information to rescue reads with errors in the barcode.");
    LOG_DEBUG("Write mode: %s.\n", wmode);
    settings.out = sam_open(argv[optind+1], wmode);
    if (settings.in == 0 || settings.out == 0)
        LOG_EXIT("fail to read/write input files\n");
    sam_hdr_write(settings.out, settings.hdr);

    bam_rsq_bookends(&settings);
    bam_hdr_destroy(settings.hdr);
    sam_close(settings.in); sam_close(settings.out);
    if(settings.fqh) fclose(settings.fqh);
    LOG_INFO("Successfully completed bmftools rsq.\n");
    return EXIT_SUCCESS;
}

}
