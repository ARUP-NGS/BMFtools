#ifndef SPLITTER_H
#define SPLITTER_H
#include <stdint.h>
#include <zlib.h>

namespace BMF {

    struct marksplit_settings_t {
        uint32_t blen:16;
        uint32_t blen1_2:16;
        char *ffq_prefix; // Final fastq prefix
        char *homing_sequence; // Homing sequence...
        int homing_sequence_length; // Length of homing sequence, should it be used.
        char *input_r1_path;
        char *input_r2_path;
        char *index_fq_path; // Make sure this is null if it's inline!
        int max_blen;
        int n_handles; // Number of handles
        int notification_interval; // How many sets of records do you want to process between progress reports?
        uint32_t n_nucs:10;
        uint32_t offset:4;
        uint32_t salt:4;
        uint32_t run_hash_dmp:1;
        uint32_t is_se:1;
        uint32_t cleanup:1; // Set to false to leave temporary files
        uint32_t to_stdout:1;
        uint32_t gzip_output:1;
        uint32_t gzip_compression:4;
        uint32_t hp_threshold:5;
        char *tmp_basename;
        char *rescaler; // Four-dimensional rescaler array. Size: [readlen, NQSCORES, 4] (length of reads, number of original quality scores, number of bases)
        char *rescaler_path; // Path to rescaler for
        int threads;
        char mode[4];
    };

    void free_marksplit_settings(marksplit_settings_t settings);

    struct mark_splitter_t {
        gzFile *tmp_out_handles_r1;
        gzFile *tmp_out_handles_r2;
        uint32_t n_nucs;
        int n_handles;
        char **fnames_r1;
        char **fnames_r2;
    };

    mark_splitter_t init_splitter(marksplit_settings_t* settings_ptr);
    void splitter_destroy(mark_splitter_t *var);

    struct splitterhash_params_t {
        char **infnames_r1;
        char **infnames_r2;
        char **outfnames_r1;
        char **outfnames_r2;
        int n; // Number of infnames and outfnames
        int paired; // 1 if paired, 0 if single-end
    };

    void splitterhash_destroy(splitterhash_params_t *params);
    splitterhash_params_t *init_splitterhash(marksplit_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);

} /* namespace BMF */

#endif /* SPLITTER_H */
