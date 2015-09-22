#ifndef MAX_N_BLENS
#define MAX_N_BLENS 5
#endif
#ifndef MAX_HOMING_SEQUENCE
#define MAX_HOMING_SEQUENCE 8
#endif

typedef struct blens {
    int max_blen; // Last value in blens
    int min_blen; // Lowest value in blens
    int blens[MAX_N_BLENS]; // Array holding blens
    int n; // Number of blens to look for
    int current_blen;
    int homing_sequence_length;
    char homing_sequence[MAX_HOMING_SEQUENCE + 1];
} blens_t;

typedef struct crms_settings {
    int hp_threshold; // The minimum length of a homopolymer run to fail a barcode.
    int n_nucs; // Number of nucleotides to split by.
    char *output_basename;
    char *input_r1_path;
    char *input_r2_path;
    int n_handles; // Number of handles
    int notification_interval; // How many sets of records do you want to process between progress reports?
    blens_t *blen_data;
    int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
    char ****rescaler; // Three-dimensional rescaler array. Size: [readlen, 39, 4] (length of reads, number of original quality scores, number of bases)
    char *rescaler_path; // Path to rescaler for
} crms_settings_t;

void free_rescaler_array(crms_settings_t settings) {
    int readlen = count_lines(settings.rescaler_path);
    for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < readlen; ++j) {
            for(int k = 0; k < 39; ++k) {
                if(settings.rescaler[i][j][k]) {
                   free(settings.rescaler[i][j][k]);
                   settings.rescaler[i][j][k] = NULL;
                }
            }
            if(settings.rescaler[i][j]) {
                free(settings.rescaler[i][j]);
                settings.rescaler[i][j] = NULL;
            }
        }
        if(settings.rescaler[i]) {
            free(settings.rescaler[i]);
            settings.rescaler[i] = NULL;
        }
    }
    free(settings.rescaler);
    settings.rescaler = NULL;
    return;
}
