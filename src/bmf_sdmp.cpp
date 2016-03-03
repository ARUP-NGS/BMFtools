#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include "bmf_dmp.h"

namespace BMF {

	void sdmp_usage(char *argv[])
	{
	        fprintf(stderr, "Usage: bmftools %s <options> -i <Index.seq> <Fq.R1.seq> <Fq.R2.seq>"
	                        "\nFlags:\n"
	                        "-t: Homopolymer failure threshold. A molecular barcode with"
	                        " a homopolymer of length >= this limit is flagged as QC fail."
	                        "Default: 10.\n"
	                        "-o: Temporary fastq file prefix.\n"
	                        "-i: Index fastq path. Required.\n"
	                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: %i.\n"
	                        "-z: Flag to optionally pipe to gzip while producing final fastqs. Default: False.\n"
	                        "-T: If unset, write uncompressed plain text temporary files. If not, use that compression level for temporary files.\n"
	                        "-g: Gzip compression ratio if writing gzipped. Default (if writing compressed): 1 (mostly to reduce I/O).\n"
	                        "-s: Number of bases from reads 1 and 2 with which to salt the barcode. Default: 0.\n"
	                        "-m: Number of bases in the start of reads to skip when salting. Default: 1.\n"
	                        "-d: Flag to run hash dmp. Default: False.\n"
	                        "-p: Number of threads to use if running hash_dmp. Default: %i.\n"
	                        "-v: Set notification interval for split. Default: 1000000.\n"
	                        "-c: Flag to optionally cat all files together in one command. Faster than sequential cats, but might break."
	                        "In addition, won't work for enormous filenames or too many arguments. Default: False.\n"
	                        "-r: Path to flat text file with rescaled quality scores. If not provided, it will not be used.\n"
	                        "-w: Flag to leave temporary files instead of deleting them, as in default behavior.\n"
	                        "-f: If running hash_dmp, this sets the Final Fastq Prefix. \n"
	                        "-S: Single-end mode. Ignores read 2.\n"
	                        "-O: Emit final fastqs to stdout in interleaved form. Ignores -f.\n",
	                        argv[0], DEFAULT_N_NUCS, DEFAULT_N_THREADS);
	}

	static mark_splitter_t *splitmark_core_rescale(marksplit_settings_t *settings)
	{
	    if(strcmp(settings->input_r1_path, settings->input_r2_path) == 0) {
	        fprintf(stderr, "[E:%s]Input read paths are the same {'R1': %s, 'R2': %s}. WTF!\n",
	                __func__, settings->input_r1_path, settings->input_r2_path);
	        exit(EXIT_FAILURE);
	    }
	    else
	        fprintf(stderr, "[%s] Path to index fq: %s.\n", __func__, settings->index_fq_path);
	    gzFile fp_read1, fp_read2, fp_index;
	    kseq_t *seq1 = NULL, *seq2 = NULL, *seq_index = NULL;
	    int l1, l2, l_index;
	    mark_splitter_t *splitter_ptr = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
	    *splitter_ptr = init_splitter(settings);
	    if(!isfile(settings->input_r1_path) ||
	       !isfile(settings->input_r2_path) ||
	       !isfile(settings->index_fq_path)) {
	        fprintf(stderr, "[E:%s] At least one input path ('%s', '%s', '%s') is not a file. Abort!\n", __func__,
	                settings->input_r1_path, settings->input_r2_path, settings->index_fq_path);
	        exit(EXIT_FAILURE);
	    }
	    // Open fastqs
	    fprintf(stderr, "[%s] Splitter now opening files R1 ('%s'), R2 ('%s'), index ('%s').\n",
	            __func__, settings->input_r1_path, settings->input_r2_path, settings->index_fq_path);
	    fp_read1 = gzopen(settings->input_r1_path, "r"), fp_read2 = gzopen(settings->input_r2_path, "r");
	    seq1 = kseq_init(fp_read1), seq2 = kseq_init(fp_read2);
	    l1 = kseq_read(seq1), l2 = kseq_read(seq2);

	    fp_index = gzopen(settings->index_fq_path, "r");
	    seq_index = kseq_init(fp_index),
	    l_index = kseq_read(seq_index);

	    uint64_t bin = 0;
	    int pass_fail = 1;
	    tmp_mseq_t *tmp = init_tm_ptr(seq1->seq.l, seq_index->seq.l + 2 * settings->salt);
	#if !NDEBUG
	    fprintf(stderr, "[D:%s] About to check for failed opening.\n", __func__);
	#endif
	    if(l1 < 0 || l2 < 0 || l_index < 0) {
	        fprintf(stderr, "[E:%s] Could not read input fastqs. Abort mission!\n", __func__);
	        exit(EXIT_FAILURE);
	    }
	    check_rescaler(settings, nqscores * 2 * 4 * seq1->seq.l);
	    mseq_t *rseq1 = mseq_init(seq1, settings->rescaler, 0); // rseq1 is initialized
	    mseq_t *rseq2 = mseq_init(seq2, settings->rescaler, 1); // rseq2 is initialized
	    memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
	    memcpy(rseq1->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
	    memcpy(rseq1->barcode + settings->salt + seq_index->seq.l, seq2->seq.s + settings->offset, settings->salt);
	    rseq1->barcode[settings->salt * 2 + seq_index->seq.l] = '\0';
	    update_mseq(rseq1, seq1, settings->rescaler, tmp, 0, 0);
	    update_mseq(rseq2, seq2, settings->rescaler, tmp, 0, 1);
	    pass_fail = test_hp(rseq1->barcode, settings->hp_threshold);
	    bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
	    mseq2fq(splitter_ptr->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode);
	    mseq2fq(splitter_ptr->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode);
	    uint64_t count = 0;
	#if !NDEBUG
	    fprintf(stderr, "[D:%s] About to start looping.\n", __func__);
	#endif
	    while(LIKELY((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2) >= 0))
	            && (l_index = kseq_read(seq_index)) >= 0) {
	        if(UNLIKELY(++count % settings->notification_interval == 0))
	            fprintf(stderr, "[%s] Number of records processed: %lu.\n", __func__, count);
	        memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
	        memcpy(rseq1->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
	        memcpy(rseq1->barcode + settings->salt + seq_index->seq.l, seq2->seq.s + settings->offset, settings->salt);
	        update_mseq(rseq1, seq1, settings->rescaler, tmp, 0, 0);
	        update_mseq(rseq2, seq2, settings->rescaler, tmp, 0, 1);
	        pass_fail = test_hp(rseq1->barcode, settings->hp_threshold);
	        bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
	        mseq2fq(splitter_ptr->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode);
	        mseq2fq(splitter_ptr->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode);
	    }
	    tm_destroy(tmp);
	    mseq_destroy(rseq1); mseq_destroy(rseq2);
	    kseq_destroy(seq1); kseq_destroy(seq2); kseq_destroy(seq_index);
	    gzclose(fp_read1); gzclose(fp_read2); gzclose(fp_index);
	    for(int j = 0; j < settings->n_handles; ++j) {
	        gzclose(splitter_ptr->tmp_out_handles_r1[j]);
	        gzclose(splitter_ptr->tmp_out_handles_r2[j]);
	        splitter_ptr->tmp_out_handles_r1[j] = splitter_ptr->tmp_out_handles_r2[j] = NULL;
	    }
	    return splitter_ptr;
	}

	static mark_splitter_t *splitmark_core_rescale_se(marksplit_settings_t *settings)
	{
	#if !NDEBUG
	    fprintf(stderr, "[D:%s] Path to index fq: %s.\n", __func__, settings->index_fq_path);
	#endif
	    gzFile fp, fp_index;
	    kseq_t *seq = NULL, *seq_index = NULL;
	    int l, l_index;
	    mark_splitter_t *splitter_ptr = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
	    *splitter_ptr = init_splitter(settings);
	    if(!isfile(settings->input_r1_path) ||
	       !isfile(settings->index_fq_path)) {
	        fprintf(stderr, "[E:%s] At least one input path ('%s', '%s') is not a file. Abort!\n", __func__,
	                settings->input_r1_path, settings->index_fq_path);
	        exit(EXIT_FAILURE);
	    }
	    // Open fastqs
	    fp = gzopen(settings->input_r1_path, "r"), fp_index = gzopen(settings->index_fq_path, "r");
	    seq = kseq_init(fp), seq_index = kseq_init(fp_index);
	    l = kseq_read(seq), l_index = kseq_read(seq_index);
	    seq_index = kseq_init(fp_index),
	    l_index = kseq_read(seq_index);

	    tmp_mseq_t *tmp = init_tm_ptr(seq->seq.l, seq_index->seq.l + settings->salt);
	    if(l < 0 || l_index < 0) {
	        fprintf(stderr, "[E:%s] Could not read input fastqs. Abort mission!\n", __func__);
	        exit(EXIT_FAILURE);
	    }
	#if !NDEBUG
	    fprintf(stderr, "[D:%s] Splitter now opening files read ('%s') and index ('%s').\n",
	            __func__, settings->input_r1_path, settings->index_fq_path);
	#endif
	    mseq_t *rseq = mseq_init(seq, settings->rescaler, 0);
	    memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
	    memcpy(rseq->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
	    rseq->barcode[settings->salt + seq_index->seq.l] = '\0';
	    update_mseq(rseq, seq, settings->rescaler, tmp, 0, 0);
	    mseq2fq(splitter_ptr->tmp_out_handles_r1[get_binner_type(rseq->barcode, settings->n_nucs, uint64_t)],
	            rseq, test_hp(rseq->barcode, settings->hp_threshold), rseq->barcode);
	    uint64_t count = 0;
	    while (LIKELY((l = kseq_read(seq)) >= 0 && (l_index = kseq_read(seq_index)) >= 0)) {
	        if(UNLIKELY(++count % settings->notification_interval == 0))
	            fprintf(stderr, "[%s] Number of records processed: %lu.\n", __func__, count);
	        memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
	        memcpy(rseq->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
	        update_mseq(rseq, seq, settings->rescaler, tmp, 0, 0);
	        mseq2fq(splitter_ptr->tmp_out_handles_r1[get_binner_type(rseq->barcode, settings->n_nucs, uint64_t)],
	                rseq, test_hp(rseq->barcode, settings->hp_threshold), rseq->barcode);
	    }
	    tm_destroy(tmp);
	    mseq_destroy(rseq);
	    kseq_destroy(seq); kseq_destroy(seq_index);
	    gzclose(fp); gzclose(fp_index);
	    for(int j = 0; j < settings->n_handles; ++j)
	        gzclose(splitter_ptr->tmp_out_handles_r1[j]);
	    // Set out handles to NULL.
	    memset(splitter_ptr->tmp_out_handles_r1, 0, settings->n_handles * sizeof(FILE *));
	    return splitter_ptr;
	}

	int sdmp_main(int argc, char *argv[])
	{
	    if(argc < 3 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
	        sdmp_usage(argv); return EXIT_FAILURE;
	    }
	    // Build settings struct
	    marksplit_settings_t settings = {0};
	    settings.hp_threshold = 10;
	    settings.n_nucs = DEFAULT_N_NUCS;
	    settings.notification_interval = 1000000;
	    settings.threads = DEFAULT_N_THREADS;
	    settings.gzip_compression = 1;
	    settings.cleanup = 1;
	    sprintf(settings.mode, "wT");

	    int c;
	    while ((c = getopt(argc, argv, "t:o:i:n:m:s:f:u:p:g:v:r:T:hdczw?S&")) > -1) {
	        switch(c) {
	            case 'd': settings.run_hash_dmp = 1; break;
	            case 'f': settings.ffq_prefix = strdup(optarg); break;
	            case 'i': settings.index_fq_path = strdup(optarg); break;
	            case 'm': settings.offset = atoi(optarg); break;
	            case 'n': settings.n_nucs = atoi(optarg); break;
	            case 'o': settings.tmp_basename = strdup(optarg);break;
	            case 'T': sprintf(settings.mode, "wb%i", atoi(optarg) % 10); break;
	            case 'p': settings.threads = atoi(optarg); break;
	            case 's': settings.salt = atoi(optarg); break;
	            case 't': settings.hp_threshold = atoi(optarg); break;
	            case 'v': settings.notification_interval = atoi(optarg); break;
	            case 'z': settings.gzip_output = 1; break;
	            case 'g': settings.gzip_compression = atoi(optarg); if(settings.gzip_compression > 9) settings.gzip_compression = 9; break;
	            case 'w': settings.cleanup = 0; break;
	            case 'r':
	                settings.rescaler_path = strdup(optarg); settings.rescaler = parse_1d_rescaler(settings.rescaler_path);
	                break;
	            case 'S': settings.is_se = 1; break;
	            case 'O': settings.to_stdout = 1; break;
	            case '?': case 'h': sdmp_usage(argv); return EXIT_SUCCESS;
	        }
	    }

        dlib::increase_nofile_limit(settings.threads);
	    omp_set_num_threads(settings.threads);

	    settings.n_handles = ipow(4, settings.n_nucs);
	    if(settings.n_handles * 3 > dlib::get_fileno_limit()) {
	        int o_fnl = dlib::get_fileno_limit();
	        dlib::increase_nofile_limit(kroundup32(settings.n_handles));
	        fprintf(stderr, "Increased nofile limit from %i to %i.\n", o_fnl,
	                kroundup32(settings.n_handles));
	    }

	    if(argc == 1) {
	        sdmp_usage(argv);
	        return EXIT_SUCCESS;
	    }

	    if(settings.is_se) {
	        if(argc != optind + 1) {
	            fprintf(stderr, "[E:%s] Precisely one input fastq required for se mode. See usage.\n", __func__);
	            sdmp_usage(argv);
	            return EXIT_FAILURE;
	        }
	        settings.input_r1_path =  strdup(argv[optind]);
	    }
	    else {
	        if(argc != optind + 2) {
	            fprintf(stderr, "[E:%s] Both read 1 and read 2 fastqs are required. See usage.\n", __func__);
	            sdmp_usage(argv);
	            return EXIT_FAILURE;
	        }
	        settings.input_r1_path =  strdup(argv[optind]);
	        settings.input_r2_path =  strdup(argv[optind + 1]);
	    }

	    if(!settings.index_fq_path) {
	        fprintf(stderr, "[E:%s] Index fastq required. See usage.\n", __func__);
	        sdmp_usage(argv);
	        return EXIT_FAILURE;
	    }
	    if(!settings.tmp_basename) {
	        settings.tmp_basename = (char *)malloc(21);
	        rand_string(settings.tmp_basename, 20);
	        fprintf(stderr, "[%s] Mark/split prefix not provided. Defaulting to random string ('%s').\n",
	                __func__, settings.tmp_basename);
	    }

	    splitterhash_params_t *params = NULL;
	    mark_splitter_t *splitter = settings.is_se ? splitmark_core_rescale_se(&settings): splitmark_core_rescale(&settings);
	    if(!settings.run_hash_dmp) {
	        fprintf(stderr, "[%s] Finished mark/split.\n", __func__);
	        goto cleanup;
	    }
	    fprintf(stderr, "[%s] Now executing hashmap-powered read collapsing and molecular demultiplexing.\n",
	                    __func__);
	    if(!settings.ffq_prefix) make_outfname(&settings);
	    params = init_splitterhash(&settings, splitter);
	    fprintf(stderr, "[%s] Running dmp block in parallel with %i threads.\n", __func__, settings.threads);

	    parallel_hash_dmp_core(&settings, params, &hash_dmp_core);
	    // Make sure that both files are empty.
	    char ffq_r1[200], ffq_r2[200];
	    sprintf(ffq_r1, settings.gzip_output ? "%s.R1.fq": "%s.R1.fq.gz", settings.ffq_prefix);
	    sprintf(ffq_r2, settings.gzip_output ? "%s.R2.fq": "%s.R2.fq.gz", settings.ffq_prefix);
	    if(settings.to_stdout)
	        call_stdout(&settings, params, ffq_r1, ffq_r2);
	    else cat_fastqs(&settings, params, ffq_r1, ffq_r2);
	    cleanup_hashdmp(&settings, params);
	    splitterhash_destroy(params);
	    fprintf(stderr, "[%s] Successfully completed bmftools sdmp.\n", __func__);

	    cleanup:
	    splitter_destroy(splitter);
	    free_marksplit_settings(settings);
	    return EXIT_SUCCESS;
	}

}
