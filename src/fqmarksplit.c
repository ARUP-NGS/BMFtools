#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/resource.h>
#include "cstr_util.h"
#include "khash_dmp_core.h"

void print_usage(char *argv[])
{
		fprintf(stderr, "Usage: %s <options> -i <Index.seq> <Fq.R1.seq> <Fq.R2.seq>"
						"\nFlags:\n"
						"-t: Homopolymer failure threshold. A molecular barcode with"
						" a homopolymer of length >= this limit is flagged as QC fail."
						"Default: 10.\n"
						"-o: Output basename. Currently required, as string "
						"manipulation in C is a bit of work and I'd rather spend my "
						"time building code than messing around with string "
						"manipulation that doesn't add to the code base.\n"
						"-i: Index fastq path. Required.\n"
						"-n: Number of nucleotides at the beginning of the barcode to use to split the output.\n"
						"-z: Flag to optionally pipe to gzip while producing final fastqs. Default: False.\n"
						"-g: Gzip compression ratio if piping to gzip (-z). Default: 1 (weak compression).\n"
						"-s: Number of bases from reads 1 and 2 with which to salt the barcode. Default: 0.\n"
						"-m: Number of bases in the start of reads to skip when salting. Default: 0. Recommended: 1.\n"
						"-d: Flag to run hash dmp. Default: False.\n"
						"-p: Number of threads to use if running hash_dmp. Default: 4.\n"
						"-v: Set notification interval for split. Default: 1000000.\n"
						"-c: Flag to optionally cat all files together in one command. Faster than sequential cats, but might break."
						"In addition, won't work for enormous filenames or too many arguments. Default: False.\n"
						"-r: Path to flat text file with rescaled quality scores. If not provided, it will not be used.\n"
						"-f: If running hash_dmp, this sets the Final Fastq Prefix. \n",
						argv[0]);
}

void print_opt_err(char *argv[], char *optarg)
{
	print_usage(argv);
	fprintf(stderr, "Unrecognized option %s. Abort!\n", optarg);
	exit(1);
}

int main(int argc, char *argv[])
{
	// Build settings struct
	const char *default_basename = "metasyntactic_var";
	mss_settings_t settings = {
		.hp_threshold = 10,
		.n_nucs = 2,
		.index_fq_path = NULL,
		.output_basename = (char *)default_basename,
		.input_r1_path = NULL,
		.input_r2_path = NULL,
		.n_handles = 0,
		.notification_interval = 1000000,
		.run_hash_dmp = 0,
		.panthera = 0,
		.gzip_output = 0,
		.ffq_prefix = NULL,
		.salt = 0,
		.offset = 0,
		.threads = 4,
		.gzip_compression = 1,
		.rescaler = NULL,
		.rescaler_path = NULL,
	};
	omp_set_dynamic(0); // Tell omp that I want to set my number of threads 4realz

	int c;
	while ((c = getopt(argc, argv, "t:o:i:n:m:s:f:u:p:g:v:r:hdcz")) > -1) {
		switch(c) {
			case 'c': settings.panthera = 1; break;
			case 'd': settings.run_hash_dmp = 1; break;
			case 'f': settings.ffq_prefix = strdup(optarg); break;
			case 'i': settings.index_fq_path = strdup(optarg); break;
			case 'm': settings.offset = atoi(optarg); break;
			case 'n': settings.n_nucs = atoi(optarg); break;
			case 'o': settings.output_basename = strdup(optarg);break;;
			case 'p': settings.threads = atoi(optarg); break;
			case 's': settings.salt = atoi(optarg); break;
			case 't': settings.hp_threshold = atoi(optarg); break;
			case 'v': settings.notification_interval = atoi(optarg); break;
			case 'z': settings.gzip_output = 1; break;
			case 'g': settings.gzip_compression = atoi(optarg); if(settings.gzip_compression > 9) settings.gzip_compression = 9; break;
			case 'r':
				fprintf(stderr, "About to parse in rescaler.\n");
				settings.rescaler_path = strdup(optarg); settings.rescaler = parse_1d_rescaler(settings.rescaler_path);
				fprintf(stderr, "Parsed rescaler.\n"); break;
			case 'h': print_usage(argv); return 0;
			default: print_opt_err(argv, optarg);
		}
	}


	increase_nofile_limit(settings.threads);
	omp_set_num_threads(settings.threads);

	settings.n_handles = ipow(4, settings.n_nucs);
	if(settings.n_handles * 3 > get_fileno_limit()) {
		int o_fnl = get_fileno_limit();
		increase_nofile_limit(kroundup32(settings.n_handles));
		fprintf(stderr, "Increased nofile limit from %i to %i.\n", o_fnl,
				kroundup32(settings.n_handles));
	}

	if(argc - 1 != optind + 1) {
		fprintf(stderr, "Both read 1 and read 2 fastqs are required. See usage.\n");
		print_usage(argv);
		return 1;
	}
	settings.input_r1_path =  strdup(argv[optind]);
	settings.input_r2_path =  strdup(argv[optind + 1]);

	if(!settings.index_fq_path) {
		fprintf(stderr, "Index fastq required. See usage.\n");
		print_usage(argv);
		return 1;
	}
	if(!settings.output_basename) {
		settings.output_basename = make_crms_outfname(settings.input_r1_path);
		fprintf(stderr, "Output basename not provided. Defaulting to variation on input: %s.\n", settings.output_basename);
	}

/*
	fprintf(stderr, "Hey, can I even read this fastq? %s, %s, %i", seq1->seq.s, seq1->qual.s, l);
	fprintf(stderr, "Hey, my basename is %s\n", settings.output_basename);
*/
	mark_splitter_t *splitter = splitmark_core_rescale(&settings);
	if(settings.run_hash_dmp) {
		fprintf(stderr, "Now executing hash dmp.\n");
		if(!settings.ffq_prefix) {
			settings.ffq_prefix = make_default_outfname(settings.input_r2_path, ".dmp.final");
		}
		// Whatever I end up putting into here.
		splitterhash_params_t *params = init_splitterhash_mss(&settings, splitter);
		for(int i = 0; i < params->n; ++i) {
			fprintf(stderr, "infnames R1 %s, R2 %s. outfnames R1 %s, R2 %s\n",
					params->infnames_r1[i], params->infnames_r2[i],
					params->outfnames_r1[i], params->outfnames_r2[i]);
		}
		fprintf(stderr, "Now running dmp block in parallel with %i threads.\n", settings.threads);
#if NOPARALLEL
#else
		#pragma omp parallel shared(params)
		{
			#pragma omp for
#endif
			for(int i = 0; i < params->n; ++i) {
				char tmpbuf[500];
				fprintf(stderr, "Now running hash dmp core on input filename %s and output filename %s.\n",
						params->infnames_r1[i], params->outfnames_r1[i]);
				khash_dmp_core(params->infnames_r1[i], params->outfnames_r1[i]);
				fprintf(stderr, "Now running hash dmp core on input filename %s and output filename %s.\n",
						params->infnames_r2[i], params->outfnames_r2[i]);
				khash_dmp_core(params->infnames_r2[i], params->outfnames_r2[i]);
				sprintf(tmpbuf, "rm %s %s", params->infnames_r1[i], params->infnames_r2[i]);
				CHECK_CALL(tmpbuf);
			}
#if NOPARALLEL
#else
		}
#endif
		// Remove temporary split files
		fprintf(stderr, "Now removing temporary files.\n");
		/*
		#pragma omp parallel for shared(splitter)
		for(int i = 0; i < splitter->n_handles; ++i) {
			int tmp_ret;
			char tmpbuf[500];
			fprintf(stderr, "Now removing temporary files %s and %s.\n",
					splitter->fnames_r1[i], splitter->fnames_r2[i]);
			sprintf(tmpbuf, "rm %s %s", splitter->fnames_r1[i], splitter->fnames_r2[i]);
			//fprintf(stderr, "Don't feel like executing command '%s' today. Eh.\n", tmpbuf);
			CHECK_CALL(tmpbuf, tmp_ret);
		}
		*/
		char del_buf[500];
		char cat_buff2[CAT_BUFFER_SIZE];
		char cat_buff1[CAT_BUFFER_SIZE];
		// Make sure that both files are empty.
		char cat_buff[CAT_BUFFER_SIZE];
		char ffq_r1[200];
		char ffq_r2[200];
		sprintf(ffq_r1, settings.gzip_output ? "%s.R1.fq.gz": "%s.R1.fq", settings.ffq_prefix);
		sprintf(ffq_r2, settings.gzip_output ? "%s.R2.fq.gz": "%s.R2.fq", settings.ffq_prefix);
		sprintf(cat_buff, "> %s", ffq_r1);
		CHECK_CALL(cat_buff);
		sprintf(cat_buff, "> %s", ffq_r2);
		CHECK_CALL(cat_buff);
		if(!settings.panthera) {
			for(int i = 0; i < settings.n_handles; ++i) {
				// Clear files if present
				if(settings.gzip_output)
					sprintf(cat_buff, "cat %s | gzip - -%i >> %s", params->outfnames_r1[i], settings.gzip_compression, ffq_r1);
				else
					sprintf(cat_buff, "cat %s >> %s", params->outfnames_r1[i], ffq_r1);
				FILE *p1 = popen(cat_buff, "w");

				if(settings.gzip_output)
					sprintf(cat_buff, "cat %s | gzip - -%i >> %s", params->outfnames_r2[i], settings.gzip_compression, ffq_r2);
				else
					sprintf(cat_buff, "cat %s >> %s", params->outfnames_r2[i], ffq_r2);
				CHECK_CALL(cat_buff);
				if(pclose(p1)) {
					fprintf(stderr, "System call failed. Command : '%s'.\n", cat_buff);
					exit(EXIT_FAILURE);
				}
			}
		}
		else {
			fprintf(stderr, "Now building cat string.\n");
			sprintf(cat_buff1, "/bin/cat ");
			sprintf(cat_buff2, "/bin/cat ");
			for(int i = 0; i < settings.n_handles; ++i) {
				strcat(cat_buff1, params->outfnames_r1[i]);
				strcat(cat_buff1, " ");
				strcat(cat_buff2, params->outfnames_r2[i]);
				strcat(cat_buff2, " ");
			}
			if(settings.gzip_output) {
				sprintf(del_buf, " | gzip - -%i", settings.gzip_compression);
				strcat(cat_buff1, del_buf);
				strcat(cat_buff2, del_buf);
			}
			strcat(cat_buff1, " > ");
			strcat(cat_buff1, ffq_r1);
			strcat(cat_buff2, " > ");
			strcat(cat_buff2, ffq_r2);
			FILE *c1_popen = popen(cat_buff1, "w");
			CHECK_CALL(cat_buff2);
			if(pclose(c1_popen)) {
				fprintf(stderr, "First cat command failed. Abort!\n");
				exit(EXIT_FAILURE);
			}
			fprintf(stderr, "Now cleaning up intermediate files.\n");
			#pragma omp parallel for shared(params)
			for(int i = 0; i < params->n; ++i) {
				char tmpbuf[500];
				sprintf(tmpbuf, "rm %s %s", params->outfnames_r1[i], params->outfnames_r2[i]);
				popen(tmpbuf, "w");
			}
			fprintf(stderr, "Finished cleaning up intermediate files.\n");
			//fprintf(stderr, "Not executing %s today. Eh.\n", cat_buff1);
		}
		splitterhash_destroy(params);
	}
	FREE_SPLITTER_PTR(splitter);
	FREE_SETTINGS(settings);
	return 0;
}
