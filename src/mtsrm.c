/*
 * Multi-threaded sort/rescue/merge
 *
 * Calls bmfsort on an input file or stdin, splits
 */

char **tmpfnames(char *split_prefix, bam_hdr_t *hdr) {
	char **ret = (char **)calloc(hdr->n_targets + 1, sizeof(char *));
	char buf[200];
	int i;
	for(i = 0; i < hdr->n_targets; ++i) {
		sprintf(buf, "%s.%s.bam", split_prefix, hdr->target_name[i]);
		ret[i] = strdup(buf);
	}
	sprintf(buf, "%s.unmapped.bam", split_prefix, hdr->target_name[i]);
	ret[i] = strdup(buf);
	return ret;
}

int main(int argc, char *argv[]) {
	return 0;
}
