#include "seq_util.h"
#include <math.h>

double UNUSED_FUNC(shannon_entropy) (char *str)
{
	khash_t(shen) *h = kh_init(shen);
	double shen = 0.0;
	khiter_t k = 0;
	int ret = 0;
	uint64_t i;
	for(i = 0; *(str + i); ++i) {
		k = kh_get(shen, h, str[i]);
		if(k == kh_end(h)) {
			k = kh_put(shen, h, str[i], &ret);
			kh_val(h, k) = 1uL;
			continue;
		}
		else
			++kh_val(h, k);
	}
	for(k = 0; k != kh_end(h); ++k) {
		if(!kh_exist(h, k))
			continue;
		const double f = (double)kh_val(h, k) / i;
		shen -= f * log(f);
	}
	kh_destroy(shen, h);
	return shen;
}

double UNUSED_FUNC(shannon_entropy_acgtn) (char *str)
{
	uint64_t counts[5] = {0uL, 0uL, 0uL, 0uL, 0uL};
	double shen = 0.0;
	uint64_t i;
	for(i = 0; *(str + i); ++i) {
		switch(str[i]) {
			case 'a':
			case 'A': ++counts[0]; break;
			case 'c':
			case 'C': ++counts[1]; break;
			case 'g':
			case 'G': ++counts[2]; break;
			case 't':
			case 'T': ++counts[3]; break;
			case 'n':
			case 'N': ++counts[4]; break;
		}
	}
	uint64_t sum = i + 1; // sum is now the total number of read bases
	for(int j = 0; j < 5; ++j) {
		const double f = (double)counts[i] / sum;
		shen -= f * log(f);
	}
	return shen;
}
