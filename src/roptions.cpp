#include "roptions.h"
#include "rutils.h"
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

static ko_longopt_t long_options[] = {
	{ (char*)"level_column",     	ko_required_argument, 	300 },
	{ (char*)"q-mid-occ",     		ko_required_argument, 	301 },
	{ (char*)"mid_occ_frac",     	ko_required_argument, 	302 },
	{ (char*)"min-events",			ko_required_argument, 	303 },
	{ (char*)"bw",					ko_required_argument, 	304 },
	{ (char*)"max-target-gap",		ko_required_argument, 	305 },
	{ (char*)"max-query-gap",		ko_required_argument, 	306 },
	{ (char*)"min-anchors",			ko_required_argument, 	307 },
	{ (char*)"min-score",			ko_required_argument, 	308 },
	{ (char*)"chain-gap-scale",		ko_required_argument, 	309 },
	{ (char*)"chain-skip-scale",	ko_required_argument, 	310 },
	{ (char*)"best-chains",			ko_required_argument, 	311 },
	{ (char*)"primary-ratio",		ko_required_argument, 	312 },
	{ (char*)"primary-length",		ko_required_argument, 	313 },
	{ (char*)"max-skips",			ko_required_argument, 	314 },
	{ (char*)"max-iterations",		ko_required_argument, 	315 },
	{ (char*)"rmq",					ko_no_argument, 	  	316 },
	{ (char*)"rmq-inner-dist",		ko_required_argument, 	317 },
	{ (char*)"rmq-size-cap",		ko_required_argument, 	318 },
	{ (char*)"bw-long",				ko_required_argument, 	319 },
	{ (char*)"max-chunks",			ko_required_argument, 	320 },
	{ (char*)"min-mapq",			ko_required_argument, 	321 },
	{ (char*)"alt-drop",			ko_required_argument, 	322 },
	{ (char*)"w-besta",				ko_required_argument, 	323 },
	{ (char*)"w-bestma",			ko_required_argument, 	324 },
	{ (char*)"w-bestq",				ko_required_argument, 	325 },
	{ (char*)"w-bestmq",			ko_required_argument, 	326 },
	{ (char*)"w-bestmc",			ko_required_argument, 	327 },
	{ (char*)"w-threshold",			ko_required_argument, 	328 },
	{ (char*)"bp-per-sec",			ko_required_argument, 	329 },
	{ (char*)"sample-rate",			ko_required_argument, 	330 },
	{ (char*)"chunk-size",			ko_required_argument, 	331 },
	{ (char*)"seg-window-length1",	ko_required_argument, 	332 },
	{ (char*)"seg-window-length2",	ko_required_argument, 	333 },
	{ (char*)"seg-threshold1",		ko_required_argument, 	334 },
	{ (char*)"seg-threshold2",		ko_required_argument, 	335 },
	{ (char*)"seg-peak-height",		ko_required_argument, 	336 },
	{ (char*)"sequence-until",     	ko_no_argument,       	337 },
	{ (char*)"threshold",			ko_required_argument, 	338 },
	{ (char*)"n-samples",			ko_required_argument, 	339 },
	{ (char*)"test-frequency",		ko_required_argument, 	340 },
	{ (char*)"min-reads",			ko_required_argument, 	341 },
	{ (char*)"occ-frac",			ko_required_argument, 	342 },
	{ (char*)"depletion",			ko_no_argument, 	  	343 },
	{ (char*)"store-sig",			ko_no_argument, 	  	344 },
	{ (char*)"sig-target",			ko_no_argument, 	  	345 },
	{ (char*)"disable-adaptive",	ko_no_argument, 	  	346 },
	{ (char*)"sig-diff",			ko_required_argument, 	347 },
	{ (char*)"align",				ko_no_argument, 	  	348 },
	{ (char*)"dtw-evaluate-chains",	ko_no_argument,		  	349 },
	{ (char*)"dtw-output-cigar",	ko_no_argument,		  	350 },
	{ (char*)"dtw-border-constraint", ko_required_argument,	351 },
	{ (char*)"dtw-log-scores",		ko_no_argument,			352 },
	{ (char*)"no-chainingscore-filtering",	ko_no_argument,	353 },
	{ (char*)"dtw-match-bonus",		ko_required_argument,	354 },
	{ (char*)"output-chains",		ko_no_argument,			355 },
	{ (char*)"dtw-fill-method",		ko_required_argument,	356 },
	{ (char*)"dtw-min-score", 		ko_required_argument, 	357 },
	{ (char*)"log-anchors",			ko_no_argument,			358 },
	{ (char*)"log-num-anchors",		ko_no_argument,			359 },
	{ (char*)"rev-collision-count", ko_required_argument, 	360 },
	{ (char*)"chn-rev-bump", 		ko_required_argument, 	361 },
	{ (char*)"rev-query",			ko_no_argument, 		362 },
	{ (char*)"r10",					ko_no_argument, 		363 },
	{ (char*)"fine-min",			ko_required_argument, 	364 },
	{ (char*)"fine-max",			ko_required_argument, 	365 },
	{ (char*)"fine-range",			ko_required_argument, 	366 },
	{ (char*)"out-quantize",		ko_no_argument,  		367 },
	{ (char*)"no-event-detection",	ko_no_argument,  		368 },
	{ (char*)"io-thread",			ko_required_argument, 	369 },
	{ (char*)"min-score2",			ko_required_argument, 	370 },
	{ (char*)"version",				ko_no_argument, 	  	370 },
	{ 0, 0, 0 }
};

int ri_set_opt(const char *preset, ri_idxopt_t *io, ri_mapopt_t *mo)
{
	if (preset == 0) {
		ri_idxopt_init(io);
		ri_mapopt_init(mo);
	} else if (strcmp(preset, "viral") == 0) {
		io->e = 6;
		mo->bw = 100; mo->max_target_gap_length = 500; mo->max_query_gap_length = 500;
		mo->max_num_chunk = 5, mo->min_chaining_score = 10; mo->chain_gap_scale = 1.2f; mo->chain_skip_scale = 0.3f;
	} else if (strcmp(preset, "sensitive") == 0) {
		//default
	} else if (strcmp(preset, "fast") == 0) {
		io->fine_range = 0.6;
		mo->min_mapq = 5, mo->min_chaining_score = 10, mo->chain_gap_scale = 0.6f;
	} else if (strcmp(preset, "faster") == 0) {
		io->e = 11; io->w = 3;
		io->fine_range = 0.6;
		mo->max_num_chunk = 5; mo->min_mapq = 5, mo->min_chaining_score = 10, mo->chain_gap_scale = 0.6f;
	} else if (strcmp(preset, "ava-viral") == 0) {
		io->e = 6;
		mo->chain_gap_scale = 1.2f; mo->chain_skip_scale = 0.3f;

		io->w = 0;
		io->diff = 0.45f;
		mo->min_chaining_score = 20;
		mo->min_chaining_score2 = 30;
		mo->min_num_anchors = 5;
		mo->min_mapq = 5;
		mo->bw = 1000;
		mo->max_target_gap_length = 2500;
		mo->max_query_gap_length = 2500;

		io->flag |= RI_I_SIG_TARGET;
		mo->flag |= RI_M_ALL_CHAINS;
		mo->flag |= RI_M_NO_ADAPTIVE;

		mo->pri_ratio = 0.0f;
	} else if (strcmp(preset, "ava") == 0) {
		//default
		io->w = 3;
		io->diff = 0.45f;
		mo->min_chaining_score = 40;
		mo->min_chaining_score2 = 75;
		mo->min_num_anchors = 5;
		mo->min_mapq = 5;
		mo->bw = 5000;
		mo->max_target_gap_length = 2500;
		mo->max_query_gap_length = 2500;

		// mo->min_mid_occ = 5000;

		io->flag |= RI_I_SIG_TARGET;
		mo->flag |= RI_M_ALL_CHAINS;
		mo->flag |= RI_M_NO_ADAPTIVE;

		mo->pri_ratio = 0.0f;
	} else if (strcmp(preset, "ava-sensitive") == 0) {
		//default
		io->w = 0;
		io->diff = 0.45f;
		mo->min_chaining_score = 75;
		mo->min_chaining_score2 = 100;
		mo->min_num_anchors = 5;
		mo->min_mapq = 5;
		mo->bw = 1000;
		mo->max_target_gap_length = 2500;
		mo->max_query_gap_length = 2500;

		// mo->min_mid_occ = 10000;

		io->flag |= RI_I_SIG_TARGET;
		mo->flag |= RI_M_ALL_CHAINS;
		mo->flag |= RI_M_NO_ADAPTIVE;

		mo->pri_ratio = 0.0f;
	} else if (strcmp(preset, "ava-large") == 0) {
		io->fine_range = 0.6;
		mo->chain_gap_scale = 0.6f;

		io->w = 5;
		io->diff = 0.45f;
		mo->min_chaining_score = 20;
		mo->min_chaining_score2 = 50;
		mo->min_num_anchors = 2;
		mo->min_mapq = 2;
		mo->bw = 5000;
		mo->max_target_gap_length = 2500;
		mo->max_query_gap_length = 2500;

		// mo->min_mid_occ = 10000;

		io->flag |= RI_I_SIG_TARGET;
		mo->flag |= RI_M_ALL_CHAINS;
		mo->flag |= RI_M_NO_ADAPTIVE;

		mo->pri_ratio = 0.0f;
	} else if (strcmp(preset, "sequence-until") == 0) {
		//default
	} else return -1;
	return 0;
}

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	return (int64_t)(x + .499);
}

static inline void yes_or_no(ri_mapopt_t *opt, int64_t flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int ri_mapopt_parse_dtw_border_constraint(ri_mapopt_t *opt, char* arg){
	if(strcmp(arg, "global") == 0){
		opt->dtw_border_constraint = RI_M_DTW_BORDER_CONSTRAINT_GLOBAL;
	}
	else if(strcmp(arg, "sparse") == 0){
		opt->dtw_border_constraint = RI_M_DTW_BORDER_CONSTRAINT_SPARSE;
	}
	else if(strcmp(arg, "local") == 0){
		opt->dtw_border_constraint = RI_M_DTW_BORDER_CONSTRAINT_LOCAL;
	}
	else{
		return -1;
	}
	return 0;
}

int ri_mapopt_parse_dtw_fill_method(ri_mapopt_t *opt, char* arg) {
	if (strcmp(arg, "banded") == 0) {
		opt->dtw_fill_method = RI_M_DTW_FILL_METHOD_BANDED;
	} else if (strcmp(arg, "full") == 0) {
		opt->dtw_fill_method = RI_M_DTW_FILL_METHOD_FULL;
	} else if (strncmp(arg, "banded=", 7) == 0) {
		opt->dtw_fill_method = RI_M_DTW_FILL_METHOD_BANDED;
		opt->dtw_band_radius_frac = std::atof(arg + 7);
	} else {
		return -1;
	}
	return 0;
}

const char* ri_maptopt_dtw_mode_to_string(uint32_t dtw_border_constraint){
	switch(dtw_border_constraint){
	case RI_M_DTW_BORDER_CONSTRAINT_GLOBAL:
		return "full";
	case RI_M_DTW_BORDER_CONSTRAINT_SPARSE:
		return "sparse";
	case RI_M_DTW_BORDER_CONSTRAINT_LOCAL:
		return "window";
	default:
		return "unknown";
	}
}

void ri_idxopt_init(ri_idxopt_t *opt)
{
	memset(opt, 0, sizeof(ri_idxopt_t));
	opt->e = 8; opt->w = 0; opt->q = 4; opt->n = 0; opt->k = 6, opt->lev_col = 1;
	opt->b = 14;
	opt->diff = 0.35f;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;

	opt->fine_min = -2.0f;
	opt->fine_max = 2.0f;
	opt->fine_range = 0.4;

	opt->window_length1 = 3; //--seg-window-length1
    opt->window_length2 = 9; //--seg-window-length2
    opt->threshold1 = 4.0f; //--seg-threshold1
    opt->threshold2 = 3.5f; //--seg-threshold2
    opt->peak_height = 0.4f; //--seg-peak-height

	// opt->window_length1 = 3; //--seg-window-length1
    // opt->window_length2 = 6; //--seg-window-length2
    // opt->threshold1 = 1.4f; //--seg-threshold1
    // opt->threshold2 = 9.0f; //--seg-threshold2
    // opt->peak_height = 0.2f; //--seg-peak-height

	opt->bp_per_sec = 450; //--bp-per-sec
	opt->sample_rate = 4000; //--sample-rate
	opt->sample_per_base = (float)opt->sample_rate / opt->bp_per_sec;
}

void ri_mapopt_init(ri_mapopt_t *opt)
{
	memset(opt, 0, sizeof(ri_mapopt_t));

	opt->bp_per_sec = 450; //--bp-per-sec
	opt->sample_rate = 4000; //--sample-rate
	opt->chunk_size = 4000; //--chunk-size
	opt->sample_per_base = (float)opt->sample_rate / opt->bp_per_sec;

	//seeding
	// opt->mid_occ_frac = 75e-4f; //--mid-occ-frac
	opt->mid_occ_frac = 1e-2f; //--mid-occ-frac
	// opt->mid_occ_frac = 5e-3f; //--mid-occ-frac
	opt->q_occ_frac = 1e-2f; //--q-occ-frac
	opt->min_mid_occ = 50; //--q-mid-occ [I1, I2]
	opt->max_mid_occ = 500000; //--q-mid-occ [I1, I2]

	opt->max_max_occ = 32767;
	opt->occ_dist = 500;

	//chaining
	opt->bw = 500; //--bw
	opt->bw_long = 0; //--bw-long
	opt->max_target_gap_length = 2500; //--max-target-gap
	opt->max_query_gap_length = 2500; //--max-query-gap
	opt->max_chain_iter = 200; //--max-iterations
	opt->max_num_skips = 5; //--max-skips
	opt->min_num_anchors = 2; //--min-anchors
	// opt->num_best_chains = 3; //--best-chains
	opt->min_chaining_score = 15; //--min-score
	opt->min_chaining_score2 = 0; //--min-score2
	opt->rmq_inner_dist = 1000; //--rmq-inner-dist
	opt->rmq_size_cap = 100000; //--rmq-size-cap
	opt->chain_gap_scale = 0.8f; //--chain-gap-scale
	opt->chain_skip_scale = 0.0f; //--chain-skip-scale

	opt->mask_level = 0.5f; //--primary-ratio
	opt->mask_len = INT_MAX; //--primary-length
	
	opt->pri_ratio = 0.3f;
	opt->best_n = 0; //--best-chains

	opt->top_n_mean = 0; //--top-n-mean

	opt->alt_drop = 0.15f; //--alt-drop

	opt->w_bestmq=0.05f; //--w-bestmq
	opt->w_bestmc=0.6f; //--w-bestmc
	opt->w_bestq=0.35f; //--w-bestq
	opt->w_besta=0.2f; //--w-besta
	opt->w_bestma=0.2f; //--w-bestma
	// opt->w_best2c=0.1f; //--w-best2c
	opt->w_threshold = 0.45f; //--w-threshold

	opt->mini_batch_size = 500000000; //-K

	opt->step_size = 1;
	opt->min_events = 50; //--min-events
	opt->max_num_chunk = 10;//--max-chunks

	opt->min_mapq = 2; //--min-mapq

	//dtw
	opt->dtw_border_constraint = RI_M_DTW_BORDER_CONSTRAINT_SPARSE;
	opt->dtw_fill_method = RI_M_DTW_FILL_METHOD_BANDED;
	opt->dtw_band_radius_frac = 0.10f;
	opt->dtw_match_bonus = 0.4f;
	opt->dtw_min_score = 20.0f;

	//Reverse complementing
	opt->rev_col_limit = 100;
	opt->chn_rev_bump = 1.0f;

	//Default options for event detection.
	opt->window_length1 = 3; //--seg-window-length1
    opt->window_length2 = 9; //--seg-window-length2
    opt->threshold1 = 4.0f; //--seg-threshold1
    opt->threshold2 = 3.5f; //--seg-threshold2
    opt->peak_height = 0.4f; //--seg-peak-height

	// opt->window_length1 = 3; //--seg-window-length1
    // opt->window_length2 = 7; //--seg-window-length2
    // opt->threshold1 = 4.0f; //--seg-threshold1
    // opt->threshold2 = 3.0f; //--seg-threshold2
    // opt->peak_height = 0.4f; //--seg-peak-height

	// opt->window_length1 = 3; //--seg-window-length1
    // opt->window_length2 = 6; //--seg-window-length2
    // opt->threshold1 = 1.4f; //--seg-threshold1
    // opt->threshold2 = 9.0f; //--seg-threshold2
    // opt->peak_height = 0.2f; //--seg-peak-height

	//TODO: RNA values:
	// opt->window_length1 = 7,
	// opt->window_length2 = 14,
	// opt->threshold1 = 2.5f,
	// opt->threshold2 = 9.0f,
	// opt->peak_height = 1.0f;

	//Sequence until parameters
	opt->t_threshold = 1.5f; //--threshold
	opt->tn_samples = 5; //--n-samples
	opt->ttest_freq = 500; //--test-frequency
	opt->tmin_reads = 500; //--min-reads
}

config_t parse_options(int argc, char *argv[]) {
	config_t config;
	const char *opt_str = "k:d:p:e:q:w:n:o:t:K:x:h";
	config.o = KETOPT_INIT;
	int c;
	char *s;
	FILE *fp_help = stderr;


	ri_verbose = 3;
	liftrlimit();
	ri_realtime0 = ri_realtime();
	ri_set_opt(0, &config.ipt, &config.opt);

	// test command line options and apply option -x/preset first
	while ((c = ketopt(&config.o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'x') {
			if (ri_set_opt(config.o.arg, &config.ipt, &config.opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", config.o.arg);
				return config;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return config;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[config.o.i - 1]);
			return config;
		}
	}
	config.o = KETOPT_INIT;

	while ((c = ketopt(&config.o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'd') config.fnw = config.o.arg;
		else if (c == 'p') config.fpore = config.o.arg;
		else if (c == 'k') config.ipt.k = atoi(config.o.arg);
		else if (c == 'e') config.ipt.e = atoi(config.o.arg);
		else if (c == 'q') config.ipt.q = atoi(config.o.arg);
		else if (c == 'w') config.ipt.w = atoi(config.o.arg);
		else if (c == 'n') config.ipt.n = atoi(config.o.arg);
		else if (c == 't') config.n_threads = atoi(config.o.arg);
		else if (c == 'v') ri_verbose = atoi(config.o.arg);
		else if (c == 'K') {config.opt.mini_batch_size = mm_parse_num(config.o.arg);}
		else if (c == 'h') fp_help = stdout;
		else if (c == 'o') {
			if (strcmp(config.o.arg, "-") != 0) {
				if (freopen(config.o.arg, "wb", stdout) == NULL) {
					fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", config.o.arg, strerror(errno));
					exit(1);
				}
			}
		}
		else if (c == 300) config.ipt.lev_col = atoi(config.o.arg);// --level_column
		else if (c == 301) { //--q-mid-occ
			config.opt.min_mid_occ = strtol(config.o.arg, &s, 10); // min
			if (*s == ',') config.opt.max_mid_occ = strtol(s + 1, &s, 10); //max
			// config.opt.q_mid_occ = atoi(config.o.arg);// --q-mid-occ
		}
		else if (c == 302) config.opt.mid_occ_frac = atof(config.o.arg);// --occ-frac
		else if (c == 303) config.opt.min_events = (uint32_t)atoi(config.o.arg); // --min-events
		else if (c == 304) config.opt.bw = atoi(config.o.arg);// --bw
		else if (c == 305) config.opt.max_target_gap_length = atoi(config.o.arg);// --max-target-gap
		else if (c == 306) config.opt.max_query_gap_length = atoi(config.o.arg);// --max-query-gap
		else if (c == 307) config.opt.min_num_anchors = atoi(config.o.arg);// --min-anchors
		else if (c == 308) config.opt.min_chaining_score = atoi(config.o.arg);// --min-score
		else if (c == 309) config.opt.chain_gap_scale = atof(config.o.arg);// --chain-gap-scale
		else if (c == 310) config.opt.chain_skip_scale = atof(config.o.arg);// --chain-skip-scale
		else if (c == 311) config.opt.best_n = atoi(config.o.arg);// --best-chains
		else if (c == 312) config.opt.mask_level = atof(config.o.arg);// --primary-ratio
		else if (c == 313) config.opt.mask_len = atoi(config.o.arg);// --primary-length
		else if (c == 314) config.opt.max_num_skips = atoi(config.o.arg);// --max-skips
		else if (c == 315) config.opt.max_chain_iter = atoi(config.o.arg);// --max-iterations
		else if (c == 316) config.opt.flag |= RI_M_RMQ; // --rmq
		else if (c == 317) config.opt.rmq_inner_dist = atoi(config.o.arg); // --rmq-inner-dist
		else if (c == 318) config.opt.rmq_size_cap = atoi(config.o.arg); // --rmq-size-cap
		else if (c == 319) config.opt.bw_long = atoi(config.o.arg);// --bw-long
		else if (c == 320) config.opt.max_num_chunk = atoi(config.o.arg);// --max-chunks
		else if (c == 321) config.opt.min_mapq = atoi(config.o.arg);// --min-mapq
		else if (c == 322) config.opt.alt_drop = atof(config.o.arg);// --alt-drop
		else if (c == 323) config.opt.w_besta = atof(config.o.arg);// --w-besta
		else if (c == 324) config.opt.w_bestma = atof(config.o.arg);// --w-bestma
		else if (c == 325) config.opt.w_bestq = atof(config.o.arg);// --w-bestq
		else if (c == 326) config.opt.w_bestmq = atof(config.o.arg);// --w-bestmq
		else if (c == 327) config.opt.w_bestmc = atof(config.o.arg);// --w-bestmc
		else if (c == 328) config.opt.w_threshold = atof(config.o.arg);// --w-threshold
		else if (c == 329) {
			config.opt.bp_per_sec = atoi(config.o.arg); config.opt.sample_per_base = (float)config.opt.sample_rate / config.opt.bp_per_sec;
			config.ipt.bp_per_sec = atoi(config.o.arg); config.ipt.sample_per_base = (float)config.ipt.sample_rate / config.ipt.bp_per_sec;
		}// --bp-per-sec
		else if (c == 330) {
			config.opt.sample_rate = atoi(config.o.arg); config.opt.sample_per_base = (float)config.opt.sample_rate / config.opt.bp_per_sec;
			config.ipt.sample_rate = atoi(config.o.arg); config.ipt.sample_per_base = (float)config.ipt.sample_rate / config.ipt.bp_per_sec;
		}// --sample-rate
		else if (c == 331) config.opt.chunk_size = atoi(config.o.arg);// --chunk-size
		else if (c == 332) {config.opt.window_length1 = atoi(config.o.arg); config.ipt.window_length1 = atoi(config.o.arg);}// --seg-window-length1
		else if (c == 333) {config.opt.window_length2 = atoi(config.o.arg); config.ipt.window_length2 = atoi(config.o.arg);}// --seg-window-length2
		else if (c == 334) {config.opt.threshold1 = atof(config.o.arg); config.ipt.threshold1 = atof(config.o.arg);}// --seg-threshold1
		else if (c == 335) {config.opt.threshold2 = atof(config.o.arg); config.ipt.threshold2 = atof(config.o.arg);}// --seg-threshold2
		else if (c == 336) {config.opt.peak_height = atof(config.o.arg); config.ipt.peak_height = atof(config.o.arg);}// --seg-peak-height
		else if (c == 337) config.opt.flag |= RI_M_SEQUENCEUNTIL;// --sequence-until
		else if (c == 338) config.opt.t_threshold = atof(config.o.arg);// --threshold
		else if (c == 339) config.opt.tn_samples = atoi(config.o.arg);// --n-samples
		else if (c == 340) config.opt.ttest_freq = atoi(config.o.arg);// --test-frequency
		else if (c == 341) config.opt.tmin_reads = atoi(config.o.arg);// --min-reads
		else if (c == 342) config.opt.mid_occ_frac = atof(config.o.arg);// --occ-frac
		else if (c == 343) { // --depletion
			config.opt.best_n = 5; config.opt.min_mapq = 10; config.opt.w_threshold = 0.50f;
			config.opt.min_num_anchors = 2; config.opt.min_chaining_score = 15; config.opt.chain_skip_scale = 0.0f;
		}
		else if (c == 344) {config.ipt.flag |= RI_I_STORE_SIG;} // --store-sig
		else if (c == 345) {config.ipt.flag |= RI_I_SIG_TARGET;} // --sig-target
		else if (c == 346) {config.opt.flag |= RI_M_NO_ADAPTIVE;} // --disable-adaptive
		else if (c == 347) {config.ipt.diff = atof(config.o.arg);} // --sig-diff
		else if (c == 348) {config.opt.flag |= RI_M_ALIGN;} // --align
		else if (c == 349) config.opt.flag |= RI_M_DTW_EVALUATE_CHAINS; // --dtw-evaluate-chains
		else if (c == 350) config.opt.flag |= RI_M_DTW_OUTPUT_CIGAR; // --dtw-output-cigar
		else if (c == 351) { //--dtw-border-constraint
			if(ri_mapopt_parse_dtw_border_constraint(&config.opt, config.o.arg) != 0){
				fprintf(stderr, "[ERROR] unknown DTW border constraint in \"%s\"\n", argv[config.o.i - 1]);
				return config;
			}
		}
		else if (c == 352) config.opt.flag |= RI_M_DTW_LOG_SCORES; // --dtw-log-scores
		else if (c == 353) config.opt.flag |= RI_M_DISABLE_CHAININGSCORE_FILTERING; // --no-chainingscore-filtering
		else if (c == 354) config.opt.dtw_match_bonus = atof(config.o.arg); // --dtw-match-bonus
		else if (c == 355) config.opt.flag |= RI_M_OUTPUT_CHAINS; // --output-chains
		else if (c == 356) { //dtw-fill-method
			if(ri_mapopt_parse_dtw_fill_method(&config.opt, config.o.arg) != 0){
				fprintf(stderr, "[ERROR] unknown DTW fill method in \"%s\"\n", argv[config.o.i - 1]);
				return config;
			}
		}
		else if (c == 357) config.opt.dtw_min_score = atof(config.o.arg); // --dtw-min-score
		else if (c == 358) config.opt.flag |= RI_M_LOG_ANCHORS; // --log-anchors
		else if (c == 359) config.opt.flag |= RI_M_LOG_NUM_ANCHORS; // --log-num-anchors
		else if (c == 360) config.opt.rev_col_limit = atoi(config.o.arg); // --rev-collision-count
		else if (c == 361) config.opt.chn_rev_bump = atof(config.o.arg); // --chn-rev-bump
		// else if (c == 362) {config.ipt.flag |= RI_I_REV_QUERY;}// --rev-query
		else if (c == 363) { // --r10
			config.ipt.k = 9;

			config.ipt.window_length1 = 3; config.ipt.window_length2 = 6;
			config.ipt.threshold1 = 6.5f; config.ipt.threshold2 = 4.0f;
			config.ipt.peak_height = 0.2f;

			config.opt.window_length1 = 3; config.opt.window_length2 = 6;
			config.opt.threshold1 = 6.5f; config.opt.threshold2 = 4.0f;
			config.opt.peak_height = 0.2f;

			config.opt.chain_gap_scale = 1.2f;

			config.opt.bp_per_sec = 400;
			config.ipt.bp_per_sec = 400;
			config.opt.sample_rate = 5000; config.opt.sample_per_base = (float)config.opt.sample_rate / config.opt.bp_per_sec;
			config.ipt.sample_rate = 5000; config.ipt.sample_per_base = (float)config.ipt.sample_rate / config.ipt.bp_per_sec;

			// io->fine_range = 0.6;
			// mo->min_mapq = 5, mo->min_chaining_score = 10, mo->chain_gap_scale = 0.6f;
		}
		else if (c == 364) {config.ipt.fine_min = atof(config.o.arg);}// --fine-min
		else if (c == 365) {config.ipt.fine_max = atof(config.o.arg);}// --fine-max
		else if (c == 366) {config.ipt.fine_range = atof(config.o.arg);}// --fine-range
		else if (c == 367) {config.ipt.flag |= RI_I_OUT_QUANTIZE; config.ipt.flag |= RI_I_SIG_TARGET;}// --out-quantize
		else if (c == 368) {config.ipt.flag |= RI_I_NO_EVENT_DETECTION;}// --no-event-detection
		else if (c == 369) {config.io_n_threads = atoi(config.o.arg);}// --io-thread
		else if (c == 370) config.opt.min_chaining_score2 = atoi(config.o.arg);// --min-score2
		else if (c == 371) {puts(RH_VERSION); return config;}// --version
		else if (c == 'V') {puts(RH_VERSION); return config;}
	}

	if (argc == config.o.ind || fp_help == stdout) {
		fprintf(fp_help, "Usage: rawhash [options] <target.fa>|<target.idx> [query.fast5] [...]\n");
		fprintf(fp_help, "Options:\n");

		fprintf(fp_help, "    --version     show version number\n");

		fprintf(fp_help, "  K-mer (pore) Model:\n");
		fprintf(fp_help, "    -p FILE      pore model FILE [].\n");
		fprintf(fp_help, "    -k INT       size of the k-mers in the pore model [%d]. This is usually 6 for R9.4 and 9 for R10.\n", config.ipt.k);
		fprintf(fp_help, "    --level_column INT       0-based column index where the mean values are stored in the pore file [%d]. This is usually 1 for both R9.4 and R10.\n", config.ipt.lev_col);

		fprintf(fp_help, "\n  Indexing:\n");
		fprintf(fp_help, "    -d FILE     [Strongly recommended to create before mapping] dump index to FILE [].\n");
		fprintf(fp_help, "    -e INT     number of events concatanated in a single hash (usually no larger than 10). Also applies during mapping [%d].\n", config.ipt.e);
		fprintf(fp_help, "    -q INT     Number of bits to use for quantization [%d]. Number of quantized buckets are created accordingly (2^INT).\n", config.ipt.q);
		fprintf(fp_help, "    -w INT     minimizer window size [%d]. Enables minimizer-based seeding in indexing and mapping (may reduce accuracy but improves the performance and memory space efficiency).\n", config.ipt.w);
		fprintf(fp_help, "    --store-sig      Stores the target signal in the index file.\n");
		fprintf(fp_help, "    --sig-target     The target sequence (reference) contains signals rather than base characters.\n");
		fprintf(fp_help, "    --sig-diff FLOAT    [Advanced] Signal value (FLOAT) difference between two consecutive events to be packed together in a single hash value [%g].\n", config.ipt.diff);

		// fprintf(fp_help, "    -n NUM     number of consecutive seeds to use for BLEND-based seeding [%d]. Enables the BLEND mechanism (may improve accuracy but reduces the performance at the moment)\n", config.ipt.n);

		fprintf(fp_help, "\n  Seeding:\n");
		fprintf(fp_help, "    --q-mid-occ INT1[,INT2]     Lower and upper bounds of k-mer occurrences [%d, %d]. The final k-mer occurrence threshold is max{INT1, min{INT2, --occ-frac}}. This option prevents excessively small or large -f estimated from the input reference.\n", config.opt.min_mid_occ, config.opt.max_mid_occ);
		// fprintf(fp_help, "    --occ-frac FLOAT     Discard a query seed if its occurrence is higher than FLOAT fraction of all query seeds [%g]. Set 0 to disable. [Note: Both --q-mid-occ and --occ-frac should be met for a seed to be discarded].\n", config.opt.q_occ_frac);

		fprintf(fp_help, "\n  Chaining Parameters:\n");
		fprintf(fp_help, "    --min-events INT     minimum number of INT events in a chunk to start chain elongation [%u].\n", config.opt.min_events);
		fprintf(fp_help, "    --bw INT     maximum INT gap length in a chain [%d].\n", config.opt.bw);
		fprintf(fp_help, "    --max-target-gap INT     maximum INT target gap length in a chain [%d].\n", config.opt.max_target_gap_length);
		fprintf(fp_help, "    --max-query-gap INT     maximum INT query gap length in a chain [%d].\n", config.opt.max_query_gap_length);
		fprintf(fp_help, "    --min-anchors INT     chain is discarded if it contains less than INT number of anchors [%d].\n", config.opt.min_num_anchors);
		fprintf(fp_help, "    --best-chains INT     best INT secondary chains to keep with their primary chains when making the mapping decisions [%d]\n", config.opt.best_n);
		fprintf(fp_help, "    --min-score INT     chain is discarded if its score is < INT [%d]\n", config.opt.min_chaining_score);
		fprintf(fp_help, "    --chain-gap-scale FLOAT     [Advanced] Determines [chain gap penalty] = FLOAT * 0.01 * e  [%g]\n", config.opt.chain_gap_scale);
		fprintf(fp_help, "    --chain-skip-scale FLOAT     [Advanced] Determines [chain skip penalty] = FLOAT * 0.01 * e  [%g]\n", config.opt.chain_skip_scale);
		// fprintf(fp_help, "    --chain-match-score INT     [Advanced] Match score (Used in MAPQ and Primary chain identification) [%d]\n", config.opt.a);
		fprintf(fp_help, "    --primary-ratio FLOAT     [Advanced] The chain is primary if its region ratio uncovered by other chains is larger than FLOAT [%g]\n", config.opt.mask_level);
		fprintf(fp_help, "    --primary-length INT     [Advanced] The chain is primary if its region length uncovered by other chains is larger than INT [%d]\n", config.opt.mask_len);
		fprintf(fp_help, "    --max-skips INT     [Advanced] stop looking for a predecessor for an anchor if the best predecessor is not updated after INT many iterations [%d]\n", config.opt.max_num_skips);
		fprintf(fp_help, "    --max-iterations INT     [Advanced] maximum INT number predecessor anchors to check to calculate the best score for an anchor [%d]\n", config.opt.max_chain_iter);
		fprintf(fp_help, "    --rmq     [Advanced] Uses RMQ-based chaining. Faster but less accurate than default (DP)\n");
		fprintf(fp_help, "    --rmq-inner-dist INT     [Advanced] RMQ inner distance [%d]\n", config.opt.rmq_inner_dist);
		fprintf(fp_help, "    --rmq-size-cap INT     [Advanced] RMQ cap size [%d]\n", config.opt.rmq_size_cap);
		fprintf(fp_help, "    --bw-long INT     [Advanced] maximum long INT gap length to re-chain the chains. Disabled by default. To enable, set it to larger than --bw [%d]\n", config.opt.bw_long);

		fprintf(fp_help, "\n  DTW Parameters (as introduced in RawAlign):\n");
		fprintf(fp_help, "    --dtw-evaluate-chains     evaluate chains using DTW. Note, the index must be built using --store-sig for this functionality to work [%s]\n", config.opt.flag & RI_M_DTW_EVALUATE_CHAINS? "yes" : "no");
		fprintf(fp_help, "    --dtw-output-cigar     output CIGAR string for DTW [%s]\n", config.opt.flag & RI_M_DTW_OUTPUT_CIGAR? "yes" : "no");
		fprintf(fp_help, "    --dtw-border-constraint STR     DTW border constraint: 'global', 'sparse' (i.e., align only between anchors), 'local' [%s]\n", ri_maptopt_dtw_mode_to_string(config.opt.dtw_border_constraint));
		fprintf(fp_help, "    --dtw-match-bonus FLOAT     DTW match bonus FLOAT [%g]\n", config.opt.dtw_match_bonus);

		fprintf(fp_help, "\n  Mapping Decisions (Mapping and sequencing is stopped after taking any of these decisions):\n");
		fprintf(fp_help, "    --max-chunks INT     stop mapping (read not mapped) after sequencing INT number of chunks [%u]\n", config.opt.max_num_chunk);
		fprintf(fp_help, "    --min-mapq INT     map the read if there is only one chain and its MAPQ > INT [%d]\n", config.opt.min_mapq);
		fprintf(fp_help, "    --disable-adaptive     Disables stopping the read early and rather lets the read to be sequenced fully to make the analysis. This is not activated by default.\n");

		fprintf(fp_help, "\n  Nanopore Parameters:\n");
		fprintf(fp_help, "    --bp-per-sec INT     DNA molecules transiting through the pore (bp per second) [%u]\n", config.opt.bp_per_sec);
		fprintf(fp_help, "    --sample-rate INT     current sample rate in Hz [%u]\n", config.opt.sample_rate);
		fprintf(fp_help, "    --chunk-size INT     current samples in a single chunk (by default set to the amount of signals sampled in 1 second) [%u]\n", config.opt.chunk_size);

		fprintf(fp_help, "    --seg-window-length1 INT     [Advanced] First window length in segmentation [%u]\n", config.opt.window_length1);
		fprintf(fp_help, "    --seg-window-length2 INT     [Advanced] Second window length in segmentation [%u]\n", config.opt.window_length2);
		fprintf(fp_help, "    --seg-threshold1 FLOAT     [Advanced] Peak value threshold for the first window in segmentation [%g]\n", config.opt.threshold1);
		fprintf(fp_help, "    --seg-threshold2 FLOAT     [Advanced] Peak value threshold for the first window in segmentation [%g]\n", config.opt.threshold2);
		fprintf(fp_help, "    --seg-peak-height FLOAT     [Advanced] Peak height than the current signal to confirm the peak point in segmentation [%g]\n", config.opt.peak_height);

		fprintf(fp_help, "\n  Sequence Until Parameters:\n");
		fprintf(fp_help, "    --sequence-until     Activates Sequence Until and performs real-time relative abundance calculations. The computation will stop as soon as an estimation with high confidence is reached without processing further reads from the set.\n");
		fprintf(fp_help, "    --threshold FLOAT     outliers are determined if cross-correlation distance > FLOAT [%g]. Sequencing will stop if there are no outliers in the sample of estimations.\n", config.opt.t_threshold);
		fprintf(fp_help, "    --n-samples INT     New estimation is tested against INT many previous estimations [%u]\n", config.opt.tn_samples);
		fprintf(fp_help, "    --test-frequency INT     Make a new estimation after every INT reads [%u]\n", config.opt.ttest_freq);
		fprintf(fp_help, "    --min-reads INT     Minimum number of reads to sequence before making the first estimation [%u]\n", config.opt.tmin_reads);

		fprintf(fp_help, "\n  Input/Output:\n");
		fprintf(fp_help, "    -o FILE     output mappings to FILE [stdout]\n");
		fprintf(fp_help, "    -t INT      number of threads [%d]\n", config.n_threads);
		fprintf(fp_help, "    --io-thread INT      number of threads allocated for IO operations (i.e., reading from a file) out of all threads that will be used for this run (-t). Only available for S/BLOW5 files for now. INT must be smaller than the overall number of threads (-t) [%d]\n", config.io_n_threads);
		fprintf(fp_help, "    -K NUM      minibatch size for mapping [500M]. Increasing this value may increase thread utilization. If there are many larger FAST5 files, it is recommended to keep this value between 500M - 5G to use less memory while utilizing threads nicely.\n");
//		fprintf(fp_help, "    -v INT     verbose level [%d]\n", ri_verbose);

		fprintf(fp_help, "\n  Experimental/Debugging Parameters:\n");
		fprintf(fp_help, "    --out-quantize     	Output the quantized values from raw signals provided as input. Mapping is not performed and the index file is not needed.\n");
		fprintf(fp_help, "    --no-event-detection  Do not perform event detection. This can be set if your raw signal is already segmented.\n");

		fprintf(fp_help, "\n  Presets:\n");
		fprintf(fp_help, "    --depletion     Should be used for quickly depleting organisms for use cases that require high precision (e.g., for contamination analysis or relative abundance estimation). Can be used with or without the -x preset (--best-chains 5 --min-mapq 10 --w-threshold 0.5 --min-anchors 2 --min-score 15 --chain-skip-scale 0).\n");
		fprintf(fp_help, "    --r10     Sets the segmentation parameters for R10.4.1. Can be used with or without the -x preset (-k9 --seg-window-length1 3 --seg-window-length2 6 --seg-threshold1 6.5 --seg-threshold2 4 --seg-peak-height 0.2 --chain-gap-scale 1.2).\n");
		fprintf(fp_help, "    -x STR     preset (always applied before other options) []\n");
		fprintf(fp_help, "                 - viral     Enables accurate mapping to very small genomes such as viral genomes (-e 6 -q 4 --max-chunks 5 --bw 100 --max-target-gap 500 --max-target-gap 500 --min-score 10 --chain-gap-scale 1.2 --chain-skip-scale 0.3).\n");
		fprintf(fp_help, "                 - sensitive     Enables sensitive mapping. Suitable when working with small genomes of size < 500M (-e 8 -q 4 --fine-range 0.4).\n");
		fprintf(fp_help, "                 - fast     Enables fast mapping with slightly reduced accuracy. Suitable when reads are mapped to large genomes of size > 500M and < 5Gb (-e 8 -q 4 --max-chunks 20).\n");
		fprintf(fp_help, "                 - faster     Enables faster mapping than '-x fast' and reduced memory space usage for indexing with slightly reduced accuracy. This mechanism uses the minimizer sketching technique and should be used when '-x fast' cannot meet the real-time requirements for a particular genome (e.g., for very large genomes > 5Gb)\n");
		fprintf(fp_help, "\n  Rawsamble Presets:\n");
		fprintf(fp_help, "                 - ava     	 All-vs-all overlapping mode (default for Rawsamble).\n");
		fprintf(fp_help, "                 - ava-sensitive     	 More sensitive All-vs-all overlapping mode. Can be slightly slower than -ava but likely to generate longer unitigs in downstream asssembly.\n");
		fprintf(fp_help, "                 - ava-viral     	 All-vs-all overlapping for very small genomes such as viral genomes.\n");
		fprintf(fp_help, "                 - ava-large     	 All-vs-all overlapping for large genomes of size > 10Gb\n");

		// fprintf(fp_help, "\nSee `man ./rawhash.1' for detailed description of these and other advanced command-line options.\n");
		return config;
	}

	if(config.n_threads < config.io_n_threads){
		fprintf(stderr, "[ERROR] The overall number of threads (-t [%d]) must NOT be smaller than the number of IO threads (--io-thread [%d).\n", config.n_threads, config.io_n_threads);
		return config;
	}

	if(config.ipt.w && config.ipt.n){
		fprintf(stderr, "[ERROR] minimizer window 'w' ('%d') and BLEND 'neighbor' ('%d') values cannot be set together. At least one of them must be zero to enable one of the seeding options: %s\n", config.ipt.w, config.ipt.n, strerror(errno));
		return config;
	}

	config.valid = true;
	return config;
}