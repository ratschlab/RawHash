#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "rawhash.h"
#include "ketopt.h"

int main(int argc, char *argv[]) {
	ri_idx_reader_t *idx_rdr;
	ri_idx_t *ri;
	config_t config = parse_options(argc, argv);
	if (!config.valid) return 1;

	idx_rdr = ri_idx_reader_open(argv[config.o.ind], &config.ipt, config.fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[config.o.ind], strerror(errno));
		return 1;
	}

	if (!idx_rdr->is_idx && config.fnw == 0 && argc - config.o.ind < 2 && !(config.ipt.flag&RI_I_OUT_QUANTIZE)) {
		fprintf(stderr, "[ERROR] missing input: please specify a query FAST5/SLOW5/POD5 file(s) to map or option -d to store the index in a file before running the mapping\n");
		ri_idx_reader_close(idx_rdr);
		return 1;
	}


	ri_pore_t pore;
	pore.pore_vals = NULL;
	pore.pore_inds = NULL;
	pore.max_val = -5000.0;
	pore.min_val = 5000.0;
	if(!(config.ipt.flag&RI_I_OUT_QUANTIZE)){
		if((!idx_rdr->is_idx && config.fpore == 0) && !(!(config.ipt.flag&RI_I_REV_QUERY) && config.ipt.flag&RI_I_SIG_TARGET)){
			fprintf(stderr, "[ERROR] missing input: please specify a pore model file with -p when generating the index from a sequence file\n");
			ri_idx_reader_close(idx_rdr);
			return 1;
		}else if(!idx_rdr->is_idx && config.fpore){
			load_pore(config.fpore, config.ipt.k, config.ipt.lev_col, &pore);
			if(!pore.pore_vals){
				fprintf(stderr, "[ERROR] cannot parse the k-mer pore model file. Please see the example k-mer model files provided in the RawHash repository.\n");
				ri_idx_reader_close(idx_rdr);
				return 1;
			}
		}
	}

	while ((ri = ri_idx_reader_read(idx_rdr, &pore, config.n_threads, config.io_n_threads)) != 0) {
		int ret;
		if (ri_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, ri_realtime() - ri_realtime0, ri_cputime() / (ri_realtime() - ri_realtime0), ri->n_seq);
		if (argc != config.o.ind + 1) ri_mapopt_update(&config.opt, ri);
		if (ri_verbose >= 3) ri_idx_stat(ri);
		if (argc - (config.o.ind + 1) == 0) {
			fprintf(stderr, "[INFO] No files to query index on. Only the index is constructed.\n");
			ri_idx_destroy(ri);
			continue; // no query files, just creating the index
		}
		ret = 0;
		// if (!(config.opt.flag & MM_F_FRAG_MODE)) { //TODO: enable frag mode directly from options
		// for (i = o.ind + 1; i < argc; ++i) {
		// 	ret = ri_map_file(ri, argv[i], &config.opt, n_threads, io_n_threads);
		// 	if (ret < 0) break;
		// }
		// }
		// else { //TODO: enable frag mode directly from options
			ret = ri_map_file_frag(ri, argc - (config.o.ind + 1), (const char**)&argv[config.o.ind + 1], &config.opt, config.n_threads, config.io_n_threads);
		// }
		ri_idx_destroy(ri);
		if (ret < 0) {
			fprintf(stderr, "ERROR: failed to map the query file\n");
			exit(EXIT_FAILURE);
		}
	}
	// n_parts = idx_rdr->n_parts;
	ri_idx_reader_close(idx_rdr);
	if(pore.pore_vals)free(pore.pore_vals);
	if(pore.pore_inds)free(pore.pore_inds);

	if (fflush(stdout) == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}

	if (ri_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, RH_VERSION);
		// fprintf(stderr, "[M::%s] CMD:", __func__);
		// for (i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, ri_realtime() - ri_realtime0, ri_cputime(), ri_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
