//
// Created by sayan on 3/2/26.
//

#include <string>
#include <vector>
#include <stdexcept>
#include "slow5/slow5.h"

#include "rawhash.h"

using namespace std;

struct Alignment {
    string ctg;
    int r_st = 0, r_en = 0, strand = 1;
    float pres_frac = 0.0f;

    Alignment() = default;

    Alignment(const char *header, bool fwd, int start, float pres_frac, int qry_len):
        ctg(header), r_st(start), r_en(start + qry_len), strand(fwd?1:-1), pres_frac(pres_frac) {}

    bool valid() const { return this->ctg != "*"; }
};

const Alignment NO_ALIGNMENT = {"*", true, 0, 0.0f, 0};

struct Index {
    config_t config;
    ri_idx_reader_t *idx_rdr;
    ri_idx_t *ri;
    ri_pore_t pore{};

    explicit Index(vector<char*> &argv) {
        config = parse_options(argv.size(), argv.data());
        if (!config.valid) throw std::invalid_argument("Config is invalid.");
        idx_rdr = ri_idx_reader_open(argv[config.o.ind], &config.ipt, config.fnw);
        if (idx_rdr == 0) {
            fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[config.o.ind], strerror(errno));
            exit(1);
        }
        pore.pore_vals = nullptr;
        pore.pore_inds = nullptr;
        pore.max_val = -5000.0;
        pore.min_val = 5000.0;
        ri = ri_idx_reader_read(idx_rdr, &pore, config.n_threads, config.io_n_threads);
        if (ri == 0) {
            fprintf(stderr, "[ERROR] Could not load index from file '%s': %s\n", argv[config.o.ind], strerror(errno));
            exit(1);
        }
        ri_mapopt_update(&config.opt, ri);
        ri_idx_stat(ri);
        fprintf(stderr, "Index is ready..\n");
    }

    Alignment query(std::vector<float> &signal) {
        ri_reg1_t *reg = map_signal(signal, ri, &config.opt);
        if (!reg) return NO_ALIGNMENT;

        Alignment result = NO_ALIGNMENT;
        if (reg->n_maps > 0 && reg->maps[0].mapped) {
            const ri_map_t &m = reg->maps[0];
            const char *header = (ri->flag & RI_I_SIG_TARGET)
                ? ri->sig[m.ref_id].name
                : ri->seq[m.ref_id].name;
            result = Alignment(header, !m.rev, (int)m.fragment_start_position,
                               (float)m.mapq / 60.0f, (int)m.read_end_position);
        }

        if (reg->maps) { free(reg->maps); reg->maps = nullptr; }
        free(reg);
        return result;
    }
};

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))

int main(int argc, char **argv) {
    // create a vector from the argv
    std::vector<char*> args(argv, argv + argc);
    auto index = Index(args);
    sleep(1);
    const char* FILE_PATH = "/data/NASExperiments/data/simulated/Zymo/signals/Sigs1_180.blow5";

    slow5_file_t *sp = slow5_open(FILE_PATH,"r");
    if(sp==NULL){
        fprintf(stderr,"Error in opening file\n");
        exit(EXIT_FAILURE);
    }

    slow5_rec_t *rec = NULL; //slow5 record to be read
    int ret=0; //for return value

    int n_aligned = 0, n_read = 0;

    //iterate through the file until end
    vector<float> signal;
    printf("%10s %10s\n", "Aligned", "Unaligned");
    while((ret = slow5_get_next(&rec,sp)) >= 0){
        // printf("%s\t",rec->read_id);
        uint64_t len_raw_signal = rec->len_raw_signal;
        signal.resize(len_raw_signal);
        for(uint64_t i=0;i<len_raw_signal;i++){ //iterate through the raw signal and print in picoamperes
            signal[i] = TO_PICOAMPS(rec->raw_signal[i],rec->digitisation,rec->offset,rec->range);
            // printf("%f ",pA);
        }
        // printf("\n");
        n_read++;
        auto alignment = index.query(signal);
        if (alignment.valid()) n_aligned++;
        if (n_read % 100 == 0) {
            printf("\r%10d %10d", n_aligned, n_read - n_aligned);
            fflush(stdout);
        }
    }
    printf("\n");

    if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
        fprintf(stderr,"Error in slow5_get_next. Error code %d\n",ret);
        exit(EXIT_FAILURE);
    }

    //free the SLOW5 record
    slow5_rec_free(rec);

    //close the SLOW5 file
    slow5_close(sp);
}