//
// Created by sayan on 3/2/26.
//

#include <string>
#include <vector>
#include <stdexcept>
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/io.h"

#include "rawhash.h"

using namespace std;
namespace py = pybind11;

static inline char* to_c_str(std::string &s) {
    auto cs = new char[s.size()+1];
    strcpy(cs, s.c_str());
    return cs;
}

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

struct Request {
    int channel = 0;
    string id;
    py::array_t<float> signal;
    Request() = default;
    Request(int channel, string &id, py::array_t<float> &signal): channel(channel), id(id), signal(signal) {}
};

struct RawRequest {
    int channel;
    string id;
    vector<float> signal;   // plain C++ — no Python reference kept
};

struct Response {
    int channel = 0;
    string id;
    Alignment alignment;
    Response() = default;
    Response(int channel, string &id, Alignment &alignment): channel(channel), id(id), alignment(alignment) {}
};

struct ResponseGenerator {
    parlay::sequence<Response> responses;
    explicit ResponseGenerator(parlay::sequence<Response> &&responses): responses(std::move(responses)) {}
    Response next() {
        if (responses.empty())
            throw py::stop_iteration();
        else {
            auto response = responses.back();
            responses.pop_back();
            return response;
        }
    }
    ResponseGenerator& iter() {
        return *this;
    }
};

struct Index {
    config_t config;
    ri_idx_reader_t *idx_rdr = nullptr;
    ri_idx_t *ri = nullptr;
    ri_pore_t pore{};

    explicit Index(vector<char*> &argv) {
        load(argv.size(), argv.data());
    }

    explicit Index(const py::list &args) {
        vector<char*> argv;
        for (auto& arg: args) {
            auto s = arg.cast<string>();
            argv.push_back(to_c_str(s));
        }
        load(argv.size(), argv.data());
    }

    void load(int argc, char *argv[]) {
        config = parse_options(argc, argv);
        if (!config.valid) throw std::invalid_argument("Config is invalid.");
        idx_rdr = ri_idx_reader_open(argv[config.o.ind], &config.ipt, config.fnw);
        if (idx_rdr == 0) {
            log_error("Failed to open file '%s': %s", argv[config.o.ind], strerror(errno));
            exit(1);
        }
        pore.pore_vals = nullptr;
        pore.pore_inds = nullptr;
        pore.max_val = -5000.0;
        pore.min_val = 5000.0;
        ri = ri_idx_reader_read(idx_rdr, &pore, config.n_threads, config.io_n_threads);
        if (ri == 0)
            log_error("Could not load index from file '%s': %s", argv[config.o.ind], strerror(errno));
        ri_mapopt_update(&config.opt, ri);
    }

    void validate() const {
        ri_idx_stat(ri);
    }

    Alignment query(py::array_t<float> &signal) {
        ri_reg1_t *reg = map_signal(signal.data(), signal.size(), ri, &config.opt);
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

    Alignment query_raw(vector<float> &signal) {
        ri_reg1_t *reg = map_signal(signal.data(), signal.size(), ri, &config.opt);
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

    ResponseGenerator query_stream(const py::iterator& reads) {
	    std::vector<RawRequest> requests;
        {
            py::gil_scoped_acquire gil;
            for (auto &read: reads) {
                auto req = read.cast<Request>();
                auto buf = req.signal.request();
                auto *ptr = static_cast<float*>(buf.ptr);
                requests.push_back({
                    req.channel,
                    req.id,
                    vector<float>(ptr, ptr + buf.size)
                });
            }
        }
        auto responses = parlay::tabulate(requests.size(), [&](size_t i) {
            auto alignment = query_raw(requests[i].signal);
            return Response(requests[i].channel, requests[i].id, alignment);
        });
        return ResponseGenerator(std::move(responses));
    }

};

PYBIND11_MODULE(_core, m) {
    py::class_<Alignment>(m, "Alignment")
            .def(py::init<>())  // Default constructor
            .def(py::init<const char*, bool, int, float, int>(),  // Parameterized constructor
                 py::arg("header"), py::arg("fwd"), py::arg("start"),
                 py::arg("pres_frac"), py::arg("qry_len"))
            .def_readonly("ctg", &Alignment::ctg)
            .def_readonly("r_st", &Alignment::r_st)
            .def_readonly("r_en", &Alignment::r_en)
            .def_readonly("strand", &Alignment::strand)
            .def_readonly("pres_frac", &Alignment::pres_frac);

    py::class_<Index>(m, "Index")
            .def(py::init<const py::list&>())
            .def("validate", &Index::validate)
            .def("query", &Index::query)
            .def("query_stream", &Index::query_stream);

    py::class_<Request>(m, "Request")
            .def(py::init<int, string&, py::array_t<float> &>(), py::arg("channel"), py::arg("id"), py::arg("signal"))
            .def_readwrite("channel", &Request::channel)
            .def_readwrite("id", &Request::id)
            .def_readwrite("signal", &Request::signal);

    py::class_<Response>(m, "Response")
            .def(py::init<int, string&, Alignment&>(), py::arg("channel"), py::arg("id"), py::arg("alignment"))
            .def_readwrite("channel", &Response::channel)
            .def_readwrite("id", &Response::id)
            .def_readwrite("alignment", &Response::alignment);

    py::class_<ResponseGenerator>(m, "ResponseGenerator")
            .def("__iter__", &ResponseGenerator::iter)
            .def("__next__", &ResponseGenerator::next);
}
