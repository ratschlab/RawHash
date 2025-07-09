//
// Created by Sayan Goswami on 27.06.2025.
//

#include "rawhash_wrapper.hpp"
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "kthread.h"
#include "numeric"
#include "execution"

namespace py = pybind11;

static inline char* to_c_str(std::string &s) {
    auto cs = new char[s.size()+1];
    strcpy(cs, s.c_str());
    return cs;
}

static std::vector<char*> kwargs_to_argv(const py::args& args, const py::kwargs& kwargs) {
    std::vector<std::string> arg_strings;

    // First argument is conventionally the program name
    arg_strings.emplace_back("program_name");

    // Handle positional arguments (flags)
    for (auto& arg : args) {
        std::string flag = "--" + arg.cast<std::string>();  // Convert to "--flag" format
        arg_strings.emplace_back(flag);
    }

    // Handle keyword arguments ("--key=value")
    for (auto& item : kwargs) {
        auto key = item.first.cast<std::string>();
        std::string value = py::str(item.second);  // Convert value to string
        if (key.size()==1) arg_strings.emplace_back("-" + key + "=" + value);
        else arg_strings.emplace_back("--" + key + "=" + value);
    }

    std::vector<char*> argv;  // Store pointers to C-style strings
    for (auto &s: arg_strings)
        argv.push_back(to_c_str(s));

    // Now `argc` and `argv.data()` can be used in a function expecting C-style args
    return argv;
}

struct Index {
    explicit Index(const py::args& args, const py::kwargs& kwargs) {
        auto argvec = kwargs_to_argv(args, kwargs);
        mapper = new RawHashMapper(argvec.size(), argvec.data());
    }

    std::vector<Alignment> query(std::vector<float> signal) {
        return mapper->map(signal.data(), signal.size());
    }

    std::vector<std::vector<Alignment>> query_batch(const py::list& signals) {
        auto nr = signals.size();
        std::vector<std::vector<Alignment>> results(nr);
        std::vector<int> id(nr);
        std::iota(id.begin(), id.end(), 0);
        std::for_each(std::execution::par, id.begin(), id.end(), [&](int i){
            auto signal = signals[i].cast<std::vector<float>>();
            results[i] = mapper->map(signal.data(), signal.size());
        });
        return results;
    }

private:
    RawHashMapper *mapper;
};

PYBIND11_MODULE(_core, m) {
    py::class_<Alignment>(m, "Alignment")
            .def(py::init<>())  // Default constructor
            .def_readonly("ctg", &Alignment::ctg)
            .def_readonly("r_st", &Alignment::r_st)
            .def_readonly("r_en", &Alignment::r_en)
            .def_readonly("strand", &Alignment::strand);

    py::class_<Index>(m, "Index")
            .def(py::init<const py::args&, const py::kwargs&>())
            .def("query", &Index::query)
            .def("query_batch", &Index::query_batch);
}