#include <pybind11/pybind11.h>
#include "litchi_eat.hpp"

namespace py = pybind11;


//Nside, rankA, rankB, curvIndex, mint, maxt, numt, smooth, linThresh, function, forceOutname, sequence, maskname, maskThresh
PYBIND11_MODULE(litchieat, m) {
    py::class_<paramStruct>(m, "paramStruct")
        .def(py::init<>())
        .def_readwrite("Nside", &paramStruct::Nside)
        .def_readwrite("rankA", &paramStruct::rankA)
        .def_readwrite("rankB", &paramStruct::rankB)
        .def_readwrite("curvIndex", &paramStruct::curvIndex)
        .def_readwrite("mint", &paramStruct::mint)
        .def_readwrite("maxt", &paramStruct::maxt)
        .def_readwrite("numt", &paramStruct::numt)
        .def_readwrite("smooth", &paramStruct::smooth)
        .def_readwrite("linThresh", &paramStruct::linThresh)
        .def_readwrite("function", &paramStruct::function)
        .def_readwrite("forceOutname", &paramStruct::forceOutname)
        .def_readwrite("sequence", &paramStruct::sequence)
        .def_readwrite("maskname", &paramStruct::maskname)
        .def_readwrite("maskThresh", &paramStruct::maskThresh)
        .def("__copy__",[](const paramStruct& self){
            return paramStruct(self);
        })
        .def("__deepcopy__",[](const paramStruct& self, py::dict){
            return paramStruct(self);
        });
    m.def("makeSingleMinkmap", &makeSingleMinkmap, "Reads file given by inname, creates Minkmap according to params, writes to outname", py::arg("inname"), py::arg("params"), py::arg("outname"));
    m.def("makeSequence", &makeSequence, "Reads file given by inname, creates threshold sequence of Minkmaps according to params, writes to outname", py::arg("inname"), py::arg("params"), py::arg("outname"));
}
