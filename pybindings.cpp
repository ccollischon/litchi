//c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) pybindings.cpp -o ../../example(python3-config --extension-suffix)
//g++-10 -O3 -Wall -shared -std=c++20 -fPIC $(python3 -m pybind11 --includes) pybindings.cpp -o ../../litchieat(python3-config --extension-suffix) -fopenmp -lhealpix_cxx
//#include <iostream>
#include <pybind11/pybind11.h>
#include "litchi_eat.hpp"

namespace py = pybind11;

//int add(int i=1, int j=2) {
//    return i + j;
//}
//
//struct Pet {
//    Pet(const std::string &name) : name(name) { }
//    void setName(const std::string &name_) { name = name_; }
//    const std::string &getName() const { return name; }
//
//    std::string name;
//};
//
//void petname(std::string name = "asdf") {
//    Pet p(name);
//    std::cout << p.getName() << std::endl;
//}
//
//PYBIND11_MODULE(example, m) {
//    m.doc() = "pybind11 example plugin"; // optional module docstring
//
//    //py::class_<Pet>(m, "Pet")
//    //    .def(py::init<const std::string &>())
//    //    .def("setName", &Pet::setName)
//    //    .def("getName", &Pet::getName)
//    //    .def_readwrite("name", &Pet::name);
//    m.def("add", &add, "A function which adds two numbers", py::arg("i") = 1, py::arg("j") = 2);
//    m.def("petname", &petname, "Creates pet with name and prints name", py::arg("name")="asdf");
//
//}

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
        .def_readwrite("maskThresh", &paramStruct::maskThresh);
    m.def("makeSingleMinkmap", &makeSingleMinkmap, "Reads file given by inname, creates Minkmap according to params, writes to outname", py::arg("inname"), py::arg("params"), py::arg("outname"));
    m.def("makeSequence", &makeSequence, "Reads file given by inname, creates threshold sequence of Minkmaps according to params, writes to outname", py::arg("inname"), py::arg("params"), py::arg("outname"));
}
