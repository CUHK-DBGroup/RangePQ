#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // To convert list of list <-> vec<vec<int>> for e.g. posting_lists
#include "segr.h"
#include "BST.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

namespace segr {
PYBIND11_MODULE(main, m) {
    py::class_<SegrCpp>(m, "SegrCpp")
        .def(py::init<>())  // required in pickle
        .def(py::init<py::array_t<float>, bool>())
        .def("reconfigure", &SegrCpp::Reconfigure)
        .def("reconfigure_2", &SegrCpp::Reconfigure_2)
        .def("add_codes", &SegrCpp::AddCodes)
        .def("delete_key", &SegrCpp::DeleteKey)
        .def("memory", &SegrCpp::Memory)
        .def("query_linear", &SegrCpp::QueryLinear,
             py::arg("query").noconvert(),  // Prohibit implicit data conversion
             py::arg("topk"),
             py::arg("target_ids").noconvert()  // Prohibit implicit data conversion
             )
        .def("query_seg_range", &SegrCpp::QuerySegRange,
            py::arg("query").noconvert(),  // Prohibit implicit data conversion
            py::arg("topk"),py::arg("range_l"),py::arg("range_r"),py::arg("L")
            )
        .def("clear", &SegrCpp::Clear)
        .def_readwrite("verbose", &SegrCpp::verbose_)
        .def_readonly("coarse_centers", &SegrCpp::coarse_centers_)
        .def_readonly("flattened_codes", &SegrCpp::flattened_codes_)
        .def_readonly("postingList", &SegrCpp::postingList)
        .def_property_readonly("N", &SegrCpp::GetN)
        .def_property_readonly("nlist", &SegrCpp::GetNumList)
        .def(py::pickle(
            [](const SegrCpp &p){
                return py::make_tuple(p.codewords_, p.verbose_,
                p.coarse_centers_, p.flattened_codes_, p.postingList);
            },
            [](py::tuple t){
                if (t.size() != 5) {
                    throw std::runtime_error("Invalid state when reading pickled item");
                }
                SegrCpp p;
                p.codewords_ = t[0].cast<std::vector<std::vector<std::vector<float>>>>();
                p.M_ = p.codewords_.size();
                p.Ks_ = p.codewords_[0].size();
                p.verbose_ = t[1].cast<bool>();
                p.coarse_centers_ = t[2].cast<std::vector<std::vector<unsigned char>>>();
                p.flattened_codes_ = t[3].cast<std::vector<unsigned char>>();
                p.postingList = t[4].cast<std::vector<BST*>>();
                return p;
            }
        ));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
} // namespace rii