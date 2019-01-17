#include <pybind11/pybind11.h>

#include "kmerizer.h"


PYBIND11_MODULE(kmerizer, m) {
    m.doc() = "Fast C++ based k-mer counter.";

    py::class_<strainge::kmerizer>(m, "kmerize_iter")
        .def(py::init<int, std::string const&>())
        .def("__iter__",
                [](strainge::kmerizer const& obj) { return py::make_iterator(obj.begin(), obj.end()); },
                py::keep_alive<0, 1>());

    m.def("kmerize", &strainge::kmerize,
            "Return a NumPy array with k-mers for a given sequence",
            py::arg("k"), py::arg("sequence"));
    m.def("kmerize_into_array", &strainge::kmerize_into_array,
            "Kmerize a sequence and store k-mers in a pre-allocated NumPy array",
            py::arg("k"), py::arg("sequence"), py::arg("out_array"), py::arg("offset"));
    m.def("merge_counts", &strainge::merge_counts,
            "Merge and sum two k-mer sets and their count arrays.",
            py::arg("kmers1"), py::arg("counts1"), py::arg("kmers2"), py::arg("counts2"));
    m.def("count_common", &strainge::count_common,
            "Count the number of common k-mers between two sets.",
            py::arg("kmers1"), py::arg("kmers2"));
    m.def("build_kmer_count_matrix", &strainge::build_kmer_count_matrix,
            "Build a matrix with all k-mer counts combined for the given list of k-mer sets. "
            "Input should be a list of (kmers, counts) tuples.",
            py::arg("kmerset_list"));

    m.def("intersect", &strainge::intersect,
            "Return the intersection of two k-mer sets.",
            py::arg("kmers1"), py::arg("kmers2"));
    m.def("intersect_ix", &strainge::intersect_ix,
            "Return a NumPy index array. This array has value True at position i "
            "if the k-mer at kmers1[i] is also in kmers2. This can be used to index "
            "other NumPy arrays the same size as kmers1.",
            py::arg("kmers1"), py::arg("kmers2"));

    m.def("diff", &strainge::diff,
            "Return the difference of two k-mer sets (kmers1 minus kmers2)",
            py::arg("kmers1"), py::arg("kmers2"));

    m.def("kmerset_in_product", &strainge::kmerset_in_product,
            "Calculate the in-product between the count vectors of two kmer-"
            "sets.",
            py::arg("kmers1"), py::arg("counts1"), py::arg("kmers2"), py::arg("counts2"));

    m.def("fnvhash_kmers", &strainge::fnvhash_kmers,
            "Hash all values in the k-mer set using FNV hash function.",
            py::arg("k"), py::arg("kmers"));

    py::register_exception<strainge::KmerizeError>(m, "KmerizeError");
}
