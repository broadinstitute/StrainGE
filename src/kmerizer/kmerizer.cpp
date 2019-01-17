/**
 * Simple k-mer counter for Python
 */

#include <vector>
#include <set>
#include "kmerizer.h"

using std::vector;
using std::set;
using std::size_t;

namespace strainge {

    static vector<kmer_t> kmerize_internal(int k, const std::string& sequence) {
        vector<kmer_t> kmers;

        int shift = 2 * (k - 1);
        kmer_t mask = (k < 32) ? ((kmer_t) 1 << (2 * k)) - 1 : -1;

        int n = 0;
        kmer_t fw = 0;
        kmer_t reverse = 0;

        for(char b : sequence) {
            b = toupper(b);

            Base value;
            switch(b) {
                case 'A':
                    value = Base::A;
                    break;
                case 'C':
                    value = Base::C;
                    break;
                case 'G':
                    value = Base::G;
                    break;
                case 'T':
                    value = Base::T;
                    break;
                default:
                    fw = reverse = n = 0;
                    continue;
            }

            fw = ((fw << 2) & mask) | static_cast<uint64_t>(value);
            reverse = ((reverse >> 2) & mask) | (static_cast<uint64_t>(rc(value)) << shift);

            if(++n >= k) {
                auto kmer = fw < reverse ? fw : reverse;
                kmers.push_back(kmer);
            }
        }

        return kmers;
    }

    kmerset_t kmerize(int k, const std::string& sequence) {
        check_k(k);

        auto kmers = kmerize_internal(k, sequence);

        // Create NumPy array from found k-mers (copies data)
        kmerset_t kmerset = kmerset_t(py::buffer_info(
                kmers.data(), sizeof(kmer_t),
                py::format_descriptor<kmer_t>::format(),
                1, // dimensions
                { kmers.size() },
                { sizeof(kmer_t) }
        ));

        return kmerset;
    }


    /**
     * This function k-merizes the given sequence and stores the k-mers in
     * a pre-allocated NumPy array.
     *
     * @param k k-mer size, can be at most 32
     * @param sequence The sequence to k-merize
     * @param array The pre-allocated numpy array (dtype=uint64)
     * @param offset The position in the array to start storing k-mers
     *
     * @return Number of k-mers stored
     */
    size_t kmerize_into_array(int k, const std::string& sequence,
            kmerset_t& out_array, unsigned int offset) {
        check_k(k);

        auto kmers = kmerize_internal(k, sequence);

        if(kmers.size() + offset > out_array.shape(0)) {
            throw KmerizeError("Number of kmers exceeds space available in NumPy array");
        }

        auto proxy = out_array.mutable_unchecked<1>();

        for(kmer_t kmer : kmers) {
            proxy(offset) = kmer;
            ++offset;
        }

        return kmers.size();
    }

    size_t count_common(const kmerset_t& kmers1,
            const kmerset_t& kmers2) {
        size_t common = 0;
        size_t size1 = kmers1.shape(0);
        size_t size2 = kmers2.shape(0);

        auto proxy1 = kmers1.unchecked<1>();
        auto proxy2 = kmers2.unchecked<1>();

        // Arrays should be sorted, walk through them in parallel
        for(size_t i1 = 0, i2 = 0; i1 < size1 && i2 < size2;) {
            kmer_t kmer1 = proxy1(i1);
            kmer_t kmer2 = proxy2(i2);

            if(kmer1 == kmer2) {
                // In both sets
                ++common;

                ++i1;
                ++i2;
            } else if(kmer1 < kmer2) {
                // First set lags
                ++i1;
            } else {
                ++i2;
            }
        }

        return common;
    }

    std::tuple<kmerset_t, kmercounts_t> merge_counts(
            const kmerset_t& kmers1,
            const kmercounts_t& counts1,
            const kmerset_t& kmers2,
            const kmercounts_t& counts2) {
        size_t size1 = kmers1.shape(0);
        size_t size2 = kmers2.shape(0);
        size_t new_size = size1 + size2 - count_common(kmers1, kmers2);

        kmerset_t new_set(new_size);
        kmercounts_t new_counts(new_size);

        // Direct access proxies
        auto proxy1 = kmers1.unchecked<1>();
        auto proxy2 = kmers2.unchecked<1>();
        auto proxyc1 = counts1.unchecked<1>();
        auto proxyc2 = counts2.unchecked<1>();

        auto proxy_new = new_set.mutable_unchecked<1>();
        auto proxy_counts = new_counts.mutable_unchecked<1>();

        size_t kcount = 0;
        size_t i1, i2;
        for(i1 = 0, i2 = 0; i1 < size1 && i2 < size2;) {
            kmer_t kmer1 = proxy1(i1);
            kmer_t kmer2 = proxy2(i2);

            if(kmer1 == kmer2) {
                // In both sets
                proxy_new(kcount) = kmer1;
                proxy_counts(kcount) = proxyc1(i1) + proxyc2(i2);
                ++i1;
                ++i2;
            } else if (kmer1 < kmer2) {
                // Only in set 1
                proxy_new(kcount) = kmer1;
                proxy_counts(kcount) = proxyc1(i1);
                ++i1;
            } else {
                proxy_new(kcount) = kmer2;
                proxy_counts(kcount) = proxyc2(i2);
                ++i2;
            }
            ++kcount;
        }

        // Check for leftovers
        while(i1 < size1) {
            proxy_new(kcount) = proxy1(i1);
            proxy_counts(kcount) = proxyc1(i1);
            ++kcount;
            ++i1;
        }
        while(i2 < size2) {
            proxy_new(kcount) = proxy2(i2);
            proxy_counts(kcount) = proxyc2(i2);
            ++kcount;
            ++i2;
        }

        return std::make_tuple(new_set, new_counts);
    }

    std::tuple<vector<kmer_t>, py::array_t<uint64_t>> build_kmer_count_matrix(
            const std::vector<kmers_with_counts_t>& kmersets) {
        // Sorted set, because we want our output matrix to be sorted too.
        set<kmer_t> all_unique_kmers;

        for(auto& elem : kmersets) {
            const kmerset_t& kmers = std::get<0>(elem);

            auto proxy = kmers.unchecked<1>();
            for(int i = 0; i < kmers.size(); ++i) {
                all_unique_kmers.insert(proxy(i));
            }
        }

        py::array_t<uint64_t> kmer_matrix({all_unique_kmers.size(), kmersets.size()});

        // Keep track where we are for each individual kmerset
        vector<size_t> positions;
        for(size_t i = 0; i < kmersets.size(); ++i) {
            positions.push_back(0);
        }

        // Loop over all k-mers and and see which of the kmer sets contain this
        // k-mer. Add the corresponding counts to the matrix
        auto proxy = kmer_matrix.mutable_unchecked<2>();
        int i = 0;
        for(auto kmer : all_unique_kmers) {
            int j = 0;
            for(auto& elem : kmersets) {
                const kmerset_t& kmers = std::get<0>(elem);
                const kmercounts_t& counts = std::get<1>(elem);

                auto set_proxy = kmers.unchecked<1>();
                auto counts_proxy = counts.unchecked<1>();
                size_t set_pos = positions[j];

                if(set_proxy(set_pos) == kmer) {
                    proxy(i, j) = counts_proxy(set_pos);
                    ++positions[j];
                } else {
                    proxy(i, j) = 0;
                }

                ++j;
            }
            ++i;
        }

        // Create a vector of all unique k-mers.
        // We do this instead of returning the set itself because the
        // ordered C++ set gets converted to a unordered Python set, which means
        // we lose the k-mer order.
        vector<kmer_t> kmer_list(all_unique_kmers.begin(), all_unique_kmers.end());

        return std::make_tuple(kmer_list, kmer_matrix);
    }


    kmerset_t intersect(const kmerset_t& kmers1,
            const kmerset_t& kmers2) {
        unsigned int common = count_common(kmers1, kmers2);
        kmerset_t new_set(common);

        size_t size1 = kmers1.shape(0);
        size_t size2 = kmers2.shape(0);

        auto proxy1 = kmers1.unchecked<1>();
        auto proxy2 = kmers2.unchecked<1>();
        auto proxy_new = new_set.mutable_unchecked<1>();

        // Arrays should be sorted, walk through them in parallel
        size_t kcount = 0;
        for(size_t i1 = 0, i2 = 0; i1 < size1 && i2 < size2;) {
            kmer_t kmer1 = proxy1(i1);
            kmer_t kmer2 = proxy2(i2);

            if(kmer1 == kmer2) {
                proxy_new(kcount) = kmer1;
                ++i1;
                ++i2;
                ++kcount;
            } else if(kmer1 < kmer2) {
                ++i1;
            } else {
                ++i2;
            }
        }

        return new_set;
    }

    py::array_t<bool> intersect_ix(const kmerset_t& kmers1,
            const kmerset_t& kmers2) {
        size_t size1 = kmers1.shape(0);
        size_t size2 = kmers2.shape(0);

        py::array_t<bool> new_set(size1);

        auto proxy1 = kmers1.unchecked<1>();
        auto proxy2 = kmers2.unchecked<1>();
        auto proxy_new = new_set.mutable_unchecked<1>();

        // Set all to false by default
        for(size_t i = 0; i < size1; ++i) {
            proxy_new(i) = false;
        }

        for(size_t i1 = 0, i2 = 0; i1 < size1 && i2 < size2;) {
            kmer_t kmer1 = proxy1(i1);
            kmer_t kmer2 = proxy2(i2);
            if(kmer1 == kmer2) {
                proxy_new(i1) = true;

                ++i1;
                ++i2;
            } else if (kmer1 < kmer2) {
                // Only in kmers1
                ++i1;
            } else if(kmer2 < kmer1) {
                // Only in kmers2
                ++i2;
            }
        }

        return new_set;
    }

    kmerset_t diff(const kmerset_t& kmers1, const kmerset_t& kmers2) {
        size_t size1 = kmers1.shape(0);
        size_t size2 = kmers2.shape(0);
        size_t common = count_common(kmers1, kmers2);

        size_t new_size = size1 - common;
        kmerset_t new_set(new_size);

        // Direct access proxies
        auto proxy1 = kmers1.unchecked<1>();
        auto proxy2 = kmers2.unchecked<1>();
        auto proxy_new = new_set.mutable_unchecked<1>();

        size_t kcount = 0, i1, i2;
        for(i1 = 0, i2 = 0; i1 < size1 && i2 < size2; ) {
            kmer_t kmer1 = proxy1(i1);
            kmer_t kmer2 = proxy2(i2);

            if(kmer1 == kmer2) {
                // In both
                ++i1;
                ++i2;
            } else if(kmer1 < kmer2) {
                // Only in kmers1
                proxy_new(kcount) = kmer1;
                ++kcount;
                ++i1;
            } else {
                ++i2;
            }
        }

        // check for leftovers
        while(i1 < size1) {
            proxy_new(kcount) = proxy1(i1);
            ++kcount;
            ++i1;
        }

        return new_set;
    }

    uint64_t kmerset_in_product(
            const kmerset_t& kmers1,
            const kmercounts_t& counts1,
            const kmerset_t& kmers2,
            const kmercounts_t& counts2) {

        size_t size1 = kmers1.shape(0);
        size_t size2 = kmers2.shape(0);

        auto proxy1 = kmers1.unchecked<1>();
        auto proxy2 = kmers2.unchecked<1>();

        auto counts_proxy1 = counts1.unchecked<1>();
        auto counts_proxy2 = counts2.unchecked<1>();

        uint64_t product = 0;
        for(size_t i1 = 0, i2 = 0; i1 < size1 && i2 < size2;) {
            kmer_t kmer1 = proxy1(i1);
            kmer_t kmer2 = proxy2(i2);

            if(kmer1 == kmer2) {
                product += counts_proxy1(i1) * counts_proxy2(i2);
                ++i1;
                ++i2;
            } else if(kmer1 < kmer2) {
                ++i1;
            } else if(kmer2 < kmer1) {
                ++i2;
            }
        }

        return product;
    }

    const uint64_t FNV_PRIME = 1099511628211u;
    const uint64_t FNV_OFFSET = 14695981039346656037u;

    kmerset_t fnvhash_kmers(int k, const kmerset_t& kmers) {
        check_k(k);

        size_t size = kmers.shape(0);
        kmerset_t hashed(size);

        int kbits = 2 * k;

        // Direct access proxies
        auto proxy = kmers.unchecked<1>();
        auto proxy_hashed = hashed.mutable_unchecked<1>();

        for(size_t i = 0; i < size; ++i) {
            kmer_t kmer = proxy(i);
            kmer_t hash = FNV_OFFSET;

            for(int b = kbits; b > 0; b -= 8) {
                hash ^= (kmer & 0xFF);
                hash *= FNV_PRIME;
                kmer >>= 8;
            }

            proxy_hashed(i) = hash;
        }

        return hashed;
    }
}
