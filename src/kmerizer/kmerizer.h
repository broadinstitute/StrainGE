/**
 * Simple k-mer counter for Python
 */

#include <tuple>
#include <string>
#include <iterator>
#include <exception>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <cstdint>

#ifndef KMERIZER_H
#define KMERIZER_H

namespace py = pybind11;

namespace strainge {
    enum class Base : uint64_t {
        A = 0,
        C = 1,
        G = 2,
        T = 3
    };

    /**
     * Reverse complement of a base
     */
    constexpr Base rc(Base base) {
        return static_cast<Base>(static_cast<uint64_t>(base) ^ 3);
    }

    typedef uint64_t kmer_t;
    typedef py::array_t<kmer_t> kmerset_t;
    typedef py::array_t<uint64_t> kmercounts_t;
    typedef std::tuple<kmerset_t, kmercounts_t> kmers_with_counts_t;

    /**
     * Iterator class to obtain all k-mers in a string without the requirement
     * to keep them in memory.
     */
    class kmerizer {
        public:
            template<typename T>
            class BaseKmerIterator {
                public:
                    using value_type = T;
                    using difference_type = std::ptrdiff_t;
                    using pointer = T*;
                    using reference = T&;
                    using iterator_category = std::forward_iterator_tag;

                    BaseKmerIterator() : k(0), fw(0), rev(0), n(0), shift(0), mask(0) { }

                    BaseKmerIterator(int k, std::string const& sequence) :
                            k(k), fw(0), rev(0), n(0),
                            shift(2 * (k - 1)),
                            mask((k < 32) ? ((kmer_t) 1 << (2 * k)) - 1 : -1),
                            pos(sequence.begin()),
                            end(sequence.end()) {

                        this->next_kmer();
                    }

                    BaseKmerIterator(int k, std::string::const_iterator const& end) :
                            k(k), fw(0), rev(0), n(0),
                            shift(2 * (k - 1)),
                            mask((k < 32) ? ((kmer_t) 1 << (2 * k)) - 1 : -1),
                            pos(end),
                            end(end) {
                    }

                    value_type operator*() const {
                        return (this->fw < this->rev) ? this->fw : this->rev;
                    }

                    // Pre-increment
                    BaseKmerIterator& operator++() {
                        this->next_kmer();

                        return *this;
                    }

                    // Post increment
                    BaseKmerIterator operator++(int) {
                        BaseKmerIterator tmp = *this;
                        this->next_kmer();

                        return tmp;
                    }

                    bool operator==(BaseKmerIterator const& rhs) {
                        return (this->pos == rhs.pos &&
                                this->fw == rhs.fw &&
                                this->rev == rhs.rev);
                    }

                    bool operator!=(BaseKmerIterator const& rhs) {
                        return !(*this == rhs);
                    }


                private:
                    int k;
                    kmer_t fw;
                    kmer_t rev;
                    int n;

                    int const shift;
                    int const mask;

                    std::string::const_iterator pos;
                    std::string::const_iterator end;

                    void next_kmer() {
                        if(this->pos == this->end) {
                            fw = rev = n = 0;
                            return;
                        }

                        do {
                            char const b = toupper(*pos);
                            ++this->pos;

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
                                    fw = rev = n = 0;
                                    continue;
                            }

                            this->fw = ((fw << 2) & mask) | static_cast<uint64_t>(value);
                            this->rev = ((rev >> 2) & mask) | (static_cast<uint64_t>(rc(value)) << shift);

                            if(this->n < k) {
                                ++this->n;
                            }
                        } while(this->n < k && this->pos != this->end);
                    }
            };

            using iterator = BaseKmerIterator<kmer_t>;
            using const_iterator = BaseKmerIterator<const kmer_t>;

            kmerizer(int k, std::string const& sequence) : k(k), sequence(sequence) { }

            iterator begin() {
                return iterator(k, sequence);
            }

            iterator end() {
                return iterator(k, sequence.end());
            }

            const_iterator begin() const {
                return const_iterator(k, sequence);
            }

            const_iterator end() const {
                return const_iterator(k, sequence.end());
            }

        private:
            int k;
            std::string sequence;
    };

    class KmerizeError : public std::runtime_error {
        public:
            KmerizeError(const std::string& msg) : std::runtime_error(msg) { }
    };


    inline void check_k(int k) {
        if(k < 1 || k > 32) {
            throw KmerizeError("k is out of range, must be in range [1, 32]");
        }
    }

    /**
     * This function returns a NumPy array of k-mers of the given sequence.
     * k-mers are encoded as a uint64_t.
     *
     * @param k k-mer size, can be at most 32
     * @param sequence The sequence to kmerize
     *
     * @return NumPy array of size len(sequence) + 1 - k with each k-mer encoded
     *     as a uint64_t
     */
    kmerset_t kmerize(int k, const std::string& sequence);

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
            kmerset_t& array, unsigned int offset);


    /**
     * Count the number of common k-mers between two k-mer sets. The k-mer sets
     * should be sorted.
     *
     * @param kmers1 The first set of k-mers as NumPy array
     * @param kmers2 The second set of k-mers as NumPy array
     *
     * @return Number of common k-mers
     */
    size_t count_common(const kmerset_t& kmers1,
            const kmerset_t& kmers2);

    /**
     * Merge two k-mer sets and their corresponding counts.
     *
     * @return A tuple with a new NumPy array containing k-mers from both sets
     *    and a separate NumPy array with corresponding updated counts
     */
    std::tuple<kmerset_t, kmercounts_t> merge_counts(
            const kmerset_t& kmers1,
            const kmercounts_t& counts1,
            const kmerset_t& kmers2,
            const kmercounts_t& counts2
    );

    /**
     * Build a k-mer count matrix for a given list of kmersets.
     *
     * This function returns a matrix with all k-mers and corresponding counts
     * combined for a given list of kmer sets. This is useful when you want to
     * use this data with scikit learn or other machine learning libraries.
     *
     * The input is a vector of tuples (kmers, counts).
     *
     * @return A tuple with as first entry a list with the actual kmers,
     *         i.e. the labels for the rows, and as second entry a
     *         KxL NumPy array with k-mer counts. K: number of unique k-mers
     *         in given list of k-mer sets, L: number of k-mer sets given. The
     *         k-mers will be sorted.
     */
    std::tuple<std::vector<kmer_t>, py::array_t<uint64_t>> build_kmer_count_matrix(
            const std::vector<kmers_with_counts_t>& kmersets);

    /**
     * Return the intersection between two k-mer sets.
     *
     * @return A new numpy array containing only k-mers that are in both sets.
     */
    kmerset_t intersect(const kmerset_t& kmers1,
            const kmerset_t& kmers2);

    /**
     * Return a numpy array containing the value True on position i if the k-mer
     * at kmers1[i] is also in kmers2. This can be used to index other NumPy
     * arrays of the same size as kmers1.
     *
     * @return A NumPy array that can be used to index other numpy arrays.
     */
    py::array_t<bool> intersect_ix(const kmerset_t& kmers1,
            const kmerset_t& kmers2);

    /**
     * Return the difference between two k-mer sets.
     *
     * @return NumPy array containing k-mers that are in k-mer set 1, but not in
     *     k-mer set 2.
     */
    kmerset_t diff(const kmerset_t& kmers1, const kmerset_t& kmers2);

    /**
     * Calculate the in-product between two k-mer sets.
     */
    uint64_t kmerset_in_product(
            const kmerset_t& kmers1,
            const kmercounts_t& counts1,
            const kmerset_t& kmers2,
            const kmercounts_t& counts2
    );

    /**
     * Hash each k-mer in a k-mer set using Fowler-Voll-No hash function. Useful
     * when calculating min-hash of a k-mer set.
     *
     * @param k k-mer size
     * @param kmers k-mer set to hash
     *
     * @return A NumPy array with each k-mer hash value
     */
    kmerset_t fnvhash_kmers(int k, const kmerset_t& kmers);
}

#endif
