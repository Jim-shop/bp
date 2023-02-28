#pragma once

#ifndef _BP_DECODER_HPP_
#define _BP_DECODER_HPP_

#include "sparse_matrix.hpp"

namespace bp_decoder
{
    class BpDecoder
    {
    public: // types
        enum class Method
        {
            MIN_SUM,
            PRODUCT_SUM
        };

    private: // vars
        Method method;
        double error_prob;
        int max_iter;

    private: // utils
        auto getParityChecks(::sparse_matrix::Mod2SparseMatrix const &matrix, ::std::vector<uint8_t> const &code_word);
        auto init(::sparse_matrix::Mod2SparseMatrix &matrix, ::std::vector<uint8_t> const &bit_error);
        void update(::sparse_matrix::Mod2SparseMatrix &matrix, ::std::vector<uint8_t> &decoding, ::std::vector<double> &log_prob_ratios, ::std::vector<uint8_t> const &syndrome, double prob_ratio_initial);
        size_t hammingWeight(::std::vector<uint8_t> const &src1, ::std::vector<uint8_t> const &src2);

    public: // apis
        BpDecoder(Method method, double error_prob, int max_iter);
        ::std::tuple<size_t, bool, std::vector<double>, std::vector<uint8_t>> run(::sparse_matrix::Mod2SparseMatrix matrix, ::std::vector<uint8_t> const &bit_error);
    };
}

#endif