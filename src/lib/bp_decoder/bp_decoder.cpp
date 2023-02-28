#include "bp_decoder.hpp"

#include <iostream>
#include <cmath>

namespace bp_decoder
{
    auto BpDecoder::getParityChecks(::sparse_matrix::Mod2SparseMatrix const &matrix, ::std::vector<uint8_t> const &code_word)
    {
        return matrix * code_word;
    }
    auto BpDecoder::init(::sparse_matrix::Mod2SparseMatrix &hx, ::std::vector<uint8_t> const &bit_error)
    {
        auto const prob_ratio_initial = this->error_prob / (1 - this->error_prob);
        for (auto const &row : hx.items_each_row)
        {
            for (auto const &item : row)
            {
                item->value.prob_rate = prob_ratio_initial;
                item->value.like_rate = 1;
            }
        }
        ::std::vector<double> log_prob_ratios(hx.col, 0.);
        ::std::vector<uint8_t> bit_syndrome{this->getParityChecks(hx, bit_error)};
        ::std::vector<uint8_t> decoding(hx.col, 0);
        return ::std::make_tuple(prob_ratio_initial, log_prob_ratios, bit_syndrome, decoding);
    }
    void BpDecoder::update(::sparse_matrix::Mod2SparseMatrix &matrix, ::std::vector<uint8_t> &decoding, ::std::vector<double> &log_prob_ratios, ::std::vector<uint8_t> const &syndrome, double prob_ratio_initial)
    {
        // Recompute likelihood ratios.
        double dl, t;
        for (auto i{0ULL}; i < matrix.row; i++)
        {
            dl = syndrome[i] ? -1 : 1;
            auto const &row = matrix.items_each_row[i];
            for (auto it{row.begin()}; it != row.end(); it++)
            {
                (*it)->value.like_rate = dl;
                dl *= 2 / (1 + (*it)->value.prob_rate) - 1;
            }
            dl = 1;
            for (auto it{row.rbegin()}; it != row.rend(); it++)
            {
                t = (*it)->value.like_rate * dl;
                (*it)->value.like_rate = (1 - t) / (1 + t);
                dl *= 2 / (1 + (*it)->value.prob_rate) - 1;
            }
        }

        // Recompute probability ratios.  Also find the next guess based on the individually most likely values.
        double pr;
        for (auto j{0ULL}; j < matrix.col; j++)
        {
            pr = prob_ratio_initial;
            auto const &col = matrix.items_each_col[j];
            for (auto it{col.begin()}; it != col.end(); it++)
            {
                (*it)->value.prob_rate = pr;
                pr *= (*it)->value.like_rate;
            }
            if (::std::isnan(pr))
                pr = 1;
            log_prob_ratios[j] = ::std::log(1 / pr);
            decoding[j] = pr >= 1;
            pr = 1;
            for (auto it{col.rbegin()}; it != col.rend(); it++)
            {
                (*it)->value.prob_rate *= pr;
                if (::std::isnan((*it)->value.prob_rate))
                    (*it)->value.prob_rate = 1;
                pr *= (*it)->value.like_rate;
            }
        }
    }
    size_t BpDecoder::hammingWeight(::std::vector<uint8_t> const &src1, ::std::vector<uint8_t> const &src2)
    {
        // assert(src1.size() == src2.size());
        auto size = src1.size();
        auto result{0ULL};
        for (auto i{0ULL}; i < size; i++)
            result += src1[i] ^ src2[i];
        return result;
    }
    BpDecoder::BpDecoder(Method method, double error_prob, int max_iter)
        : method{method}, error_prob{error_prob}, max_iter{max_iter} {}
    ::std::tuple<size_t, bool, std::vector<double>, std::vector<uint8_t>> BpDecoder::run(::sparse_matrix::Mod2SparseMatrix matrix, ::std::vector<uint8_t> const &bit_error)
    {
        // @todo
        if (this->method == Method::MIN_SUM)
            ::std::cerr << "MIN_SUM not implemented yet.";
        // setup
        auto &&[prob_ratio_initial, log_prob_ratios, bit_syndrome, decoding] = this->init(matrix, bit_error);
        auto best_hamming_weight{SIZE_MAX};
        ::std::vector<double> best_log_prob_ratios;
        auto has_decreased{false};
        // run
        for (auto it{0}; it < max_iter; it++)
        {
            this->update(matrix, decoding, log_prob_ratios, bit_syndrome, prob_ratio_initial);
            auto candidate_synd{this->getParityChecks(matrix, decoding)};
            auto hamming_weight = this->hammingWeight(bit_syndrome, candidate_synd);
            if (hamming_weight == 0)
                return std::make_tuple(it, true, log_prob_ratios, decoding);
            if (hamming_weight < best_hamming_weight)
            {
                best_hamming_weight = hamming_weight;
                best_log_prob_ratios = log_prob_ratios;
                has_decreased = true;
            }
        }
        if (has_decreased)
            return std::make_tuple(max_iter, false, best_log_prob_ratios, decoding);
        else
            return std::make_tuple(max_iter, false, log_prob_ratios, decoding);
    }
}