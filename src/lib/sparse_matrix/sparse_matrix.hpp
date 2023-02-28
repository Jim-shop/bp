#pragma once

#ifndef _SPARSE_MATRIX_HPP_
#define _SPARSE_MATRIX_HPP_

#include <vector>
#include <memory>
#include <iostream>
#include <string>
#include <algorithm>

using ::std::operator""s;
using ::std::operator""sv;

namespace sparse_matrix
{
    template <typename T>
    class SparseMatrix
    {
    private: // types
        struct Item
        {
            size_t row_index, col_index;
            T value;
        };

    public: // members
        size_t row, col;
        ::std::vector<::std::vector<::std::shared_ptr<Item>>> items_each_row, items_each_col;

    public: // apis
        SparseMatrix() : row{0}, col{0} {}
        SparseMatrix(SparseMatrix const &that)
        {
            this->row = that.row;
            this->col = that.col;
            this->items_each_row.resize(this->row);
            this->items_each_col.resize(this->col);
            for (auto i{0ULL}; i < this->row; i++)
            {
                for (auto const &item : that.items_each_row[i])
                {
                    auto copy_item{::std::make_shared<Item>(*item)};
                    this->items_each_row[i].push_back(copy_item);
                    this->items_each_col[item->col_index].push_back(copy_item);
                }
            }
        }
        // Read alist file. See http://www.inference.org.uk/mackay/codes/alist.html
        friend ::std::istream &operator>>(::std::istream &stream, SparseMatrix &me)
        {
            if (!stream)
                throw ::std::runtime_error("Could not open alist file."s);
            stream >> me.row >> me.col;

            size_t max_weight_forall_rows, max_weight_forall_cols;
            stream >> max_weight_forall_rows >> max_weight_forall_cols;
            if (max_weight_forall_rows > me.row || max_weight_forall_cols > me.col)
                throw ::std::runtime_error("Unexpected EOF or invalid max_weight."s);

            ::std::vector<size_t> weights_each_row(me.row), weights_each_col(me.col);
            ::std::vector<::std::vector<size_t>> indexs_each_row, indexs_each_col;
            for (auto &row_weight : weights_each_row)
            {
                stream >> row_weight;
                if (row_weight > max_weight_forall_rows)
                    throw ::std::runtime_error("Unexpected EOF or invalid row weight."s);
                indexs_each_row.emplace_back(::std::vector<size_t>(row_weight));
            }
            for (auto &col_weight : weights_each_col)
            {
                stream >> col_weight;
                if (col_weight > max_weight_forall_cols)
                    throw ::std::runtime_error("Unexpected EOF or invalid col weight."s);
                indexs_each_col.emplace_back(::std::vector<size_t>(col_weight));
            }
            for (auto &row : indexs_each_row)
            {
                for (auto &index : row)
                {
                    stream >> index; // start from 1
                    index--;
                    if (index >= me.col)
                        throw ::std::runtime_error("Unexpected EOF or invalid index."s);
                }
            }
            for (auto &col : indexs_each_col)
            {
                for (auto &index : col)
                {
                    stream >> index;
                    index--;
                    if (index >= me.row)
                        throw ::std::runtime_error("Unexpected EOF or invalid index."s);
                }
            }

            me.items_each_row.resize(me.row);
            me.items_each_col.resize(me.col);
            for (size_t i{0}; i < me.row; i++)
            {
                for (size_t j : indexs_each_row[i])
                {
                    auto item_ptr{::std::make_shared<Item>(i, j)};
                    me.items_each_row[i].push_back(item_ptr);
                    me.items_each_col[j].push_back(item_ptr);
                }
            }
            return stream;
        }
        friend ::std::ostream &operator<<(::std::ostream &stream, SparseMatrix const &me)
        {
            ::std::vector<size_t> weights_each_row, weights_each_col;
            for (auto const &row : me.items_each_row)
                weights_each_row.push_back(row.size());
            for (auto const &col : me.items_each_col)
                weights_each_col.push_back(col.size());
            auto max_weight_forall_rows = *::std::max_element(weights_each_row.begin(), weights_each_row.end()),
                 max_weight_forall_cols = *::std::max_element(weights_each_col.begin(), weights_each_col.end());

            stream << me.row << ' ' << me.col << '\n'
                   << max_weight_forall_rows << ' ' << max_weight_forall_cols << '\n';
            for (auto row_weight : weights_each_row)
                stream << row_weight << ' ';
            stream << '\n';
            for (auto col_weight : weights_each_col)
                stream << col_weight << ' ';
            stream << '\n';
            for (auto const &row : me.items_each_row)
            {
                for (auto const &item : row)
                    stream << (item->col_index + 1) << ' ';
                stream << '\n';
            }
            for (auto const &col : me.items_each_col)
            {
                for (auto const &item : col)
                    stream << (item->row_index + 1) << ' ';
                stream << '\n';
            }
            return stream;
        }
    };
    struct Prob
    {
        double prob_rate, like_rate;
    };
    class Mod2SparseMatrix : public SparseMatrix<Prob>
    {
    public: // extra api
        ::std::vector<uint8_t> operator*(::std::vector<uint8_t> const &vec) const;
    };
}

#endif
