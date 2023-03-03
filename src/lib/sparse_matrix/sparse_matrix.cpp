/** Copyright (c) 2023 Jim-shop
 * bp is licensed under Mulan PubL v2.
 * You can use this software according to the terms and conditions of the Mulan PubL v2.
 * You may obtain a copy of Mulan PubL v2 at:
 *          http://license.coscl.org.cn/MulanPubL-2.0
 * THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
 * EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
 * MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
 * See the Mulan PubL v2 for more details.
 */
#include "sparse_matrix.hpp"

namespace sparse_matrix
{
    ::std::vector<uint8_t> Mod2SparseMatrix::operator*(::std::vector<uint8_t> const &vec) const
    {
        if (this->col != vec.size())
            throw ::std::runtime_error("Vec length mismatch matrix col."s);
        ::std::vector<uint8_t> result(this->row, 0);
        for (auto j{0}; j < this->col; j++)
            if (vec[j])
                for (auto const &item : this->items_each_col[j])
                    result[item->row_index] ^= 1;
        return result;
    }
}
