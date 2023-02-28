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
