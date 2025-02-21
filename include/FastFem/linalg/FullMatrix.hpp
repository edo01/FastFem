#ifndef FULLMATRIX_HPP
#define FULLMATRIX_HPP

#include <vector>
#include "FastFem/types/CommonTypes.hpp"

namespace fastfem{
namespace linalg{

using fastfem::types::ff_index;

class FullMatrix
{
public:
    FullMatrix(ff_index n_rows, ff_index n_cols) : n_rows(n_rows), n_cols(n_cols), data(n_rows * n_cols) {}
    FullMatrix(ff_index dim) : FullMatrix(dim, dim) {}
    inline const double &operator()(ff_index i, ff_index j) const { return data[i * n_cols + j]; }
    inline double &operator()(ff_index i, ff_index j) { return data[i * n_cols + j]; }

    inline void set_to_zero() { std::fill(data.begin(), data.end(), 0.0); }

    inline ff_index get_n_rows() const { return n_rows; }
    inline ff_index get_n_cols() const { return n_cols; }
private:
    ff_index n_rows;
    ff_index n_cols;
    std::vector<double> data;
};

} // namespace linalg
} // namespace FastFem

#endif // FULLMATRIX_HPP