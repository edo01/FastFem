#include "FastFem/linalg/sparseMatrices/SkylineMatrix.hpp"

namespace fastfem{
namespace linalg{

SkylineMatrix::SkylineMatrix(size_t n_cols, const std::vector<size_t>& skyline) :
  SparseMatrix(skyline.size() - 1, n_cols),
    values(skyline.back()),
    skyline(std::make_shared<std::vector<size_t>>(skyline))
    {
        if(n_cols == 0 || skyline.size() < 2){
            throw std::invalid_argument("SkylineMatrix::SkylineMatrix(): invalid dimensions");
        }
    }

const double &SkylineMatrix::get_entry(size_t i, size_t j) const
{
    if (j > i) {
    std::swap(i, j);
    }

    size_t row_start = (*skyline)[i];
    size_t row_end = (*skyline)[i + 1];
    size_t row_length = row_end - row_start;

    // Compute the first column stored in row i
    size_t first_col_in_row = i - row_length + 1;

    // Check if A(i, j) is an implicit zero
    if (j < first_col_in_row) {
        static double zero = 0.0;
        return zero;
    }

    // Find position of A(i, j) in values array
    size_t position = row_start + (j - first_col_in_row);
    return values[position]; 
}

Vector SkylineMatrix::gemv(const Vector& x) const
{
    size_t row_start, row_end, row_length, first_col_in_row, col_index;

    if (n_cols != x.size()) {
        throw std::invalid_argument("SkylineMatrix::gemv(): incompatible dimensions");
    }

    Vector y(n_rows);

    for (size_t i = 0; i < n_rows; ++i)
    {
        row_start = skyline->at(i);
        row_end = skyline->at(i + 1);
        row_length = row_end - row_start;

        // First column stored in row i
        first_col_in_row = i - row_length + 1;

        col_index = first_col_in_row;
        for (size_t k = row_start; k < row_end; ++k)
        {
            y[i] += values[k] * x[col_index];

            // Since the matrix is symmetric, update y[col_index] (except for diagonal elements)
            if (i != col_index) {
                y[col_index] += values[k] * x[i];
            }

            ++col_index;
        }
    }

    return y;
}

void SkylineMatrix::add_entry(size_t index, double value)
{
    if (index >= skyline->back()) {
        throw std::out_of_range("SkylineMatrix::add_entry(): index out of range");
    }
    values[index] = value;
}

void SkylineMatrix::print_pattern() const
{
    for (size_t i = 0; i < n_rows; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            bool found = false;
            for (size_t k = (*skyline)[i]; k < (*skyline)[i + 1]; ++k)
            {
                if (j == i - ((*skyline)[i + 1] - k - 1)) 
                {
                    std::cout << "1 ";
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;
    }
}


} // namespace linalg
} // namespace FastFem