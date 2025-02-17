#ifndef FASTFEM_SKYLINEMATRIX_HPP
#define FASTFEM_SKYLINEMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <vector>
#include <memory>

namespace fastfem{
namespace linalg{

/**
 * @brief Stores a symmetric matrix in skyline format, storing only the lower triangular part
 */
class SkylineMatrix : public SparseMatrix {
protected:
    std::vector<double> values;
    std::shared_ptr<std::vector<size_t>> skyline;

public:
    SkylineMatrix(size_t n_cols, const std::vector<size_t>& skyline);

    Vector gemv(const Vector& x) const override;

    inline size_t nnz() const override { return skyline->back(); } 

    void add_entry(size_t index, double value);

    void print_pattern() const;

    inline bool is_symmetric() const override { return true; }

    void cholesky_factorize();

    Vector cholesky_solve(const Vector& b) const;

    void insert_entry(size_t i, size_t j, double value);

private:
    const double &get_entry(size_t i, size_t j) const override;
};

} // namespace linalg
} // namespace FastFem
    

#endif // FASTFEM_SKYLINEMATRIX_HPP