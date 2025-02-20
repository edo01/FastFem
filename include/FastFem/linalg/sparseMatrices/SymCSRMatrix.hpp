#ifndef FASTFEM_SYMCSRMATRIX_HPP
#define FASTFEM_SYMCSRMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>

namespace fastfem{
namespace linalg{

/**
 * @brief Symmetric CSR matrix exploting a symmetric CSRPattern and storing only the upper triangular part
 */
class SymCSRMatrix : public CSRMatrix
{
public:
    SymCSRMatrix(size_t n_cols, const CSRPattern& pattern);
    SymCSRMatrix(const CSRMatrix& A) : CSRMatrix(A) {}

    Vector gemv(const Vector& x) const override;
    void accumulate_entry(size_t i, size_t j, double value) override;

    inline bool is_symmetric() const override { return true; }

private:
    const double &get_entry(size_t i, size_t j) const override;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_SYMCSRMATRIX_HPP