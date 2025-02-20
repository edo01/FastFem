#ifndef FASTFEM_SYMCSRMATRIX_HPP
#define FASTFEM_SYMCSRMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>

namespace fastfem{
namespace linalg{

using types::ff_index;

/**
 * @brief Symmetric CSR matrix exploting a symmetric CSRPattern and storing only the upper triangular part
 */
class SymCSRMatrix : public CSRMatrix
{
public:
    SymCSRMatrix(ff_index n_cols, const CSRPattern& pattern);
    SymCSRMatrix(const CSRMatrix& A) : CSRMatrix(A) {}

    Vector gemv(const Vector& x) const override;
    void accumulate_entry(ff_index i, ff_index j, double value) override;

    inline bool is_symmetric() const override { return true; }
    void operator=(const double& value) override;

private:
    const double &get_entry(ff_index i, ff_index j) const override;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_SYMCSRMATRIX_HPP