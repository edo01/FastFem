#ifndef SYMCSRMATRIX_HPP
#define SYMCSRMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>
// #include <vector>
// #include <memory>

namespace fastfem{
namespace linalg{

class SymCSRMatrix : public CSRMatrix
{
public:
    SymCSRMatrix(size_t n_cols, const CSRPattern& pattern);

    Vector gemv(const Vector& x) const override;

    void add_entry(size_t index, double value) override;

    inline bool is_symmetric() const override { return true; }

private:
    const double &get_entry(size_t i, size_t j) const override;
};

} // namespace linalg
} // namespace FastFem

#endif // SYMCSRMATRIX_HPP