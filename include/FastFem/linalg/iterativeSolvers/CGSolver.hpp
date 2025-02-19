#ifndef FASTFEM_CGSOLVER_HPP
#define FASTFEM_CGSOLVER_HPP

#include <FastFem/linalg/iterativeSolvers/IterativeSolver.hpp>

namespace fastfem{
namespace linalg{

class CGSolver : public IterativeSolver
{
public:
    CGSolver(int maxIter = 1000, double tol = 1e-6)
        : IterativeSolver(maxIter, tol) {}

    ~CGSolver() = default;

private:
    Vector p;
    Vector Ap;
    double alpha;
    double beta;
    double r2;
    double b2;

    double initialize(const SparseMatrix& A, const Vector& b, Vector& x) override;
    double iterate(const SparseMatrix& A, const Vector& b, Vector& x) override;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_CGSOLVER_HPP