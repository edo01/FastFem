#ifndef FASTFEM_CGSOLVER_HPP
#define FASTFEM_CGSOLVER_HPP

#include <FastFem/linalg/iterativeSolvers/IterativeSolver.hpp>

namespace fastfem{
namespace linalg{

class CGSolver : public IterativeSolver
{
    using IterativeSolver::IterativeSolver;

private:
    Vector p;
    Vector Ap;
    double alpha;
    double beta;

    double iterate(const SparseMatrix& A, const Vector& b, Vector& x) override;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_CGSOLVER_HPP