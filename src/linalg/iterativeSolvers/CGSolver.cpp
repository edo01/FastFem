#include "FastFem/linalg/iterativeSolvers/CGSolver.hpp"

namespace fastfem{
namespace linalg{

double CGSolver::iterate(const SparseMatrix& A, const Vector& b, Vector& x)
{
    Ap = A * p;

    alpha = r2 / p.dot(Ap);

    x.axpy(alpha, p);

    r.axpy(-alpha, Ap);

    double r2_new = r.dot(r);

    beta = r2_new / r2;

    p.axpby(beta, 1.0, r);

    r2 = r2_new;

    return std::sqrt(r2 / b2);
}

} // namespace linalg
} // namespace FastFem
