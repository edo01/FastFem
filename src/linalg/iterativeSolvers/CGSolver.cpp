#include "FastFem/linalg/iterativeSolvers/CGSolver.hpp"

namespace fastfem{
namespace linalg{

double CGSolver::initialize(const SparseMatrix& A, const Vector& b, Vector& x)
{

    b2 = b.dot(b);
    if(b2 < tolerance){
        return 0.0;
    }

    // r_0 = b - A * x_0
    r = b - A * x;

    // p_0 = r_0
    p = r;

    r2 = r.dot(r);

    return std::sqrt(r2 / b2);
}

double CGSolver::iterate(const SparseMatrix& A, const Vector& b, Vector& x)
{
    Ap = A * p;

    double pAp = p.dot(Ap);

    alpha = r2 / pAp;

    Vector::axpy(alpha, p, x);

    Vector::axpy(-alpha, Ap, r);

    double r2_new = r.dot(r);

    beta = r2_new / r2;

    Vector::axpby(1.0, r, beta, p);

    r2 = r2_new;

    return std::sqrt(r2 / b2);
}

} // namespace linalg
} // namespace fastfem
