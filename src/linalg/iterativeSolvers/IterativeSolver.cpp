#include "FastFem/linalg/iterativeSolvers/IterativeSolver.hpp"

namespace fastfem{
namespace linalg{

Vector IterativeSolver::solve(const SparseMatrix& A, const Vector& b)
{
    if(A.get_n_cols() != b.size()){
        throw std::invalid_argument("IterativeSolver::solve(): incompatible dimensions");
    }

    if(A.get_n_rows() != A.get_n_cols()){
        throw std::invalid_argument("IterativeSolver::solve(): matrix must be square");
    }

    Vector x(b.size());
    r = b - A * x;

    r2 = r.dot(r);
    b2 = b.dot(b);

    unsigned iter = 0;
    double error;
    do {
        error = iterate(A, b, x);
        ++iter;
    } while(iter < maxIterations && error > tolerance);

    if(iter == maxIterations){
        throw std::runtime_error("IterativeSolver::solve(): maximum number of iterations reached");
    }

    return x;
}

} // namespace linalg
} // namespace FastFem