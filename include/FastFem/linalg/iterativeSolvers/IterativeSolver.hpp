#ifndef FASTFEM_ITERATIVESOLVER_HPP
#define FASTFEM_ITERATIVESOLVER_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>

namespace fastfem {
namespace linalg {

class IterativeSolver {
protected:
    unsigned int maxIterations;
    double tolerance;

    Vector r;
    unsigned int lastStep;
    double error;

public:
    IterativeSolver(unsigned int maxIter = 1000, double tol = 1e-6)
        : maxIterations(maxIter), tolerance(tol) {}

    virtual ~IterativeSolver() = default;

    virtual Vector solve(const SparseMatrix& A, const Vector& b);

    inline double get_error() const { return error; }
    inline unsigned int get_last_step() const { return lastStep; }

private:
    virtual double initialize(const SparseMatrix& A, const Vector& b, Vector& x) = 0;
    virtual double iterate(const SparseMatrix& A, const Vector& b, Vector& x) = 0;

};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_ITERATIVESOLVER_HPP