#ifndef NAVIER_STOKES_SOLVER_HPP
#define NAVIER_STOKES_SOLVER_HPP

#include <vector>
#include <functional>
#include "mesh/mesh.hpp"
#include "linalg/sparse_matrix.hpp"
#include "linalg/system.hpp"

namespace fem::navier_stokes
{

class Navier_stokes_solver
{
    public: 
        Navier_stokes_solver(const Mesh& m, double T, double dt);

        void build_fem_matrices();
        void solve();
    
    private:
        const SparseMatrix S;
        const SparseMatrix M;
        const Mesh&         mesh;
        double dt;
        double T;
        int N_DOFs;
        std::vector<double> omega;
        std::vector<double> Momega;
        std::vector<double> psi;
        std::vector<double> M_col_sum;
        double volume; // volume of the domain


        // conjugate gradient buffers

        void step();
        auto omega_initializer = [](double x, double y, double z) -> double {
            return x*x+y*y+z*z;
        }

}

} // fem::navier_stokes

#endif // NAVIER_STOKES_SOLVER_GPP