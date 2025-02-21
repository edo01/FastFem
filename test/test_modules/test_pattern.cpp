#include <iostream>
#include <vector>
#include <set>
#include <iomanip>
#include <FastFem/dof/DofHandler.hpp>
#include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>
#include <FastFem/linalg/sparseMatrices/SkylineMatrix.hpp>
#include <FastFem/fe/FESimplexP.hpp>
#include <FastFem/mesh/MeshMaker.hpp>
#include <FastFem/mesh/Mesh.hpp>

using namespace fastfem::mesh;
using namespace fastfem::fe;
using namespace fastfem::dof;

/**
 * @brief This test visually compares the skyline pattern with the CSR pattern 
*/
int main() {

    SquareMaker square(1);
    Mesh<2, 2> mesh = square.make_mesh();

    FESimplexP3<2, 2> fe(1);
    DoFHandler<2, 2> dof_handler(mesh);

    dof_handler.distribute_dofs(std::make_shared<FESimplexP3<2, 2>>(fe));

    fastfem::linalg::CSRPattern csr_pattern = fastfem::linalg::CSRPattern::create_from_dof_handler(dof_handler);
    fastfem::linalg::CSRPattern symmetric_pattern = fastfem::linalg::CSRPattern::create_symmetric_from_dof_handler(dof_handler);
    fastfem::linalg::SkylinePattern skyline_pattern = fastfem::linalg::SkylinePattern::create_from_dof_handler(dof_handler);
  
    
    size_t n_dofs = csr_pattern.row_ptr.size() - 1;

    fastfem::linalg::CSRMatrix csr_matrix(n_dofs, csr_pattern);
    fastfem::linalg::CSRMatrix symmetric_csr_matrix(n_dofs, symmetric_pattern);
    fastfem::linalg::SkylineMatrix skyline_matrix(n_dofs, skyline_pattern);

    std::cout << "\nCSR Matrix (Standard Pattern):\n";
    csr_matrix.print_pattern();

    std::cout << "\n";
    std::cout << "\nCSR Matrix (Symmetric Pattern):\n";
    symmetric_csr_matrix.print_pattern();

    std::cout << "\n";
    std::cout << "\nSkyline Matrix:\n";
    skyline_matrix.print_pattern(false);

}
