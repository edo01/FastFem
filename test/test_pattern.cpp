// #include <iostream>
// #include <vector>
// #include <set>
// #include <FastFem/dof/DofHandler.hpp>  // Include the DofHandler header
// #include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>  // Include CSRPattern

// // Print function for DofHandler
// template <unsigned int dim, unsigned int spacedim>
// void print_dof_handler(const fastfem::dof::DofHandler<dim, spacedim>& dof_handler) {
//     std::cout << "DofHandler with " << dof_handler.get_n_elements() << " elements:\n";
//     for (unsigned int i = 0; i < dof_handler.get_n_elements(); ++i) {
//         const auto& dofs = dof_handler.get_element_dofs(i);
//         std::cout << "Element " << i << ": ";
//         for (unsigned int dof : dofs) {
//             std::cout << dof << " ";
//         }
//         std::cout << "\n";
//     }
// }

// // Print function for CSRPattern
// void print_csr_pattern(const fastfem::linalg::CSRPattern& csr_pattern) {
//     std::cout << "CSR Pattern:\n";
    
//     std::cout << "Row Pointers: ";
//     for (size_t rp : csr_pattern.row_ptr) {
//         std::cout << rp << " ";
//     }
//     std::cout << "\n";

//     std::cout << "Column Indices: ";
//     for (size_t ci : csr_pattern.col_indices) {
//         std::cout << ci << " ";
//     }
//     std::cout << "\n";
// }

// int main() {
//     const unsigned int num_elements = 15; // Define number of elements

//     // Create a DofHandler with random DOFs per element (dim=2)
//     fastfem::dof::DofHandler<2,2> dof_handler(num_elements);

//     // Print the DofHandler structure
//     print_dof_handler(dof_handler);

//     // Generate and print CSR Pattern (Full)
//     std::cout << "\nCSR Pattern (Full):\n";
//     // Specify template parameters when calling the function
//     fastfem::linalg::CSRPattern csr_pattern = fastfem::linalg::CSRPattern::create_from_dof_handler(dof_handler);
//     print_csr_pattern(csr_pattern);

//     // // Generate and print Symmetric CSR Pattern (Lower Triangular)
//     // std::cout << "\nCSR Pattern (Symmetric, Lower Triangular):\n";
//     // // Specify template parameters when calling the function
//     // fastfem::linalg::CSRPattern csr_symmetric = fastfem::linalg::CSRPattern::create_symmetric_from_dof_handler<2, 2>(dof_handler);
//     // print_csr_pattern(csr_symmetric);

//     return 0;
// }


// #include <iostream>
// #include <vector>
// #include <set>
// #include <cassert>
// #include <FastFem/dof/DofHandler.hpp>
// #include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>

// // Function to validate the CSR pattern
// bool test_csr_pattern_correctness(const fastfem::dof::DofHandler<2, 2>& dof_handler, const fastfem::linalg::CSRPattern& csr_pattern) {
//     const unsigned int num_elements = dof_handler.get_n_elements();
//     const size_t num_dofs = csr_pattern.row_ptr.size() - 1;  // Number of DOFs (rows in CSR)

//     // 1. Check row pointer consistency (non-decreasing order)
//     for (size_t i = 1; i < csr_pattern.row_ptr.size(); ++i) {
//         if (csr_pattern.row_ptr[i] < csr_pattern.row_ptr[i - 1]) {
//             std::cerr << "Error: Row pointers are not non-decreasing at index " << i << "\n";
//             return false;
//         }
//     }

//     // 2. Check that column indices are within valid range
//     for (size_t ci : csr_pattern.col_indices) {
//         if (ci >= num_dofs) {
//             std::cerr << "Error: Column index " << ci << " out of bounds.\n";
//             return false;
//         }
//     }

//     // 3. Verify the number of nonzeros per row
//     std::vector<std::set<unsigned int>> expected_rows(num_dofs);  // Expected structure
//     for (unsigned int elem = 0; elem < num_elements; ++elem) {
//         const auto& dofs = dof_handler.get_element_dofs(elem);
//         for (unsigned int i : dofs) {
//             for (unsigned int j : dofs) {
//                 expected_rows[i].insert(j);
//             }
//         }
//     }

//     // Convert expected structure to CSR format and compare
//     std::vector<size_t> computed_row_ptr(num_dofs + 1, 0);
//     std::vector<size_t> computed_col_indices;

//     for (size_t i = 0; i < num_dofs; ++i) {
//         computed_row_ptr[i + 1] = computed_row_ptr[i] + expected_rows[i].size();
//         computed_col_indices.insert(computed_col_indices.end(), expected_rows[i].begin(), expected_rows[i].end());
//     }

//     // Compare computed values with CSR pattern
//     if (csr_pattern.row_ptr != computed_row_ptr) {
//         std::cerr << "Error: Row pointers do not match expected values.\n";
//         return false;
//     }

//     if (csr_pattern.col_indices != computed_col_indices) {
//         std::cerr << "Error: Column indices do not match expected values.\n";
//         return false;
//     }

//     std::cout << "CSR pattern is correct!\n";
//     return true;
// }

// int main() {
//     const unsigned int num_elements = 15;
//     fastfem::dof::DofHandler<2,2> dof_handler(num_elements);
//     fastfem::linalg::CSRPattern csr_pattern = fastfem::linalg::CSRPattern::create_from_dof_handler(dof_handler);


//     // Validate CSR pattern
//     if (!test_csr_pattern_correctness(dof_handler, csr_pattern)) {
//         return EXIT_FAILURE;
//     }

//     return EXIT_SUCCESS;
// }


// #include <iostream>
// #include <vector>
// #include <iomanip>
// #include <FastFem/dof/DofHandler.hpp>
// #include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>

// // Function to print a matrix in dense format
// void print_dense_matrix_from_dofhandler(const fastfem::dof::DofHandler<2, 2>& dof_handler, size_t n_dofs) {
//     std::vector<std::vector<double>> dense_matrix(n_dofs, std::vector<double>(n_dofs, 0.0));

//     const unsigned int num_elements = dof_handler.get_n_elements();

//     for (unsigned int elem = 0; elem < num_elements; ++elem) {
//         const auto& dofs = dof_handler.get_element_dofs(elem);
//         for (unsigned int i : dofs) {
//             for (unsigned int j : dofs) {
//                 dense_matrix[i][j] += 1.0; // Example value, can be replaced with real computations
//             }
//         }
//     }

//     std::cout << "Dense Matrix from DofHandler:\n";
//     for (size_t i = 0; i < n_dofs; ++i) {
//         for (size_t j = 0; j < n_dofs; ++j) {
//             std::cout << std::setw(3) << dense_matrix[i][j] << " ";
//         }
//         std::cout << "\n";
//     }
// }

// int main() {
//     const unsigned int num_elements = 10;  // Keep it small for testing (< 25 DOFs)
//     fastfem::dof::DofHandler<2,2> dof_handler(num_elements);

//     fastfem::linalg::CSRPattern csr_pattern = fastfem::linalg::CSRPattern::create_from_dof_handler(dof_handler);
//     size_t n_dofs = csr_pattern.row_ptr.size() - 1;  // Number of DOFs

//     // Ensure the matrix size is within the limit for testing
//     if (n_dofs >= 25) {
//         std::cerr << "Error: Too many DOFs for testing. Adjust the number of elements.\n";
//         return EXIT_FAILURE;
//     }

//     // Create a CSR matrix
//     fastfem::linalg::CSRMatrix csr_matrix(n_dofs, csr_pattern);

//     // Fill the matrix with some test values
//     for (size_t i = 0; i < csr_matrix.nnz(); ++i) {
//         csr_matrix.add_entry(i, 1.0);  // Assigning a test value of 1.0
//     }

//     // Print CSR matrix pattern
//     std::cout << "\nCSR Matrix (stored in CSR format):\n";
//     csr_matrix.print_pattern();
//     std::cout << "\n";
//     std::cout << "\n";
//     csr_matrix.print();

//     // Print matrix as dense representation deduced from the DofHandler
//     print_dense_matrix_from_dofhandler(dof_handler, n_dofs);

//     return EXIT_SUCCESS;
// }



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


// // Function to assemble the CSRMatrix from the DofHandler
// void assemble_csr_matrix(fastfem::linalg::CSRMatrix& matrix, const fastfem::dof::DofHandler<2, 2>& dof_handler) {
//     const auto& csr_pattern = *matrix.base_pattern;
//     const unsigned int num_elements = dof_handler.get_n_elements();

//     for (unsigned int elem = 0; elem < num_elements; ++elem) {
//         const auto& dofs = dof_handler.get_element_dofs(elem);
        
//         for (size_t i = 0; i < dofs.size(); ++i) {
//             for (size_t j = 0; j < dofs.size(); ++j) {
//                 size_t row = dofs[i];
//                 size_t col = dofs[j];

//                 double value = static_cast<double>(row + col) / 2.0;

//                 size_t start = csr_pattern.row_ptr[row];
//                 size_t end = csr_pattern.row_ptr[row + 1];

//                 for (size_t k = start; k < end; ++k) {
//                     if (csr_pattern.col_indices[k] == col) {
//                         matrix.add_entry(k, value);
//                         break;
//                     }
//                 }
//             }
//         }
//     }
// }

// // Function to assemble a symmetric CSRMatrix
// void assemble_symmetric_csr_matrix(fastfem::linalg::CSRMatrix& matrix, const fastfem::dof::DofHandler<2, 2>& dof_handler) {
//     const auto& csr_pattern = *matrix.base_pattern;
//     const unsigned int num_elements = dof_handler.get_n_elements();

//     for (unsigned int elem = 0; elem < num_elements; ++elem) {
//         const auto& dofs = dof_handler.get_element_dofs(elem);
        
//         for (size_t i = 0; i < dofs.size(); ++i) {
//             for (size_t j = i; j < dofs.size(); ++j) {  // Ensure symmetric storage
//                 size_t row = std::min(dofs[i], dofs[j]);
//                 size_t col = std::max(dofs[i], dofs[j]);

//                 double value = static_cast<double>(row + col) / 2.0;

//                 size_t start = csr_pattern.row_ptr[row];
//                 size_t end = csr_pattern.row_ptr[row + 1];

//                 for (size_t k = start; k < end; ++k) {
//                     if (csr_pattern.col_indices[k] == col) {
//                         matrix.add_entry(k, value);
//                         break;
//                     }
//                 }
//             }
//         }
//     }
// }

// Function to print a dense matrix from the DofHandler for validation
void print_dense_matrix_from_dofhandler(const fastfem::dof::DoFHandler<2, 2>& dof_handler, size_t n_dofs) {
    std::vector<std::vector<double>> dense_matrix(n_dofs, std::vector<double>(n_dofs, 0.0));

    const unsigned int num_elements = dof_handler.get_n_elements();

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        const auto& elem = *it;

        const auto& dofs = dof_handler.get_unordered_dofs_on_element(elem);
        
        for (size_t i = 0; i < dofs.size(); ++i) {
            for (size_t j = 0; j < dofs.size(); ++j) {
                size_t row = dofs[i];
                size_t col = dofs[j];

                dense_matrix[row][col] += static_cast<double>(row + col) / 2.0;
            }
        }
    }

    std::cout << "Dense Matrix from DofHandler:\n";
    for (size_t i = 0; i < n_dofs; ++i) {
        for (size_t j = 0; j < n_dofs; ++j) {
            std::cout << std::setw(6) << std::fixed << std::setprecision(1) << dense_matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

using namespace fastfem::mesh;
using namespace fastfem::fe;
using namespace fastfem::dof;

int main() {
    const unsigned int num_elements = 4;

    SquareMaker square(1);
    Mesh<2, 2> mesh = square.make_mesh();

    FESimplexP3<2, 2> fe(1);
    DoFHandler<2, 2> dof_handler(mesh, std::make_unique<FESimplexP3<2, 2>>(fe));

    dof_handler.distribute_dofs();

    //print dof handler
/*     std::cout << "DofHandler with " << dof_handler.get_n_elements() << " elements:\n";
    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        const auto& elem = *it;

        const auto& dofs = dof_handler.get_ordered_dofs_on_element(elem);

        std::cout << "Element " << elem.get_vertex(0) << ", " << elem.get_vertex(1) << ", " << elem.get_vertex(2) << ": "; 
        for (unsigned int dof : dofs) {
            std::cout << dof << " ";
        }
        std::cout << "\n";
    } */

    fastfem::linalg::CSRPattern csr_pattern = fastfem::linalg::CSRPattern::create_from_dof_handler(dof_handler);
    fastfem::linalg::CSRPattern symmetric_pattern = fastfem::linalg::CSRPattern::create_symmetric_from_dof_handler(dof_handler);
    fastfem::linalg::SkylinePattern skyline_pattern = fastfem::linalg::SkylinePattern::create_from_dof_handler(dof_handler);

    //print csr pattern
    // std::cout << "CSR Pattern:\n";
    // std::cout << "Row Pointers: ";
    // for (size_t rp : csr_pattern.row_ptr) {
    //     std::cout << rp << " ";
    // }
    // std::cout << "\n";
    // std::cout << "Column Indices: ";
    // for (size_t ci : csr_pattern.col_indices) {
    //     std::cout << ci << " ";
    // }
    // std::cout << "\n";    
    
    size_t n_dofs = csr_pattern.row_ptr.size() - 1;
    
/*     if (n_dofs >= 10) {
        std::cerr << "Error: Too many DOFs for testing. Adjust the number of elements.\n";
        return EXIT_FAILURE;
    } */

    fastfem::linalg::CSRMatrix csr_matrix(n_dofs, csr_pattern);
    fastfem::linalg::CSRMatrix symmetric_csr_matrix(n_dofs, symmetric_pattern);
    fastfem::linalg::SkylineMatrix skyline_matrix(n_dofs, skyline_pattern);

    //assemble_csr_matrix(csr_matrix, dof_handler);
    //assemble_symmetric_csr_matrix(symmetric_csr_matrix, dof_handler);

    std::cout << "\nCSR Matrix (Standard Pattern):\n";
    csr_matrix.print_pattern();

    std::cout << "\n";
    std::cout << "\nCSR Matrix (Symmetric Pattern):\n";
    symmetric_csr_matrix.print_pattern();

    std::cout << "\n";
    std::cout << "\nSkyline Matrix:\n";
    skyline_matrix.print_pattern();


    // std::cout << "\n";
    // csr_matrix.print();

    // std::cout<<"Add A[0,1] = 1.0\n";
    // csr_matrix.insert_entry(0, 1, 1.0);
    // csr_matrix.print();

    // std::cout << "\n";
    // std::cout << "Add A[4,3] = 2.75\n";
    // csr_matrix.insert_entry(4, 3, 2.75);
    // csr_matrix.print();

    // std::cout << "\n";
    // std::cout << "Add A[5,6] = 3.5\n";
    // csr_matrix.insert_entry(5, 6, 3.5);
    // csr_matrix.print();

    // std::cout << "\n";

    // std::cout << "\nCSR Matrix (Symmetric Pattern):\n";
    // //symmetric_csr_matrix.print_pattern();

    // print_dense_matrix_from_dofhandler(dof_handler, n_dofs);

    // return EXIT_SUCCESS;
}
