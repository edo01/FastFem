#include "FastFem/fe/FESimplexP.hpp"

#include <vector>
#include <stdexcept>

namespace fastfem{
namespace fe{

inline double SQUARE(double x) { return x * x; }

inline double compute_den(double xa, double ya, double xb, double yb, double xc, double yc) {
    return SQUARE(xc * (-ya + yb) + xb * (ya - yc) + xa * (-yb + yc));
}

template<>
void FESimplexP2<2, 2>::
compute_stiffness_loc(const mesh::Simplex<2, 2> &elem, linalg::FullMatrix& matrix) const{

    double xa = elem.get_vertex(0).coords[0];
    double ya = elem.get_vertex(0).coords[1];
    double xb = elem.get_vertex(1).coords[0];
    double yb = elem.get_vertex(1).coords[1];
    double xc = elem.get_vertex(2).coords[0];
    double yc = elem.get_vertex(2).coords[1];

    matrix(0,0) = (SQUARE(xb - xc) + SQUARE(yb - yc)) / 2.;
    matrix(0,1) = (-2 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix(0,2) = ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc)) / 6.;
    matrix(0,3) = 0;
    matrix(0,4) = (-1 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 6.;
    matrix(0,5) = (2 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;
    matrix(1,1) = (4 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - (ya + yb) * yc + SQUARE(yc))) / 3.;
    matrix(1,2) = (-2 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix(1,3) = (4 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;
    matrix(1,4) = 0;
    matrix(1,5) = (-4 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix(2,2) = (SQUARE(xa - xc) + SQUARE(ya - yc)) / 2.;
    matrix(2,3) = (-2 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix(2,4) = (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc)) / 6.;
    matrix(2,5) = 0;
    matrix(3,3) = (4 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - (ya + yb) * yc + SQUARE(yc))) / 3.;
    matrix(3,4) = (-2 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix(3,5) = (-4 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix(4,4) = (SQUARE(xa - xb) + SQUARE(ya - yb)) / 2.;
    matrix(4,5) = (2 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;
    matrix(5,5) = (4 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - (ya + yb) * yc + SQUARE(yc))) / 3.;

    // multiply by jacobian
    double den = compute_den(xa, ya, xb, yb, xc, yc);
    double jacobian = elem.volume() * 2;

    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = i; j < 6; ++j)
        {
            matrix(i,j) *= jacobian/den;
            
            if (i != j) // fill the lower part of the matrix
                matrix(j,i) = matrix(i,j);
        }
    }   
}

template<>
void FESimplexP2<1, 1>::
compute_stiffness_loc(const mesh::Simplex<1, 1> &, linalg::FullMatrix&) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

template<>
void FESimplexP2<1, 2>::
compute_stiffness_loc(const mesh::Simplex<1, 2> &, linalg::FullMatrix &) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

template<>
void FESimplexP2<1, 3>::
compute_stiffness_loc(const mesh::Simplex<1, 3> &, linalg::FullMatrix&) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

template<>
void FESimplexP2<2, 3>::
compute_stiffness_loc(const mesh::Simplex<2, 3> &, linalg::FullMatrix&) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}
/**
 * @brief Constructor for the 1D case, in any space dimension.
 */
template<>
FESimplexP2<1, 1>::
FESimplexP2(const unsigned int n_components) :  FESimplexP<1, 1>(n_components) {

    if(n_components != 1){
        throw std::invalid_argument("FESimplexP2: only scalar functions are supported");
    }

    using local_vertex_index = fastfem::types::local_vertex_index;
    using local_edge_index = fastfem::types::local_edge_index;

    this->n_dofs_per_element = 6*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     *  
     * (0)  [0]--[1]--[2]  (1)
     * 
     */

    // fill the vertex dofs (no dof on faces or cells)
    vertex_dofs[(local_vertex_index){0}] = {0};
    vertex_dofs[(local_vertex_index){1}] = {2};
    edge_dofs[(local_edge_index){0,1}] = {1};

    mesh::Point<1> p0 = {0};
    mesh::Point<1> p1 = {0.5};
    mesh::Point<1> p2 = {1};

    this->reference_simplex = mesh::Simplex<1, 1>({p0, p2});
    
    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
}

/**
 * @brief Constructor for the 1D case, in 2D space.
 */
template<>
FESimplexP2<1, 2>::
FESimplexP2(const unsigned int n_components) :  FESimplexP<1, 2>(n_components) {

    if(n_components != 1){
        throw std::invalid_argument("FESimplexP2: only scalar functions are supported");
    }

    using local_vertex_index = fastfem::types::local_vertex_index;
    using local_edge_index = fastfem::types::local_edge_index;

    this->n_dofs_per_element = 6*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     *  
     * (0)  [0]--[1]--[2]  (1)
     * 
     */

    // fill the vertex dofs (no dof on faces or cells)
    vertex_dofs[(local_vertex_index){0}] = {0};
    vertex_dofs[(local_vertex_index){1}] = {2};
    edge_dofs[(local_edge_index){0,1}] = {1};

    mesh::Point<2> p0 = {0, 0};
    mesh::Point<2> p1 = {0.5, 0};
    mesh::Point<2> p2 = {1, 0};

    this->reference_simplex = mesh::Simplex<1, 2>({p0, p2});
    
    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
}

/**
 * @brief Constructor for the 1D case, in 3D space.
 */
template<>
FESimplexP2<1, 3>::
FESimplexP2(const unsigned int n_components) :  FESimplexP<1, 3>(n_components) {

    if(n_components != 1){
        throw std::invalid_argument("FESimplexP2: only scalar functions are supported");
    }

    using local_vertex_index = fastfem::types::local_vertex_index;
    using local_edge_index = fastfem::types::local_edge_index;

    this->n_dofs_per_element = 6*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     *  
     * (0)  [0]--[1]--[2]  (1)
     * 
     */

    // fill the vertex dofs (no dof on faces or cells)
    vertex_dofs[(local_vertex_index){0}] = {0};
    vertex_dofs[(local_vertex_index){1}] = {2};
    edge_dofs[(local_edge_index){0,1}] = {1};

    mesh::Point<3> p0 = {0, 0, 0};
    mesh::Point<3> p1 = {0.5, 0, 0};
    mesh::Point<3> p2 = {1, 0, 0};
    
    this->reference_simplex = mesh::Simplex<1, 3>({p0, p2});
    
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
}  

/**
 * @brief Constructor for the 2D case, in 2D space.
 */
template<>
FESimplexP2<2, 2>::
FESimplexP2(const unsigned int n_components) :  FESimplexP<2, 2>(n_components) {

    if(n_components != 1){
        throw std::invalid_argument("FESimplexP2: only scalar functions are supported");
    }

    using local_vertex_index = fastfem::types::local_vertex_index;
    using local_edge_index = fastfem::types::local_edge_index;

    this->n_dofs_per_element = 6*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     *  [4] (2)
     *  |  \
     *  |   \
     *  |    \
     *  [5]  [3]
     *  |       \
     *  |        \
     *  |         \
     *  [0]--[1]--[2]
     * (0)          (1)
     */

    // fill the vertex dofs (no dof on faces or cells)
    vertex_dofs[(local_vertex_index){0}] = {0};
    vertex_dofs[(local_vertex_index){1}] = {2};
    vertex_dofs[(local_vertex_index){2}] = {4};
    edge_dofs[(local_edge_index){0,2}] = {5};
    edge_dofs[(local_edge_index){0,1}] = {1};
    edge_dofs[(local_edge_index){1,2}] = {3};

    mesh::Point<2> p0 = {0, 0};
    mesh::Point<2> p1 = {0.5, 0};
    mesh::Point<2> p2 = {1, 0};
    mesh::Point<2> p3 = {0.5, 0.5};
    mesh::Point<2> p4 = {0, 1};
    mesh::Point<2> p5 = {0, 0.5};

    this->reference_simplex = mesh::Simplex<2, 2>({p0, p2, p4});
    
    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);
    this->dofs.push_back(p4);
    this->dofs.push_back(p5);
}

/**
 * @brief Constructor for the 2D case, in 3D space.
 */
template<>
FESimplexP2<2, 3>::
FESimplexP2(const unsigned int n_components) :  FESimplexP<2, 3>(n_components) {

    if(n_components != 1){
        throw std::invalid_argument("FESimplexP2: only scalar functions are supported");
    }

    using local_vertex_index = fastfem::types::local_vertex_index;
    using local_edge_index = fastfem::types::local_edge_index;

    this->n_dofs_per_element = 6*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     *  [4] (2)
     *  |  \
     *  |   \
     *  |    \
     *  [5]  [3]
     *  |       \
     *  |        \
     *  |         \
     *  [0]--[1]--[2]
     * (0)          (1)
     */

    // fill the vertex dofs (no dof on faces or cells)
    vertex_dofs[(local_vertex_index){0}] = {0};
    vertex_dofs[(local_vertex_index){1}] = {2};
    vertex_dofs[(local_vertex_index){2}] = {4};
    edge_dofs[(local_edge_index){0,2}] = {5};
    edge_dofs[(local_edge_index){0,1}] = {1};
    edge_dofs[(local_edge_index){1,2}] = {3};

    mesh::Point<3> p0 = {0, 0, 0};
    mesh::Point<3> p1 = {0.5, 0, 0};
    mesh::Point<3> p2 = {1, 0, 0};
    mesh::Point<3> p3 = {0.5, 0.5, 0};
    mesh::Point<3> p4 = {0, 1, 0};
    mesh::Point<3> p5 = {0, 0.5, 0};

    this->reference_simplex = mesh::Simplex<2, 3>({p0, p2, p4});
    
    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);
    this->dofs.push_back(p4);
    this->dofs.push_back(p5);
}

} // namespace fe
} // namespace fastfem