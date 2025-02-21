#include "FastFem/fe/FESimplexP.hpp"

#include <vector>
#include <stdexcept>


inline double SQUARE(double x) { return x * x; }

inline double compute_den(double xa, double ya, double xb, double yb, double xc, double yc) {
    return SQUARE(xc * (-ya + yb) + xb * (ya - yc) + xa * (-yb + yc));
}

namespace fastfem{
namespace fe{

template<>
void FESimplexP3<2, 2>::
compute_stiffness_loc(const mesh::Simplex<2, 2> &elem, linalg::FullMatrix &matrix) const{

    double xa = elem.get_vertex(0).coords[0];
    double ya = elem.get_vertex(0).coords[1];
    double xb = elem.get_vertex(1).coords[0];
    double yb = elem.get_vertex(1).coords[1];
    double xc = elem.get_vertex(2).coords[0];
    double yc = elem.get_vertex(2).coords[1];

 
    matrix(0,0) = (17 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 40.0;
    matrix(0,1) = (3 * (SQUARE(xb) - 19 * xa * (xb - xc) + 17 * xb * xc - 18 * SQUARE(xc) - 19 * ya * yb + SQUARE(yb) + 19 * ya * yc + 17 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix(0,2) = (3 * (SQUARE(xb) + 8 * xa * (xb - xc) - 10 * xb * xc + 9 * SQUARE(xc) + 8 * ya * yb + SQUARE(yb) - 8 * ya * yc - 10 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix(0,3) = (-7 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix(0,4) = (-3 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 80.0;
    matrix(0,5) = (-3 * (SQUARE(xb) - 2 * xb * xc + SQUARE(xc) + SQUARE(yb - yc))) / 80.0;
    matrix(0,6) = (7 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix(0,7) = (3 * (-8 * xa * xb + 9 * SQUARE(xb) + 8 * xa * xc - 10 * xb * xc + SQUARE(xc) - 8 * ya * yb + 9 * SQUARE(yb) + 8 * ya * yc - 10 * yb * yc + SQUARE(yc))) / 80.0;
    matrix(0,8) = (3 * (-18 * SQUARE(xb) + 19 * xa * (xb - xc) + 17 * xb * xc + SQUARE(xc) + 19 * ya * yb - 18 * SQUARE(yb) - 19 * ya * yc + 17 * yb * yc + SQUARE(yc))) / 80.0;
    matrix(0,9) = 0;
    matrix(1,1) = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix(1,2) = (-27 * (SQUARE(xa) + SQUARE(xb) + 2 * xa * (xb - 2 * xc) - 4 * xb * xc + 4 * SQUARE(xc) + SQUARE(ya) + 2 * ya * yb + SQUARE(yb) - 4 * ya * yc - 4 * yb * yc + 4 * SQUARE(yc))) / 80.0;
    matrix(1,3) = (3 * (SQUARE(xa) + 8 * xa * xb - 10 * xa * xc - 8 * xb * xc + 9 * SQUARE(xc) + SQUARE(ya) + 8 * ya * yb - 10 * ya * yc - 8 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix(1,4) = (-27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 80.0;
    matrix(1,5) = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix(1,6) = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix(1,7) = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix(1,8) = (-27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 16.0;
    matrix(1,9) = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;
    matrix(2,2) = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix(2,3) = (3 * (SQUARE(xa) - 19 * xa * xb + 17 * xa * xc + 19 * xb * xc - 18 * SQUARE(xc) + SQUARE(ya) - 19 * ya * yb + 17 * ya * yc + 19 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix(2,4) = (27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 16.0;
    matrix(2,5) = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix(2,6) = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix(2,7) = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix(2,8) = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix(2,9) = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;
    matrix(3,3) = (17 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 40.0;
    matrix(3,4) = (-3 * (18 * SQUARE(xa) - 19 * xa * xb - 17 * xa * xc + 19 * xb * xc - SQUARE(xc) + 18 * SQUARE(ya) - 19 * ya * yb - 17 * ya * yc + 19 * yb * yc - SQUARE(yc))) / 80.0;
    matrix(3,5) = (3 * (9 * SQUARE(xa) + 8 * xb * xc + SQUARE(xc) - 2 * xa * (4 * xb + 5 * xc) + 9 * SQUARE(ya) - 8 * ya * yb - 10 * ya * yc + 8 * yb * yc + SQUARE(yc))) / 80.0;
    matrix(3,6) = (-7 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix(3,7) = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix(3,8) = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix(3,9) = 0;
    matrix(4,4) = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix(4,5) = (-27 * (4 * SQUARE(xa) + SQUARE(xb) + 2 * xb * xc + SQUARE(xc) - 4 * xa * (xb + xc) + 4 * SQUARE(ya) - 4 * ya * yb + SQUARE(yb) - 4 * ya * yc + 2 * yb * yc + SQUARE(yc))) / 80.0;
    matrix(4,6) = (3 * (9 * SQUARE(xa) + SQUARE(xb) + 8 * xb * xc - 2 * xa * (5 * xb + 4 * xc) + (ya - yb) * (9 * ya - yb - 8 * yc))) / 80.0;
    matrix(4,7) = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix(4,8) = (27 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 80.0;
    matrix(4,9) = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;
    matrix(5,5) = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix(5,6) = (-3 * (xa - xb) * (18 * xa + xb - 19 * xc) - 3 * (ya - yb) * (18 * ya + yb - 19 * yc)) / 80.0;
    matrix(5,7) = (-27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 16.0;
    matrix(5,8) = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix(5,9) = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;
    matrix(6,6) = (17 * (SQUARE(xa - xb) + SQUARE(ya - yb))) / 40.0;
    matrix(6,7) = (3 * (SQUARE(xa) - 18 * SQUARE(xb) + xa * (17 * xb - 19 * xc) + 19 * xb * xc + (ya - yb) * (ya + 18 * yb - 19 * yc))) / 80.0;
    matrix(6,8) = (3 * (SQUARE(xa) + 9 * SQUARE(xb) - 8 * xb * xc + xa * (-10 * xb + 8 * xc) + (ya - yb) * (ya - 9 * yb + 8 * yc))) / 80.0;
    matrix(6,9) = 0;
    matrix(7,7) = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix(7,8) = (-27 * (SQUARE(xa - 2 * xb + xc) + SQUARE(ya - 2 * yb + yc))) / 80.0;
    matrix(7,9) = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;
    matrix(8,8) = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix(8,9) = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;
    matrix(9,9) = (81 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 20.0;

    double den = compute_den(xa, ya, xb, yb, xc, yc);
    double jacobian = elem.volume() * 2;
    for (size_t i = 0; i < 10; ++i){
        for (size_t j = i; j < 10; ++j){
            matrix(i,j) *= jacobian/den;
            if (i != j)
                matrix(j,i) = matrix(i,j);
        }
    }
}

template<>
void FESimplexP3<1, 1>::
compute_stiffness_loc(const mesh::Simplex<1, 1>& /*elem*/, linalg::FullMatrix& /*matrix*/) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

template<>
void FESimplexP3<1, 2>::
compute_stiffness_loc(const mesh::Simplex<1, 2>& /*elem*/, linalg::FullMatrix& /*matrix*/) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

template<>
void FESimplexP3<1, 3>::
compute_stiffness_loc(const mesh::Simplex<1, 3>& /*elem*/, linalg::FullMatrix& /*matrix*/) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

template<>
void FESimplexP3<2, 3>::
compute_stiffness_loc(const mesh::Simplex<2, 3>& /*elem*/, linalg::FullMatrix& /*matrix*/) const{
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

/**
 * @brief Constructor for the 1D case, in 1D space.
 */
template<>
FESimplexP3<1, 1>::
FESimplexP3(const unsigned int n_components) : FESimplexP<1, 1>(n_components) {

    if (n_components != 1) throw std::runtime_error("FESimplexP3 only supports scalar fields");

    // compute the number of degrees of freedom
    this->n_dofs_per_element = 10 * n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 2 * n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     * (0)  [0]---[1]---[2]---[3]  (1)
     *
     */
    mesh::Point<1> p0 = {0};
    mesh::Point<1> p1 = {1.0 / 3.0};
    mesh::Point<1> p2 = {2.0 / 3.0};
    mesh::Point<1> p3 = {1};
    this->reference_simplex = mesh::Simplex<1, 1>({p0, p3});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);

    // fill the tables
    this->vertex_dofs[{0}] = {0};
    this->vertex_dofs[{1}] = {3};
    this->edge_dofs[{0, 1}] = {1, 2};
}

/**
 * @brief Constructor for the 1D case, in 2D space.
 */
template<>
FESimplexP3<1, 2>::
FESimplexP3(const unsigned int n_components) : FESimplexP<1, 2>(n_components) {

    if (n_components != 1) throw std::runtime_error("FESimplexP3 only supports scalar fields");

    // compute the number of degrees of freedom
    this->n_dofs_per_element = 10 * n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 2 * n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     * (0)  [0]---[1]---[2]---[3]  (1)
     *
     */
    mesh::Point<2> p0 = {0, 0};
    mesh::Point<2> p1 = {1.0 / 3.0, 0};
    mesh::Point<2> p2 = {2.0 / 3.0, 0};
    mesh::Point<2> p3 = {1, 0};
    this->reference_simplex = mesh::Simplex<1, 2>({p0, p3});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);

    // fill the tables
    this->vertex_dofs[{0}] = {0};
    this->vertex_dofs[{1}] = {3};
    this->edge_dofs[{0, 1}] = {1, 2};
}

/**
 * @brief Constructor for the 1D case, in 3D space.
 */
template<>
FESimplexP3<1, 3>::
FESimplexP3(const unsigned int n_components) : FESimplexP<1, 3>(n_components) {

    if (n_components != 1) throw std::runtime_error("FESimplexP3 only supports scalar fields");

    // compute the number of degrees of freedom
    this->n_dofs_per_element = 10 * n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 2 * n_components;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     * (0)  [0]---[1]---[2]---[3]  (1)
     *
     */
    mesh::Point<3> p0 = {0, 0, 0};
    mesh::Point<3> p1 = {1.0 / 3.0, 0, 0};
    mesh::Point<3> p2 = {2.0 / 3.0, 0, 0};
    mesh::Point<3> p3 = {1, 0, 0};
    this->reference_simplex = mesh::Simplex<1, 3>({p0, p3});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);

    // fill the tables
    this->vertex_dofs[{0}] = {0};
    this->vertex_dofs[{1}] = {3};
    this->edge_dofs[{0, 1}] = {1, 2};
}

/**
 * @brief Constructor for the 2D case, in 2D space.
 */
template<>
FESimplexP3<2, 2>::
FESimplexP3(const unsigned int n_components) : FESimplexP<2, 2>(n_components) {

    if (n_components != 1) throw std::runtime_error("FESimplexP3 only supports scalar fields");

    // compute the number of degrees of freedom
    this->n_dofs_per_element = 10 * n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 2 * n_components;
    this->n_dofs_per_face = n_components;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * (2) [6]
     *     |  \
     *     |   \
     *     |    \
     *     [7]   [5]
     *     |      \
     *     |       \
     *     |        \
     *     [8] [9]  [4]
     *     |          \
     *     |           \
     *     |            \
     *     |             \
     *     |              \
     * (0) [0]--[1]--[2]--[3] (1)
     * 
     */
    mesh::Point<2> p0 = {0, 0};
    mesh::Point<2> p1 = {1.0 / 3.0, 0};
    mesh::Point<2> p2 = {2.0 / 3.0, 0};
    mesh::Point<2> p3 = {1, 0};
    mesh::Point<2> p4 = {2.0 / 3.0, 1.0 / 3.0};
    mesh::Point<2> p5 = {1.0 / 3.0, 2.0 / 3.0};
    mesh::Point<2> p6 = {0, 1.0};
    mesh::Point<2> p7 = {0, 2.0 / 3.0};
    mesh::Point<2> p8 = {0, 1.0 / 3.0};
    mesh::Point<2> p9 = {1.0 / 3.0, 1.0 / 3.0};
    this->reference_simplex = mesh::Simplex<2>({p0, p3, p6});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);
    this->dofs.push_back(p4);
    this->dofs.push_back(p5);
    this->dofs.push_back(p6);
    this->dofs.push_back(p7);
    this->dofs.push_back(p8);
    this->dofs.push_back(p9);

    // fill the tables
    this->vertex_dofs[{0}] = {0};
    this->vertex_dofs[{1}] = {3};
    this->vertex_dofs[{2}] = {6};
    this->edge_dofs[{0, 1}] = {1, 2};
    this->edge_dofs[{1, 2}] = {4, 5};
    this->edge_dofs[{0, 2}] = {7, 8};
    this->face_dofs[{0, 1, 2}] = {9};

}
    

/**
 * @brief Constructor for the 2D case, in 3D space.
 */
template<>
FESimplexP3<2, 3>::
FESimplexP3(const unsigned int n_components) : FESimplexP<2, 3>(n_components) {

    if (n_components != 1) throw std::runtime_error("FESimplexP3 only supports scalar fields");

    // compute the number of degrees of freedom
    this->n_dofs_per_element = 10 * n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 2 * n_components;
    this->n_dofs_per_face = n_components;
    this->n_dofs_per_cell = 0;

    // set the reference simplex
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * (2) [6]
     *     |  \
     *     |   \
     *     |    \
     *     [7]   [5]
     *     |      \
     *     |       \
     *     |        \
     *     [8] [9]  [4]
     *     |          \
     *     |           \
     *     |            \
     *     |             \
     *     |              \
     * (0) [0]--[1]--[2]--[3] (1)
     * 
     */
    mesh::Point<3> p0 = {0, 0, 0};
    mesh::Point<3> p1 = {1.0 / 3.0, 0, 0};
    mesh::Point<3> p2 = {2.0 / 3.0, 0, 0};
    mesh::Point<3> p3 = {1, 0, 0};
    mesh::Point<3> p4 = {2.0 / 3.0, 1.0 / 3.0, 0};
    mesh::Point<3> p5 = {1.0 / 3.0, 2.0 / 3.0, 0};
    mesh::Point<3> p6 = {0, 1, 0};
    mesh::Point<3> p7 = {0, 2.0 / 3.0, 0};
    mesh::Point<3> p8 = {0, 1.0 / 3.0, 0};
    mesh::Point<3> p9 = {1.0 / 3.0, 1.0 / 3.0, 0};
    this->reference_simplex = mesh::Simplex<2, 3>({p0, p3, p6});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);
    this->dofs.push_back(p4);
    this->dofs.push_back(p5);
    this->dofs.push_back(p6);
    this->dofs.push_back(p7);
    this->dofs.push_back(p8);
    this->dofs.push_back(p9);

    // fill the tables
    this->vertex_dofs[{0}] = {0};
    this->vertex_dofs[{1}] = {3};
    this->vertex_dofs[{2}] = {6};
    this->edge_dofs[{0, 1}] = {1, 2};
    this->edge_dofs[{1, 2}] = {4, 5};
    this->edge_dofs[{0, 2}] = {7, 8};
    this->face_dofs[{0, 1, 2}] = {9};
}


}// namespace fe
}// namespace fastfem