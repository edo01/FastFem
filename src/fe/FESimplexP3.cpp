#include "FastFem/fe/FESimplexP.hpp"

#include <vector>
#include <stdexcept>

namespace fastfem{
namespace fe{

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