#include "FastFem/fe/FESimplexP.hpp"

#include <vector>
#include <stdexcept>

namespace fastfem{
namespace fe{

/**
 * @brief Constructor for the 1D case, in any space dimension.
 */
template<>
FESimplexP2<1, 1>::
FESimplexP2(const unsigned int n_components) :  FESimplexP<1, 1>(n_components) {

    if(n_components != 1){
        throw std::invalid_argument("FESimplexP2: only scalar functions are supported");
    }

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;
    using local_edge_id = fastfem::types::local_edge_id;

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
    vertex_dofs[(local_vertex_id){0}] = {0};
    vertex_dofs[(local_vertex_id){1}] = {2};
    edge_dofs[(local_edge_id){0,1}] = {1};

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

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;
    using local_edge_id = fastfem::types::local_edge_id;

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
    vertex_dofs[(local_vertex_id){0}] = {0};
    vertex_dofs[(local_vertex_id){1}] = {2};
    edge_dofs[(local_edge_id){0,1}] = {1};

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

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;
    using local_edge_id = fastfem::types::local_edge_id;

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
    vertex_dofs[(local_vertex_id){0}] = {0};
    vertex_dofs[(local_vertex_id){1}] = {2};
    edge_dofs[(local_edge_id){0,1}] = {1};

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

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;
    using local_edge_id = fastfem::types::local_edge_id;

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
    vertex_dofs[(local_vertex_id){0}] = {0};
    vertex_dofs[(local_vertex_id){1}] = {2};
    vertex_dofs[(local_vertex_id){2}] = {4};
    edge_dofs[(local_edge_id){0,2}] = {5};
    edge_dofs[(local_edge_id){0,1}] = {1};
    edge_dofs[(local_edge_id){1,2}] = {3};

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

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;
    using local_edge_id = fastfem::types::local_edge_id;

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
    vertex_dofs[(local_vertex_id){0}] = {0};
    vertex_dofs[(local_vertex_id){1}] = {2};
    vertex_dofs[(local_vertex_id){2}] = {4};
    edge_dofs[(local_edge_id){0,2}] = {5};
    edge_dofs[(local_edge_id){0,1}] = {1};
    edge_dofs[(local_edge_id){1,2}] = {3};

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