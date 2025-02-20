#include "FastFem/fe/FESimplexP.hpp"

namespace fastfem{
namespace fe{


/**
 * @brief Constructor for the 1D case, in 1D space.
 */
template<>
FESimplexP1<1, 1>::
FESimplexP1(const unsigned int n_components) : FESimplexP<1, 1>(n_components) {

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;

    this->n_dofs_per_element = 2*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 0;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex and the local numbering of the dofs
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     * (0) [0]--------[1] (1)
     * 
     */

    // fill the vertex dofs (no dof on edges, faces or cells)
    for(int v=0; v<2; ++v){
        std::vector<local_dof_index_t> vertex_dofs_on_vertex(n_components);
        for(int j=0; j<n_components; ++j){
            vertex_dofs_on_vertex[j] = j*2 + v;
        }
        this->vertex_dofs[(local_vertex_id){v}] = vertex_dofs_on_vertex;
    }

    // set the reference simplex
    mesh::Point<1> p0 = {0};
    mesh::Point<1> p1 = {1};
    this->reference_simplex = mesh::Simplex<1>({p0, p1});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
}

/**
 * @brief Constructor for the 1D case, in 2D space.
 */
template<>
FESimplexP1<1, 2>::
FESimplexP1(const unsigned int n_components) : FESimplexP<1, 2>(n_components) {

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;

    this->n_dofs_per_element = 2*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 0;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex and the local numbering of the dofs
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     * (0) [0]--------[1] (1)
     * 
     */

    // fill the vertex dofs (no dof on edges, faces or cells)
    for(int v=0; v<2; ++v){
        std::vector<local_dof_index_t> vertex_dofs_on_vertex(n_components);
        for(int j=0; j<n_components; ++j){
            vertex_dofs_on_vertex[j] = j*2 + v;
        }
        this->vertex_dofs[(local_vertex_id){v}] = vertex_dofs_on_vertex;
    }

    // set the reference simplex
    mesh::Point<2> p0 = {0, 0};
    mesh::Point<2> p1 = {1, 0};
    this->reference_simplex = mesh::Simplex<1, 2>({p0, p1});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);

}

/**
 * @brief Constructor for the 1D case, in 3D space.
 */
template<>
FESimplexP1<1, 3>::
FESimplexP1(const unsigned int n_components) : FESimplexP<1, 3>(n_components) {

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;

    this->n_dofs_per_element = 2*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 0;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex and the local numbering of the dofs
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * 
     * (0) [0]--------[1] (1)
     * 
     */

    // fill the vertex dofs (no dof on edges, faces or cells)
    for(int v=0; v<2; ++v){
        std::vector<local_dof_index_t> vertex_dofs_on_vertex(n_components);
        for(int j=0; j<n_components; ++j){
            vertex_dofs_on_vertex[j] = j*2 + v;
        }
        this->vertex_dofs[(local_vertex_id){v}] = vertex_dofs_on_vertex;
    }

    // set the reference simplex
    mesh::Point<3> p0 = {0, 0, 0};
    mesh::Point<3> p1 = {0, 1, 0};
    this->reference_simplex = mesh::Simplex<1, 3>({p0, p1});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
}

/**
 * @brief Constructor for the 2D case, in 2D space.
 */
template<>
FESimplexP1<2, 2>::
FESimplexP1(const unsigned int n_components) : FESimplexP<2, 2>(n_components) {

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;

    this->n_dofs_per_element = 3*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 0;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex and the local numbering of the dofs
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * (2) [2]
     *     |  \
     *     |   \
     *     |    \
     *     |     \
     * (0) [0]---[1] (1)
     */

    // fill the vertex dofs (no dof on edges, faces or cells)
    for(int v=0; v<3; ++v){
        std::vector<local_dof_index_t> vertex_dofs_on_vertex(n_components);
        for(int j=0; j<n_components; ++j){
            vertex_dofs_on_vertex[j] = j*3 + v;
        }
        this->vertex_dofs[(local_vertex_id){v}] = vertex_dofs_on_vertex;
    }

    // set the reference simplex
    mesh::Point<2> p0 = {0, 0};
    mesh::Point<2> p1 = {1, 0};
    mesh::Point<2> p2 = {0, 1};
    this->reference_simplex = mesh::Simplex<2,2>({p0, p1, p2});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
}

/**
 * @brief Constructor for the 2D case, in 3D space.
 */
template<>
FESimplexP1<2, 3>::
FESimplexP1(const unsigned int n_components) : FESimplexP<2, 3>(n_components) {

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;

    this->n_dofs_per_element = 3*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 0;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex and the local numbering of the dofs
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     * (2) [2]
     *     |  \
     *     |   \
     *     |    \
     *     |     \
     * (0) [0]---[1] (1)
     */

    // fill the vertex dofs (no dof on edges, faces or cells)
    for(int v=0; v<3; ++v){
        std::vector<local_dof_index_t> vertex_dofs_on_vertex(n_components);
        for(int j=0; j<n_components; ++j){
            vertex_dofs_on_vertex[j] = j*3 + v;
        }
        this->vertex_dofs[(local_vertex_id){v}] = vertex_dofs_on_vertex;
    }

    // set the reference simplex
    mesh::Point<3> p0 = {0, 0, 0};
    mesh::Point<3> p1 = {1, 0, 0};
    mesh::Point<3> p2 = {0, 1, 0};
    this->reference_simplex = mesh::Simplex<2,3>({p0, p1, p2});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
}

/**
 * @brief Constructor for the 3D case, in 3D space.
 */
template<>
FESimplexP1<3, 3>::
FESimplexP1(const unsigned int n_components) : FESimplexP<3, 3>(n_components) {

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    using local_vertex_id = fastfem::types::local_vertex_id;

    this->n_dofs_per_element = 4*n_components;
    this->n_dofs_per_vertex = n_components;
    this->n_dofs_per_edge = 0;
    this->n_dofs_per_face = 0;
    this->n_dofs_per_cell = 0;

    // set the reference simplex and the local numbering of the dofs
    /**
     * We use modulo class of equivalence to define the numbering of the dofs. 
     *     (3) [3]
     *        / | \
     *       /  |  \
     *      /   |   \
     *     /    |    \
     * (0)[0]---[1]---[2] (2)  
     *          (1)     
     */

    // fill the vertex dofs (no dof on edges, faces or cells)
    for(int v=0; v<4; ++v){
        std::vector<local_dof_index_t> vertex_dofs_on_vertex(n_components);
        for(int j=0; j<n_components; ++j){
            vertex_dofs_on_vertex[j] = j*4 + v;
        }
        this->vertex_dofs[(local_vertex_id){v}] = vertex_dofs_on_vertex;
    }

    // set the reference simplex
    mesh::Point<3> p0 = {0, 0, 0};
    mesh::Point<3> p1 = {1, 0, 0};
    mesh::Point<3> p2 = {0, 1, 0};
    mesh::Point<3> p3 = {0, 0, 1};
    
    this->reference_simplex = mesh::Simplex<3,3>({p0, p1, p2, p3});

    // set the dofs
    this->dofs.push_back(p0);
    this->dofs.push_back(p1);
    this->dofs.push_back(p2);
    this->dofs.push_back(p3);

}


} // namespace fe  
} // namespace fastfem
