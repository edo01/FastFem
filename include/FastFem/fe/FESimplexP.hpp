#ifndef FESIMPLEXP_HPP
#define FESIMPLEXP_HPP

#include <cassert>
#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/types/CommonTypes.hpp"
#include "FastFem/linalg/FullMatrix.hpp"

namespace fastfem{
namespace fe{

/**
 * @brief Base class for Simplicial Lagrange finite elements.
 */
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP
{
    static_assert(dim > 0, "The dimension must be greater than 0");
    static_assert(dim <= spacedim, "The dimension of the FE must be less or equal to the space it lives in");
    static_assert(spacedim <= 3, "The space dimension must be less or equal to 3");

    using global_vertex_index  = fastfem::types::global_vertex_index;
    using global_edge_index    = fastfem::types::global_edge_index;
    using global_face_index    = fastfem::types::global_face_index;
    using global_cell_index    = fastfem::types::global_cell_index;

    using local_vertex_index   = fastfem::types::local_vertex_index;
    using local_edge_index     = fastfem::types::local_edge_index;
    using local_face_index     = fastfem::types::local_face_index;
    using local_cell_index     = fastfem::types::local_cell_index;

    using local_dof_index = fastfem::types::local_dof_index;
    
public:


    FESimplexP() : n_components(1) {}
    FESimplexP(int n_components) : n_components(n_components) {} 
    virtual ~FESimplexP(){};

    /**
     * @brief Get a local numbering of the degrees of freedom on a vertex, edge, face or cell.
     * 
     * @param T The simplex on which the degrees of freedom are defined.
     * @param v The global index of the sub-simplex
     * @return std::vector<local_dof_index> The local numbering of the degrees of freedom on the sub-simplex.
     */
    std::vector<local_dof_index> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_index v) const;
    std::vector<local_dof_index> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_index e) const;
    std::vector<local_dof_index> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_face_index f) const;
    std::vector<local_dof_index> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_index c) const;

    // total number of degrees of freedom
    unsigned int get_n_dofs_per_element() const { return n_dofs_per_element; }
    
    // number of degrees of freedom per cell, not including dofs on lower dimensional objects. 
    unsigned int get_n_dofs_per_cell() const { return n_dofs_per_cell; }
    
    // number of degrees of freedom per face, not including dofs on lower dimensional objects.
    unsigned int get_n_dofs_per_face() const { return n_dofs_per_face; }
    
    // number of degrees of freedom per edge, not including dofs on lower dimensional objects.
    unsigned int get_n_dofs_per_edge() const { return n_dofs_per_edge; }

    // number of degrees of freedom per vertex, not including dofs on lower dimensional objects.
    unsigned int get_n_dofs_per_vertex() const { return n_dofs_per_vertex; }

    mesh::Simplex<dim, spacedim> get_reference_simplex() const { return reference_simplex; }

    //implements an affine map from the reference simplex to the physical simplex
    mesh::Point<spacedim> get_dof_coords(mesh::Simplex<dim, spacedim> T, local_dof_index dof) const;

    virtual unsigned int get_degree() const = 0;

    virtual void compute_stiffness_loc(const mesh::Simplex<dim, spacedim> &elem, linalg::FullMatrix &matrix) const;

protected:    
    unsigned int n_components;    // number of components of the finite element
    unsigned int n_dofs_per_element; // total number of degrees of freedom of the finite element
    unsigned int n_dofs_per_cell; // number of degrees of freedom per cell
    unsigned int n_dofs_per_face; // number of degrees of freedom per face
    unsigned int n_dofs_per_edge; // number of degrees of freedom per edge
    unsigned int n_dofs_per_vertex; // number of degrees of freedom per vertex
    
    mesh::Simplex<dim, spacedim> reference_simplex; // reference simplex
    std::vector<mesh::Point<spacedim>> dofs; // points on the reference simplex that correspond to the dofs

    /**
     * @brief Tables that map a local sub-simplex to the local numbering of the degrees of freedom.
     */
    fastfem::types::local_vertex_dof_table vertex_dofs;
    fastfem::types::local_edge_dof_table edge_dofs;
    fastfem::types::local_face_dof_table face_dofs;
    fastfem::types::local_cell_dof_table cell_dofs;

    /**
     * @brief Map a global index of a sub-simplex to a local index.
     */
    local_vertex_index map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_index v) const;
    local_edge_index map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_index e) const;
    local_face_index map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_face_index f) const;
    local_cell_index map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_index c) const;
};

/**
 * @brief Basic simplician P1 elements.
 */
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP1 : public FESimplexP<dim, spacedim>
{
    
public:

    FESimplexP1(const unsigned int n_components);
    FESimplexP1() : FESimplexP1(1) {} 

    unsigned int get_degree() const override { return 1; }

    void compute_stiffness_loc(const mesh::Simplex<dim, spacedim> &elem, linalg::FullMatrix &matrix) const override;

};


// Basic simplician P2 elements
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP2 : public FESimplexP<dim, spacedim>
{
    static_assert(dim <= 2, "FESimplexP2 does not support 3D elements");

public:

    FESimplexP2(const unsigned int n_components);
    FESimplexP2() : FESimplexP2(1) {} 

    unsigned int get_degree() const override { return 2; }

    void compute_stiffness_loc(const mesh::Simplex<dim, spacedim> &elem, linalg::FullMatrix &matrix) const override;
};


// Basic simplician P3 elements 
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP3 : public FESimplexP<dim, spacedim>
{
    static_assert(dim <= 2, "FESimplexP3 does not support 3D elements");
    
public:

    FESimplexP3(const unsigned int n_components);
    FESimplexP3() : FESimplexP3(1){};

    unsigned int get_degree() const override { return 3; }

    void compute_stiffness_loc(const mesh::Simplex<dim, spacedim> &elem, linalg::FullMatrix &matrix) const override;
};

} // namespace fe
} // namespace FastFem

#endif // FESIMPLEXP_HPP