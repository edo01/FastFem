/**
 * Implementation of a scalar Lagrange finite element Pr that yields the 
 * finite element space of continuous, piecewise polynomials of degree r
 * defined on a simplex mesh.
 * 
 * The DOFs have a local internal numbering.
 */
#ifndef FESIMPLEX_HPP
#define FESIMPLEX

#include <FastFem/mesh/Mesh.hpp>

namespace fastfem{
namespace fe{

template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexBase
{
public:

    typedef local_index_t unsigned int;

    virtual ~FESimplexBase() = default;

    unsigned int get_n_dofs_per_cell() const { return n_dofs_per_cell; }
    unsigned int get_n_dofs_per_edge() const { return n_dofs_per_edge; }
    unsigned int get_n_dofs_per_face() const { return n_dofs_per_face; }
    unsigned int get_n_dofs_per_vertex() const { return n_dofs_per_vertex; }

    // returns the point in the mesh that corresponds to the local i-th dof on the reference simplex
    mesh::Point<spacedim> get_dofs_on_simplex(local_index_t i, mesh::Simplex<dim> ) const = 0;

    mesh::Simplex<dim> get_reference_simplex() const { return reference_simplex; }

    

private:    
    unsigned int n_dofs_per_cell; // number of degrees of freedom per cell
    unsigned int n_dofs_per_edge; // number of degrees of freedom per edge
    unsigned int n_dofs_per_face; // number of degrees of freedom per face
    unsigned int n_dofs_per_vertex; // number of degrees of freedom per vertex
    mesh::Simplex<dim> reference_simplex; // reference simplex

    mesh::Vertex<spacedim> dofs[dim+1]; // degrees of freedom on the reference simplex

};


/* template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP : public FESimplexBase<dim, spacedim>
{
public:
    FESimplexP(unsigned int r) : degree(r)
    {
        
    }

private:
    unsigned int degree; // polynomial degree
}; */

} // namespace fe
} // namespace FastFem

#endif // FE_BASE_HPP