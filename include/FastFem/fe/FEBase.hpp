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
#include <cassert>

namespace fastfem{
namespace fe{

template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexBase
{
public:

    typedef local_index_t unsigned int;

    virtual FESimplexBase() : n_components(1) {}
    virtual FESimplexBase(int n_components) : n_components(n_components) {}
    virtual ~FESimplexBase() = default;

    unsigned int get_n_dofs_per_cell() const { return n_dofs_per_cell; }
    unsigned int get_n_dofs_per_edge() const { return n_dofs_per_edge; }
    unsigned int get_n_dofs_per_face() const { return n_dofs_per_face; }
    unsigned int get_n_dofs_per_vertex() const { return n_dofs_per_vertex; }

    // returns the point in the mesh that corresponds to the local i-th dof on the reference simplex
    virtual mesh::Point<spacedim> get_dofs_on_simplex(local_index_t i, mesh::Simplex<dim> ) const = 0;

    mesh::Simplex<dim, spacedim> get_reference_simplex() const { return reference_simplex; }

    

protected:    
    unsigned int n_components; // number of components of the finite element
    unsigned int n_dofs_per_cell; // number of degrees of freedom per cell
    unsigned int n_dofs_per_edge; // number of degrees of freedom per edge
    unsigned int n_dofs_per_face; // number of degrees of freedom per face
    unsigned int n_dofs_per_vertex; // number of degrees of freedom per vertex
    mesh::Simplex<dim> reference_simplex; // reference simplex

    mesh::Vertex<spacedim> dofs[dim+1]; // degrees of freedom on the reference simplex

};

// Basic simplician P1 elements in 2D space
class FESimplexP1 : public FESimplexBase<2>
{
public:

    FESimplexP1(const unsigned int n_components) : FESimplexBase<2>(n_components) {
        // compute the number of degrees of freedom
        n_dofs_per_cell = 3*n_components;
        n_dofs_per_edge = 0;
        n_dofs_per_face = 0;
        n_dofs_per_vertex = n_components;

        // set the reference simplex
        /**
         * We use modulo class of equivalence to define the numbering of the dofs. 
         * [2]
         * |  \
         * |   \
         * |    \
         * |     \
         * [0]---[1]
         */
        mesh::Point<2> p0 = {0, 0};
        mesh::Point<2> p1 = {1, 0};
        mesh::Point<2> p2 = {0, 1};
        reference_simplex = mesh::Simplex<2>({p0, p1, p2});
    }

    FESimplexP1() : FESimplexP1(1) {} 
    
    mesh::Point<2> get_dofs_on_simplex(local_index_t i, mesh::Simplex<2> s) const override {
        assert(i < n_dofs_per_cell);
        return s.get_vertex(i%n_components);
    }
    
    
};

} // namespace fe
} // namespace FastFem

#endif // FE_BASE_HPP