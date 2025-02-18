/**
 * Implementation of a scalar Lagrange finite element Pr that yields the 
 * finite element space of continuous, piecewise polynomials of degree r
 * defined on a simplex mesh.
 * 
 * The DOFs have a local internal numbering.
 */
#ifndef FESIMPLEXP_HPP
#define FESIMPLEXP_HPP

#include <FastFem/mesh/Mesh.hpp>
#include <cassert>

namespace fastfem{
namespace fe{

template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP
{
    
public:

    typedef unsigned int local_index_t;

    FESimplexP() : n_components(1) {}
    FESimplexP(int n_components) : n_components(n_components) {} 
    virtual ~FESimplexP(){};

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

    // returns the point in the mesh that corresponds to the local i-th dof on the reference simplex
    virtual mesh::Point<spacedim> get_dofs_on_simplex(local_index_t i, mesh::Simplex<dim> ) const = 0;

    mesh::Simplex<dim, spacedim> get_reference_simplex() const { return reference_simplex; }

    //std::vector<local_index_t> get_ref_dofs_on_cell(

protected:    
    unsigned int n_components;    // number of components of the finite element
    unsigned int n_dofs_per_element; // total number of degrees of freedom of the finite element
    unsigned int n_dofs_per_cell; // number of degrees of freedom per cell
    unsigned int n_dofs_per_edge; // number of degrees of freedom per edge
    unsigned int n_dofs_per_face; // number of degrees of freedom per face
    unsigned int n_dofs_per_vertex; // number of degrees of freedom per vertex
    mesh::Simplex<dim> reference_simplex; // reference simplex

    std::vector<mesh::Point<spacedim>> dofs; // points on the reference simplex that correspond to the dofs






};

// Basic simplician P1 elements in 2D space
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP1 : public FESimplexP<dim, spacedim>
{
public:

    FESimplexP1(const unsigned int n_components) : FESimplexP<2>(n_components) {
        // compute the number of degrees of freedom
        this->n_dofs_per_element = 3*n_components;
        this->n_dofs_per_vertex = n_components;
        this->n_dofs_per_edge = 0;
        this->n_dofs_per_face = 0;
        this->n_dofs_per_cell = 0;

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
        this->reference_simplex = mesh::Simplex<2>({p0, p1, p2});

        // set the dofs
        this->dofs.push_back(p0);
        this->dofs.push_back(p1);
        this->dofs.push_back(p2);
    }

    FESimplexP1() : FESimplexP1(1) {} 
    
    mesh::Point<2> get_dofs_on_simplex( typename FESimplexP<2>::local_index_t i, mesh::Simplex<2> s) const override {
        assert(i < this->n_dofs_per_cell);
        return this->dofs[i%this->n_components];
    } 
};


// Basic simplician P2 elements in 2D space
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP2 : public FESimplexP<dim, spacedim>
{
public:

    FESimplexP2(const unsigned int n_components) : FESimplexP<2>(n_components) {
        // compute the number of degrees of freedom
        this->n_dofs_per_element = 6*n_components;
        this->n_dofs_per_vertex = n_components;
        this->n_dofs_per_edge = n_components;
        this->n_dofs_per_face = 0;
        this->n_dofs_per_cell = 0;

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
        mesh::Point<2> p3 = {0.5, 0};
        mesh::Point<2> p4 = {0.5, 0.5};
        mesh::Point<2> p5 = {0, 0.5};
        this->reference_simplex = mesh::Simplex<2>({p0, p1, p2});

        // set the dofs
        this->dofs.push_back(p0);
        this->dofs.push_back(p1);
        this->dofs.push_back(p2);
        this->dofs.push_back(p3);
        this->dofs.push_back(p4);
        this->dofs.push_back(p5);
    }

    FESimplexP2() : FESimplexP2(1) {} 
    
    mesh::Point<2> get_dofs_on_simplex( typename FESimplexP<2>::local_index_t i, mesh::Simplex<2> s) const override {
        assert(i < this->n_dofs_per_cell);
        return this->dofs[i%this->n_components];
    } 
};


// Basic simplician P3 elements in 2D space
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP3 : public FESimplexP<dim, spacedim>
{
public:

    FESimplexP3(const unsigned int n_components) : FESimplexP<2>(n_components) {
        // compute the number of degrees of freedom
        this->n_dofs_per_element = 10*n_components;
        this->n_dofs_per_vertex  = n_components;
        this->n_dofs_per_edge    = 2*n_components;
        this->n_dofs_per_face    = n_components;
        this->n_dofs_per_cell    = 0;

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
        mesh::Point<2> p3 = {1.0/3.0, 0};
        mesh::Point<2> p4 = {2.0/3.0, 0};
        mesh::Point<2> p5 = {2.0/3.0, 1.0/3.0};
        mesh::Point<2> p6 = {1.0/3.0, 2.0/3.0};
        mesh::Point<2> p7 = {0, 2.0/3.0};
        mesh::Point<2> p8 = {0, 1.0/3.0};
        mesh::Point<2> p9 = {1.0/3.0, 1.0/3.0};
        this->reference_simplex = mesh::Simplex<2>({p0, p1, p2});

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
    }

    FESimplexP3() : FESimplexP3(1) {} 
    
    mesh::Point<2> get_dofs_on_simplex( typename FESimplexP<2>::local_index_t i, mesh::Simplex<2> s) const override {
        assert(i < this->n_dofs_per_cell);
        return this->dofs[i%this->n_components];
    } 
};


} // namespace fe
} // namespace FastFem

#endif // FESIMPLEXP_HPP