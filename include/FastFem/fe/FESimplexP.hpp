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

    using global_vertex_id = fastfem::types::global_vertex_id;
    using global_edge_id = fastfem::types::global_edge_id;
    using global_face_id = fastfem::types::global_face_id;
    using global_cell_id = fastfem::types::global_cell_id;

    using local_vertex_id = fastfem::types::local_vertex_id;
    using local_edge_id   = fastfem::types::local_edge_id;
    using local_face_id   = fastfem::types::local_face_id;
    using local_cell_id   = fastfem::types::local_cell_id; 

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    
public:


    FESimplexP() : n_components(1) {}
    FESimplexP(int n_components) : n_components(n_components) {} 
    virtual ~FESimplexP(){};

    std::vector<local_dof_index_t> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_id v) {
        std::vector<local_dof_index_t> local_dofs_on_vertex(n_dofs_per_vertex); //è un array !!

        // we retrieve the local vertex id in order to get the dofs associated to the vertex in the table
        local_vertex_id local_vertex = map_global_simplex_to_local(T, v);

        // now with the local_vertex we can get the dofs associated to the vertex 
        // preloaded in the constructor
        const auto& dofs = vertex_dofs.at(local_vertex);
        std::copy(dofs.begin(), dofs.end(), local_dofs_on_vertex.begin());
               
        return local_dofs_on_vertex;
    }

    std::vector<local_dof_index_t> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_id e) {
        std::vector<local_dof_index_t> local_dofs_on_edge(n_dofs_per_edge);
        local_edge_id local_edge = map_global_simplex_to_local(T, e);

        // now with the local_edge we can get the dofs associated to the edge 
        // preloaded in the constructor
        const auto& dofs = edge_dofs.at(local_edge);
        std::copy(dofs.begin(), dofs.end(), local_dofs_on_edge.begin());
               
        return local_dofs_on_edge;
    }

    std::vector<local_dof_index_t> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_face_id f) {
        std::vector<local_dof_index_t> local_dofs_on_face(n_dofs_per_face);
        local_face_id local_face = map_global_simplex_to_local(T, f);

        // now with the local_face we can get the dofs associated to the face 
        // preloaded in the constructor
        const auto& dofs = face_dofs.at(local_face);
        std::copy(dofs.begin(), dofs.end(), local_dofs_on_face.begin());
               
        return local_dofs_on_face;
    }

    std::vector<local_dof_index_t> get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_id c) {
        std::vector<local_dof_index_t> local_dofs_on_cell(n_dofs_per_cell);
        local_cell_id local_cell = map_global_simplex_to_local(T, c);

        // now with the local_cell we can get the dofs associated to the cell 
        // preloaded in the constructor
        const auto& dofs = cell_dofs.at(local_cell);
        std::copy(dofs.begin(), dofs.end(), local_dofs_on_cell.begin());
               
        return local_dofs_on_cell;
    }


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

protected:    
    unsigned int n_components;    // number of components of the finite element
    unsigned int n_dofs_per_element; // total number of degrees of freedom of the finite element
    unsigned int n_dofs_per_cell; // number of degrees of freedom per cell
    unsigned int n_dofs_per_edge; // number of degrees of freedom per edge
    unsigned int n_dofs_per_face; // number of degrees of freedom per face
    unsigned int n_dofs_per_vertex; // number of degrees of freedom per vertex
    mesh::Simplex<dim> reference_simplex; // reference simplex

    std::vector<mesh::Point<spacedim>> dofs; // points on the reference simplex that correspond to the dofs

    fastfem::types::local_dof_table<0> vertex_dofs;
    fastfem::types::local_dof_table<1> edge_dofs;
    fastfem::types::local_dof_table<2> face_dofs;
    fastfem::types::local_dof_table<3> cell_dofs;

    local_vertex_id map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_id v) {
        constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
        for(unsigned int i=0; i<n_vertices; ++i){
            if(T.get_vertex(i) == v[0]) return {i}; // sappiamo che darà warning
        } 
        assert(false);
    }

    local_edge_id map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_id e) {
        constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;

        for(int i=0; i<n_vertices; ++i){
            for(int j=i+1; j<n_vertices; ++j){

                global_edge_id e_simplex = {T.get_vertex(i), T.get_vertex(j)};
                std::sort(e_simplex.begin(), e_simplex.end());

                if(e_simplex==e) {
                    return {i,j}; // sappiamo che darà warning
                }
            }
        }  
        assert(false);
    }

    local_face_id map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_face_id f) {
        constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;

        for(int i=0; i<n_vertices; ++i){
            for(int j=i+1; j<n_vertices; ++j){
                for(int k=j+1; k<n_vertices; ++k){
                    global_face_id f_simplex = {T.get_vertex(i), T.get_vertex(j), T.get_vertex(k)};
                    std::sort(f_simplex.begin(), f_simplex.end());
                    if(f_simplex==f) {
                        return {i,j,k}; // sappiamo che darà warning
                    }
                }
            }
        } 
        assert(false);
    }

    local_cell_id map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_id c) {
        constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;

        for(int i=0; i<n_vertices; ++i){
            for(int j=i+1; j<n_vertices; ++j){
                for(int k=j+1; k<n_vertices; ++k){
                    for(int l=k+1; l<n_vertices; ++l){
                        global_cell_id c_simplex = {T.get_vertex(i), T.get_vertex(j), T.get_vertex(k), T.get_vertex(l)};
                        std::sort(c_simplex.begin(), c_simplex.end());
                        if(c_simplex==c) {
                            return {i,j,k,l}; // sappiamo che darà warning
                        }
                    }
                }
            }
        } 
        assert(false);
    }
};

// Basic simplician P1 elements in 2D space
template <unsigned int dim, unsigned int spacedim=dim>
class FESimplexP1 : public FESimplexP<dim, spacedim>
{

    using local_dof_index_t = fastfem::types::local_dof_index_t;
    
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

        // fill the tables
        for(int i=0; i<dim+1; ++i){
            // fill the vertex dofs
            std::vector<local_dof_index_t> vertex_dofs_on_vertex(n_components);
            for(int j=0; j<n_components; ++j){
                vertex_dofs_on_vertex[j] = i*n_components + j;
            }
            this->vertex_dofs[{i}] = vertex_dofs_on_vertex;
        }
    }

    FESimplexP1() : FESimplexP1(1) {} 

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
         * [4]
         * |  \
         * |   \
         * |    \
         * [5]  [3]
         * |       \
         * |        \
         * |         \
         * [0]--[1]--[2]
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

        // fill the tables
        this->vertex_dofs[{0}] = {0};
        this->vertex_dofs[{1}] = {2};
        this->vertex_dofs[{2}] = {4};
        this->edge_dofs[{0,1}] = {1};
        this->edge_dofs[{1,2}] = {3};
        this->edge_dofs[{0,2}] = {5};
    }

    FESimplexP2() : FESimplexP2(1) {} 

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
         * [6]
         * |  \
         * |   \
         * |    \
         * [7]   [5]
         * |      \
         * |       \
         * |        \
         * [8] [9]  [4]
         * |          \
         * |           \
         * |            \
         * |             \
         * |              \
         * [0]--[1]--[2]--[3]
         * 
         
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

        // fill the tables
        this->vertex_dofs[{0}] = {0};
        this->vertex_dofs[{1}] = {3};
        this->vertex_dofs[{2}] = {6};
        this->edge_dofs[{0,1}] = {1, 2};
        this->edge_dofs[{1,2}] = {4, 5};
        this->edge_dofs[{0,2}] = {7, 8};
        this->face_dofs[{0,1,2}] = {9};
    }

    FESimplexP3() : FESimplexP3(1) {} 
};


} // namespace fe
} // namespace FastFem

#endif // FESIMPLEXP_HPP