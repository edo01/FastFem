/** 
 * This file defines the basic definitions for a mesh, including vertices and elements.
 * 
 * We also assume to work with simplicial conforming meshes.
 * 
 */
#ifndef MESH_HPP
#define MESH_HPP

constexpr inline size_t binom(size_t n, size_t k) noexcept
{
    return
      (        k> n  )? 0 :          // out of range
      (k==0 || k==n  )? 1 :          // edge
      (k==1 || k==n-1)? n :          // first
      (     k+k < n  )?              // recursive:
      (binom(n-1,k-1) * n)/k :       //  path to k=1   is faster
      (binom(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
}

#define get_m_faces_on_n_simplex(m, n) (binom(n+1, m+1))

#include <stdint.h>

#include "FastFem/mesh/Geometry.hpp"
#include "FastFem/mesh/VertexHasher.hpp"
#include "FastFem/common/hash_table.h"
#include "FastFem/types/CommonTypes.hpp"

namespace fastfem {
namespace mesh {

/* A Vertex is simply a point in R^3 */
template <unsigned int spacedim>
struct Vertex 
{
    Point<spacedim> point;

    bool operator==(const Vertex<spacedim> &v) const 
    {
        for (unsigned int i = 0; i < spacedim; ++i)
        {
            if (point.coords[i] != v.point.coords[i])
                return false;
        }
        return true;
    }
};

/**
 * We separate the definition of a simplex from the mesh. A MeshSimplex is a utility class that is used to represent a simplex 
 * in a mesh. It is a container for the indices of the vertices of the simplex and 
 */
template <unsigned int dim, unsigned int spacedim=dim>
class MeshSimplex
{
    /*
     * Some static assertions to make sure that the template arguments are correct
     * based on the following constraints:
     *  - we don't want to have a simplex with a dimension greater than the space it lives in
     *  - the dimension of the simplex must be greater than 0
     *  - we don't handle spaces with dimension greater than 3
     */
    static_assert(dim > 0, "The dimension of the simplex must be greater than 0");
    static_assert(dim <= spacedim, "The dimension of the simplex must be less or equal to the space it lives in");
    static_assert(spacedim <= 3, "The space dimension must be less or equal to 3");


public:

    static constexpr unsigned int n_vertices = dim + 1;
    static constexpr unsigned int n_edges = get_m_faces_on_n_simplex(1, dim);
    static constexpr unsigned int n_faces = get_m_faces_on_n_simplex(2, dim);
    static constexpr unsigned int n_cells = get_m_faces_on_n_simplex(3, dim);

    MeshSimplex(const size_t v[dim + 1])
    {
        std::copy(v, v + dim + 1, vertices.begin());
    }

    /**
     * Get the indices of the vertices of the simplex. 
     */
    const std::array<fastfem::types::simplex_index<0>, n_vertices> get_vertex_indices() const {
        std::array<fastfem::types::simplex_index<0>, n_vertices> v_indices;

        unsigned int vertex_count = 0;

        for (size_t i = 0; i < dim + 1; ++i)
        {
            std::array<size_t, 1> v = {vertices[i]};
            v_indices[vertex_count++] = v;
        }
        return v_indices;
    }

    /**
     * Get the indices of the edges of the simplex. Each edge is represented by an
     * ordered pair of indices of the vertices.
     */
    const std::array<fastfem::types::simplex_index<1>, n_edges> get_edges_indices() const {
        std::array<fastfem::types::simplex_index<1>, n_edges> e_indices;
        unsigned int edge_count = 0;

        for (size_t i = 0; i < dim + 1; ++i)
        {
            for (size_t j = i + 1; j < dim + 1; ++j)
            {
                // always store the indices in ascending order
                if(vertices[i] < vertices[j])
                    e_indices[edge_count++] = {vertices[i], vertices[j]};
                else
                    e_indices[edge_count++] = {vertices[j], vertices[i]};
            }
        }
        
        assert(edge_count == n_edges);

        return e_indices;
    }    

    /**
     * Get the indices of the faces of the simplex. Each face is represented by an
     * ordered triple of indices of the vertices. 
     */
    const std::array<fastfem::types::simplex_index<2>, n_faces> get_faces_indices() const {
        if(n_faces == 0) return {}; // in 1D there are no faces

        std::array<fastfem::types::simplex_index<2>, n_faces> f_indices;

        unsigned int face_count = 0;

        for (size_t i = 0; i < dim + 1; ++i){
            for (size_t j = i + 1; j < dim + 1; ++j){
                for (size_t k = j + 1; k < dim + 1; ++k){
                    fastfem::types::simplex_index<2> face = {vertices[i], vertices[j], vertices[k]};
                    std::sort(face.begin(), face.end());
                    // copy starting from the face_count index
                    f_indices[face_count++] = face;
                }
            }
        }

        assert(face_count == n_faces);

        return f_indices;
    }

    const std::array<fastfem::types::simplex_index<3>, n_cells> get_cell_indices() const {
        if(n_cells == 0) return {}; // in 1D and 2D there are no cells
        
        std::array<fastfem::types::simplex_index<3>, n_cells> c_indices;

        unsigned int cell_count = 0;

        for (size_t i = 0; i < dim + 1; ++i){
            for (size_t j = i + 1; j < dim + 1; ++j){
                for (size_t k = j + 1; k < dim + 1; ++k){
                    for (size_t l = k + 1; l < dim + 1; ++l){
                        fastfem::types::simplex_index<3> cell = {vertices[i], vertices[j], vertices[k], vertices[l]};
                        std::sort(cell.begin(), cell.end());
                        c_indices[cell_count++] = cell;                        
                    }
                }
            }
        }

        assert(cell_count == n_cells);

        return c_indices;
    }

    // just check if the vertices are the same
    bool operator==(const MeshSimplex<dim, spacedim> &s) const
    {
        for (unsigned int i = 0; i < dim + 1; ++i)
        {
            if (vertices[i] != s.vertices[i])
                return false;
        }
        return true;
    }

    void set_vertex(size_t i, size_t v) { vertices[i] = v; }
    size_t get_vertex(size_t i) const { return vertices[i]; }
    int vertex_count() const { return dim + 1; }

private:
    std::array<size_t, dim + 1> vertices;
};

/**
 *  This domain, and the mesh that covers it, represents a _dim-dimensional_ manifold
 *  and lives in _spacedim_ spatial dimensions, where dim and spacedim are the template 
 *  arguments of this class. For example the surface of a sphere is a 2-dimensional manifold
 *  that lives in 3-dimensional space.
 */
template <unsigned int dim, unsigned int spacedim = dim>
class Mesh
{

public:
    size_t vtx_count() const { return vertices.size(); }
    size_t elem_count() const { return elements.size(); }

    void add_vertex(const Vertex<spacedim> &v) { vertices.push_back(v); }
    void add_element(const MeshSimplex<dim, spacedim> &e) { elements.push_back(e); }

    void set_vertex(size_t i, const Vertex<spacedim> &v) { vertices[i] = v; }
    void set_element(size_t i, const MeshSimplex<dim, spacedim> &e) { elements[i] = e; }

    Vertex<spacedim> &get_vertex(size_t i) { return vertices[i]; }
    MeshSimplex<dim, spacedim> &get_mesh_element(size_t i) { return elements[i]; }

    void resize_vertices(size_t n) { vertices.resize(n); }

    void reserve_vertices(size_t n) { vertices.reserve(n); }
    void reserve_elements(size_t n) { elements.reserve(n); }

    auto vtx_begin() { return vertices.begin(); }
    auto vtx_end() { return vertices.end(); }

    auto elem_begin() { return elements.begin(); }
    auto elem_end() { return elements.end(); }

    auto elem_begin() const { return elements.begin(); }
    auto elem_end() const { return elements.end(); }

    void add_boundary(size_t i) { boundary.push_back(i); }
    size_t boundary_count() const { return boundary.size(); }
    Vertex<spacedim> get_boundary_vertex(size_t i) const { return vertices[boundary[i]]; }
    
    Simplex<dim, spacedim> get_Simplex(MeshSimplex<dim, spacedim> s) const
    {
        Point<spacedim> p[dim + 1];
        for (size_t i = 0; i < dim + 1; ++i)
        {
            p[i] = vertices[s.get_vertex(i)];
        }
        return Simplex<dim, spacedim>(p);
    }

 

    private:
    std::vector<Vertex<spacedim>> vertices;
    std::vector<MeshSimplex<dim, spacedim>> elements;

    // Pointer to the boundary vertices
    std::vector<size_t> boundary;

};

}// mesh
}// fastfem

#endif // MESH_HPP