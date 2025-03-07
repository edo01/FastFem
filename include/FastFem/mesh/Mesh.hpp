/** 
 * This file defines the basic definitions for a mesh, including vertices and elements.
 * 
 * We also assume to work with simplicial conforming meshes.
 * 
 */
#ifndef MESH_HPP
#define MESH_HPP

#include <stdint.h>
#include <unordered_map>

#include "FastFem/mesh/Geometry.hpp"
#include "FastFem/common/hash_table.h"
#include "FastFem/types/CommonTypes.hpp"

namespace fastfem {
namespace mesh {

/* A Vertex is simply a point in R^3 */
template <unsigned int spacedim>
struct Vertex 
{
    Point<spacedim> point;

    bool operator==(const Vertex<spacedim> &v) const;
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
    static_assert(dim >= 0, "The dimension of the simplex must be equal or greater than 0");
    static_assert(dim <= spacedim, "The dimension of the simplex must be less or equal to the space it lives in");
    static_assert(spacedim <= 3, "The space dimension must be less or equal to 3");

    using global_vertex_index  = fastfem::types::global_vertex_index;
    using global_edge_index    = fastfem::types::global_edge_index;
    using global_face_index    = fastfem::types::global_face_index;
    using global_cell_index    = fastfem::types::global_cell_index;

    using local_vertex_index   = fastfem::types::local_vertex_index;

public:

    static constexpr unsigned int n_vertices = dim + 1;
    static constexpr unsigned int n_edges = num_m_subsimplex_on_n_simplex(1, dim);
    static constexpr unsigned int n_faces = num_m_subsimplex_on_n_simplex(2, dim);
    static constexpr unsigned int n_cells = num_m_subsimplex_on_n_simplex(3, dim);

    MeshSimplex(const global_vertex_index v[dim + 1]);

    /**
     * Get the indices of the vertices of the simplex. 
     */
    const std::array<global_vertex_index, n_vertices> get_vertex_indices() const;

    /**
     * Get the indices of the edges of the simplex. Each edge is represented by an
     * ordered pair of indices of the vertices.
     */
    const std::array<global_edge_index, n_edges> get_edges_indices() const;

    /**
     * Get the indices of the faces of the simplex, if present. Each face is represented by an
     * ordered triple of indices of the vertices. 
     */
    const std::array<global_face_index, n_faces> get_faces_indices() const; 

    /**
     * Get the indices of the cells of the simplex, if present. Each cell is represented by an
     * ordered quadruple of indices of the vertices.
     */
    const std::array<global_cell_index, n_cells> get_cell_indices() const;

    bool operator==(const MeshSimplex<dim, spacedim> &s) const;

    inline void set_vertex(local_vertex_index i, global_vertex_index v) { vertices[i] = v; }
    inline global_vertex_index get_vertex(local_vertex_index i) const { return vertices[i]; }
    inline local_vertex_index vertex_count() const { return dim + 1; }

private:
    std::array<global_vertex_index, dim + 1> vertices;
};

/**
 *  This domain, and the mesh that covers it, represents a _dim-dimensional_ manifold
 *  and lives in _spacedim_ spatial dimensions, where dim and spacedim are the template 
 *  arguments of this class. For example the surface of a sphere is a 2-dimensional manifold
 *  that lives in 3-dimensional space.
 * 
 */
template <unsigned int dim, unsigned int spacedim = dim>
class Mesh
{
    static_assert(dim > 0, "The dimension of the mesh must be greater than 0");
    static_assert(dim <= spacedim, "The dimension of the mesh must be less or equal to the space it lives in");
    static_assert(spacedim <= 3, "The space dimension must be less or equal to 3");

    using global_vertex_index  = fastfem::types::global_vertex_index;
    using global_element_index = fastfem::types::global_element_index;
    using boundary_index       = fastfem::types::boundary_index;
    using local_vertex_index   = fastfem::types::local_vertex_index;

public:

    // given a MeshSimplex in the mesh, we return the corresponding Simplex in the space
    Simplex<dim, spacedim> get_Simplex(MeshSimplex<dim, spacedim> s) const
    {
        Point<spacedim> p[dim + 1];
        for (local_vertex_index i = 0; i < dim + 1; ++i)
        {
            //p[i] = vertices[s.get_vertex(i)];

            p[i] = vertices[s.get_vertex(i)].point;

        }
        return Simplex<dim, spacedim>(p);
    }

    inline global_vertex_index vtx_count() const { return vertices.size(); }
    inline global_element_index elem_count() const { return elements.size(); }

    inline void add_vertex(const Vertex<spacedim> &v) { vertices.push_back(v); }
    inline void add_element(const MeshSimplex<dim, spacedim> &e) { elements.push_back(e); }

    inline void set_vertex(global_vertex_index i, const Vertex<spacedim> &v) { vertices[i] = v; }
    inline void set_element(global_element_index i, const MeshSimplex<dim, spacedim> &e) { elements[i] = e; }

    inline Vertex<spacedim> &get_vertex(global_vertex_index i) { return vertices[i]; }
    inline MeshSimplex<dim, spacedim> &get_mesh_element(global_element_index i) { return elements[i]; }

    inline void resize_vertices(size_t n) { vertices.resize(n); }
    inline void reserve_vertices(size_t n) { vertices.reserve(n); }
    inline void reserve_elements(size_t n) { elements.reserve(n); }

    inline void add_boundary_element(boundary_index tag, const MeshSimplex<dim - 1, spacedim> &e) { map_boundary_elements[tag].push_back(e); }
    inline MeshSimplex<dim - 1, spacedim> &get_boundary_element(boundary_index tag, global_element_index i) { return map_boundary_elements[tag][i]; }
    inline void reserve_boundary_elements(boundary_index tag, size_t n) { map_boundary_elements[tag].reserve(n); }
    inline size_t boundary_elem_count(boundary_index tag) const { return map_boundary_elements.at(tag).size(); }
    

    /**
     * ITERATORS
     */
    inline auto vtx_begin() { return vertices.begin(); }
    inline auto vtx_end() { return vertices.end(); }

    inline auto elem_begin() { return elements.begin(); }
    inline auto elem_end() { return elements.end(); }

    inline auto elem_begin() const { return elements.begin(); }
    inline auto elem_end() const { return elements.end(); }

    // iterator on the boundaries
    inline auto boundary_begin() const { return map_boundary_elements.begin(); }
    inline auto boundary_end() const { return map_boundary_elements.end(); }

private:
    std::vector<Vertex<spacedim>> vertices;
    std::vector<MeshSimplex<dim, spacedim>> elements;

    /*
     * We now store the elements of the boundary. The boundary is partitioned in different tags, 
     * each tag is associated to a set of simplices of dimension dim-1.
     */
    std::unordered_map<boundary_index, std::vector<MeshSimplex<dim - 1, spacedim>>> map_boundary_elements;

};

}// mesh
}// fastfem

#endif // MESH_HPP