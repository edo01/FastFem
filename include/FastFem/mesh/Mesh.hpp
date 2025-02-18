/** 
 * This file defines the basic definitions for a mesh, including vertices and elements.
 * 
 * We also assume to work with simplicial conforming meshes.
 * 
 */

#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
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

public:
    MeshSimplex(const size_t v[dim + 1])
    {
        std::copy(v, v + dim + 1, vertices.begin());
    }

    // could be implemented with std::array
    const std::vector<fastfem::types::simplex_index<0>> get_vertex_indices() const {
        std::vector<fastfem::types::simplex_index<0>> v_indices;
        for (size_t i = 0; i < dim + 1; ++i)
        {
            std::array<size_t, 1> v = {vertices[i]};
            v_indices.push_back(v);
        }
        return v_indices;
    }

    const std::vector<fastfem::types::simplex_index<1>> get_edges_indices() const {
        std::vector<fastfem::types::simplex_index<1>> e_indices;
        for (size_t i = 0; i < dim + 1; ++i)
        {
            for (size_t j = i + 1; j < dim + 1; ++j)
            {
                // always store the indices in ascending order
                if(vertices[i] < vertices[j])
                    e_indices.push_back({vertices[i], vertices[j]});
                else
                    e_indices.push_back({vertices[j], vertices[i]});
            }
        }
        return e_indices;
    }    

    // not tested
    const std::vector<fastfem::types::simplex_index<2>> get_faces_indices() const {
        std::vector<fastfem::types::simplex_index<2>> f_indices;
        for (size_t i = 0; i < dim + 1; ++i)
        {
            for (size_t j = i + 1; j < dim + 1; ++j)
            {
                for (size_t k = j + 1; k < dim + 1; ++k)
                {
                    // always store the indices in ascending order
                    if(vertices[i] < vertices[j] && vertices[j] < vertices[k])
                        f_indices.push_back({vertices[i], vertices[j], vertices[k]});
                    else if(vertices[i] < vertices[k] && vertices[k] < vertices[j])
                        f_indices.push_back({vertices[i], vertices[k], vertices[j]});
                    else if(vertices[j] < vertices[i] && vertices[i] < vertices[k])
                        f_indices.push_back({vertices[j], vertices[i], vertices[k]});
                    else if(vertices[j] < vertices[k] && vertices[k] < vertices[i])
                        f_indices.push_back({vertices[j], vertices[k], vertices[i]});
                    else if(vertices[k] < vertices[i] && vertices[i] < vertices[j])
                        f_indices.push_back({vertices[k], vertices[i], vertices[j]});
                    else if(vertices[k] < vertices[j] && vertices[j] < vertices[i])
                        f_indices.push_back({vertices[k], vertices[j], vertices[i]});
                }
            }
        }
        return f_indices;
    }

    const fastfem::types::simplex_index<dim> get_element_indices() const {
        fastfem::types::simplex_index<dim> e_indices;
        std::copy(vertices, vertices + dim + 1, e_indices.begin());
        std::sort(e_indices.begin(), e_indices.end());
        return e_indices;
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