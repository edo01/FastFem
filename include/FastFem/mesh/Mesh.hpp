/** 
 * This file defines the basic definitions for a mesh, including vertices and elements.
 * 
 * We also assume to work with simplicial conforming meshes.
 * 
 */

#ifndef MESH_HPP
#define MESH_HPP

#include <FastFem/mesh/Geometry.hpp>
#include <vector>

namespace fastfem {
namespace mesh {

/* A Vertex is simply a point in R^3 */
template <unsigned int spacedim>
struct Vertex 
{
    Point<spacedim> point;

    bool operator==(const Vertex<spacedim> &v) const 
    {
        for (int i = 0; i < spacedim; ++i)
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
        std::copy(v, v + dim + 1, vertices);
    }

    // just check if the vertices are the same
    bool operator==(const Simplex<dim> &s) const
    {
        for (size_t i = 0; i < dim + 1; ++i)
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
    size_t vertices[dim + 1];
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