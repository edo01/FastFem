/** 
 * This file defines the basic definitions for a mesh, including vertices and elements.
 * 
 * We also assume to work with simplicial conforming meshes.
 * 
 */

#ifndef MESH_HPP
#define MESH_HPP

#include <vector>

namespace fastfem {
namespace mesh {

/* A Vertex is simply a point in R^3 */
template <unsigned int spacedim>
struct Vertex 
{
    double coords[spacedim];

    bool operator==(const Vertex<spacedim> &v) const 
    {
        for (int i = 0; i < spacedim; ++i)
        {
            if (coords[i] != v.coords[i])
                return false;
        }
        return true;
    }
};


template <unsigned int dim>
class Simplex
{
public:
    size_t vertex_count() const { return dim + 1; }
    size_t face_count() const { return dim + 1; }

    size_t get_vertex(unsigned int i) const { return vertices[i]; }
    void set_vertex(unsigned int i, unsigned int v) { vertices[i] = v; }

    Simplex(const int v[])
    {
        for (size_t i = 0; i < dim + 1; ++i)
        {
            vertices[i] = v[i];
        }
        
        // build faces
        /* for (size_t i = 0; i < dim + 1; ++i)
        {
            Vertex<spacedim> face_vertices[dim];
            for (size_t j = 0; j < dim; ++j)
            {
                if (j < i)
                    face_vertices[j] = v[j];
                else
                    face_vertices[j] = v[j + 1];
            }
            faces[i] = Simplex<dim-1>(face_vertices);
        } */
    }

private:
    //Simplex<dim-1> faces[dim];
    size_t vertices[dim + 1];
};

// Specialization for 1-dimensional simplex (a line segment)
/* template <>
class Simplex<1>
{
public:
    size_t vertex_count() const { return 2; }
    size_t face_count() const { return 2; }

private:
    size_t vertices[2];
}; */

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
    void add_element(const Simplex<dim> &e) { elements.push_back(e); }

    void set_vertex(size_t i, const Vertex<spacedim> &v) { vertices[i] = v; }
    void set_element(size_t i, const Simplex<dim> &e) { elements[i] = e; }

    Vertex<spacedim> &get_vertex(size_t i) { return vertices[i]; }
    Simplex<dim> &get_element(size_t i) { return elements[i]; }

    void resize_vertices(size_t n) { vertices.resize(n); }

    void reserve_vertices(size_t n) { vertices.reserve(n); }
    void reserve_elements(size_t n) { elements.reserve(n); }

private:
    std::vector<Vertex<spacedim>> vertices;
    std::vector<Simplex<dim>> elements;
};

}// mesh
}// fastfem

#endif // MESH_HPP