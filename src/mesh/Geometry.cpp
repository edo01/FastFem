#include "FastFem/mesh/Geometry.hpp"

#include <stdexcept>

namespace fastfem{
namespace mesh{


/**
 * ------------------------ Point class ------------------------
 */

template <unsigned int spacedim>
bool Point<spacedim>::operator==(const Point<spacedim> &p) const
{
    for (size_t i = 0; i < spacedim; ++i)
    {
        if (coords[i] != p.coords[i])
            return false;
    }
    return true;
}


/**
 * ------------------------ Simplex class ------------------------
 */

template <unsigned int dim, unsigned int spacedim>
Simplex<dim, spacedim>::Simplex(Point<spacedim> v[dim + 1])
{
    std::copy(v, v + dim + 1, vertices);
}

template <unsigned int dim, unsigned int spacedim>
double Simplex<dim, spacedim>::volume() const
{
    if(dim != spacedim)
    {
        throw std::runtime_error("volume(): The dimension of the simplex must be equal to the space it lives in");
    }
    if(dim==0) return 0;
    if(dim==1) return vertices[1].coords[0] - vertices[0].coords[0];
    if(dim==2) return 0.5 * ((vertices[1].coords[0] - vertices[0].coords[0]) * (vertices[2].coords[1] - vertices[0].coords[1]) - (vertices[2].coords[0] - vertices[0].coords[0]) * (vertices[1].coords[1] - vertices[0].coords[1]));
    else return 0; // not implemented
}

template <unsigned int dim, unsigned int spacedim>
Point<spacedim> Simplex<dim, spacedim>::get_centroid() const
{
    Point<spacedim> centroid;
    
    //add a constructor to Point class or make coords a std::array
    for(size_t i=0; i<spacedim; ++i)
    {
        centroid.coords[i] = 0;
    }

    for (size_t i = 0; i < dim + 1; ++i)
    {
        for (size_t j = 0; j < spacedim; ++j)
        {
            centroid.coords[j] += vertices[i].coords[j];
        }
    }
    for (size_t j = 0; j < spacedim; ++j)
    {
        centroid.coords[j] /= dim + 1;
    }
    return centroid;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<Point<spacedim>> Simplex<dim, spacedim>::get_vertices() const
{
    return std::vector<Point<spacedim>>(vertices, vertices + dim + 1);
}

template <unsigned int dim, unsigned int spacedim>
Point<spacedim> Simplex<dim, spacedim>::get_vertex(size_t i) const
{
    return vertices[i];
}

template <unsigned int dim, unsigned int spacedim>
std::vector<Point<spacedim>> Simplex<dim, spacedim>::intersection(const Simplex<dim, spacedim> &s) const
{
    std::vector<Point<spacedim>> result;
    for (size_t i = 0; i < dim + 1; ++i)
    {
        for (size_t j = 0; j < dim + 1; ++j)
        {
            if (vertices[i] == s.get_vertex(j))
            {
                result.push_back(vertices[i]);
                break;
            }
        }
    }
    return result;
}

template <unsigned int dim, unsigned int spacedim>
bool Simplex<dim, spacedim>::operator==(const Simplex<dim, spacedim> &s) const
{
    for (size_t i = 0; i < dim + 1; ++i)
    {
        if (vertices[i] != s.get_vertex(i))
            return false;
    }
    return true;
}

template <unsigned int dim, unsigned int spacedim>
int Simplex<dim, spacedim>::n_vertices() const { return num_vertices; }

template <unsigned int dim, unsigned int spacedim>
int Simplex<dim, spacedim>::n_edges() const { return num_edges; }

// Explicit instantiation

template struct Point<1>;
template struct Point<2>;
template struct Point<3>;

template class Simplex<0, 1>;
template class Simplex<1, 1>;
template class Simplex<1, 2>;
template class Simplex<1, 3>;
template class Simplex<2, 2>;
template class Simplex<2, 3>;
template class Simplex<3, 3>;


} // namespace mesh
} // namespace fastfem