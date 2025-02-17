#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>


namespace fastfem{
namespace mesh{

template <unsigned int spacedim>
struct Point
{
    double coords[spacedim];
    
};

template <unsigned int dim, unsigned int spacedim=dim>
class Simplex
{
public:
    Simplex(const Point<spacedim> v[dim + 1])
    {
        std::copy(v, v + dim + 1, vertices);
    }

    template <typename... Args>
    Simplex(Args... args)
    {
        static_assert(sizeof...(args) == dim + 1, "Number of arguments must be equal to dim + 1");
        Point<spacedim> temp[] = {args...};
        std::copy(temp, temp + dim + 1, vertices);
    }

    // compute the volume of the simplex
    double volume() const
    {
        if(dim==0) return 0;
        if(dim==1) return vertices[1].coords[0] - vertices[0].coords[0];
        if(dim==2) return 0.5 * ((vertices[1].coords[0] - vertices[0].coords[0]) * (vertices[2].coords[1] - vertices[0].coords[1]) - (vertices[2].coords[0] - vertices[0].coords[0]) * (vertices[1].coords[1] - vertices[0].coords[1]));
        else return 0; // not implemented
    }
    
    Point<spacedim> get_centroid() const
    {
        Point<spacedim> centroid;
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

    std::vector<Point<spacedim>> get_vertices() const
    {
        return std::vector<Point<spacedim>>(vertices, vertices + dim + 1);
    }

    Point<spacedim> get_vertex(size_t i) const
    {
        return vertices[i];
    }

    // An intersection between two simplices is a simplex of lower dimension.
    // Since a simplex is completely defined by its vertices, we can simply return the vertices that are common 
    // to the two simplices.
    std::vector<Point<spacedim>> intersection(const Simplex<dim> &s) const
    {
        std::vector<Point<spacedim>> result;
        for (size_t i = 0; i < dim + 1; ++i)
        {
            for (size_t j = 0; j < dim + 1; ++j)
            {
                if (vertices[i] == s.vertices[j])
                {
                    result.push_back(vertices[i]);
                    break;
                }
            }
        }
        return result;
    }
    
    bool operator==(const Simplex<dim> &s) const
    {
        for (size_t i = 0; i < dim + 1; ++i)
        {
            if (vertices[i] != s.vertices[i])
                return false;
        }
        return true;
    }

    int n_vertices() const { return num_vertices; }
    int n_edges() const { return num_edges; }

private:
    Point<spacedim> vertices[dim + 1];
    int num_vertices = dim + 1;
    int num_edges = (dim+1)*dim/2;
};

} // namespace mesh
} // namespace fastfem

#endif // GEOMETRY_HPP