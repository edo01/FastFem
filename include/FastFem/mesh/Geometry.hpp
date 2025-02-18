#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>


namespace fastfem{
namespace mesh{

template <unsigned int spacedim>
struct Point
{
    double coords[spacedim];

    bool operator==(const Point<spacedim> &p) const;

    inline bool operator!=(const Point<spacedim> &p) const { return !(*this == p); }
    
};

template <unsigned int dim, unsigned int spacedim=dim>
class Simplex
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

public:
    Simplex() = default;
    
    Simplex(const Point<spacedim> v[dim + 1]);

    template <typename... Args>
    Simplex(Args... args)
    {
        static_assert(sizeof...(args) == dim + 1, "Number of arguments must be equal to dim + 1");
        Point<spacedim> temp[] = {args...};
        std::copy(temp, temp + dim + 1, vertices);
    }
    

    double volume() const;
    
    /**
     * Returns the centroid of the simplex.
     */
    Point<spacedim> get_centroid() const;

    std::vector<Point<spacedim>> get_vertices() const;

    Point<spacedim> get_vertex(size_t i) const;

    /**
     * Returns the intersection of this simplex with another simplex, which
     * is identified by its vertices.
     */
    std::vector<Point<spacedim>> intersection(const Simplex<dim, spacedim> &s) const;
    
    bool operator==(const Simplex<dim, spacedim> &s) const;

    int n_vertices() const;
    int n_edges() const;

private:
    Point<spacedim> vertices[dim + 1];
    int num_vertices = dim + 1;
    int num_edges = (dim+1)*dim/2;
};

} // namespace mesh
} // namespace fastfem

#endif // GEOMETRY_HPP