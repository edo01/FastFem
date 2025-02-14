#ifndef MESHMAKER_HPP
#define MESHMAKER_HPP

#include "FastFem/mesh/Mesh.hpp"

namespace fastfem {
namespace mesh {

template <unsigned int dim, unsigned int spacedim = dim>
class MeshMaker
{
public:
    MeshMaker() = default;
    virtual ~MeshMaker() = default;

    virtual Mesh<dim, spacedim> make_mesh() const = 0;
};


class CubeSurfaceMaker : public MeshMaker<2, 3>
{
public:
    /**
     * @param N Number of subdivisions per side of the cube. Tot number of divisions will be 6 * N^2
     */
    CubeSurfaceMaker(int N) : N(N) {}
    virtual ~CubeSurfaceMaker() = default;
    virtual Mesh<2, 3> make_mesh() const override;


private:
    int N;

    void build_cube_vertices(Mesh<2, 3> &mesh) const;
    void build_cube_triangles(Mesh<2, 3> &mesh) const;
};


class makeSphereSurface : public MeshMaker<2, 3>
{
public:
    makeSphereSurface(int N) : N(N) {}
    virtual ~makeSphereSurface() = default;

private:
    void sendPointsToSphere() const;

    int N;
};

} // namespace mesh
} // namespace fastfem

#endif // MESHMAKER_HPP