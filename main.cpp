#include "vtk.h"
#include "geom.h"

#include <unordered_map>

struct SurfaceMesh {
//    std::unordered_map<edge, double> pts;
    std::vector<triedge> tri;
};

struct VoxelMesh {
    int nx, ny, nz;
    std::vector<double> val;
    const point ll, ur, h;
    VoxelMesh(int nx, int ny, int nz, const point &ll, const point &ur)
        : nx(nx), ny(ny), nz(nz), val(nx * ny * nz), ll(ll), ur(ur),
        h((ur.x - ll.x) / nx, (ur.y - ll.y) / ny, (ur.z - ll.z) / nz)
    {
    }
    double &operator()(int i, int j, int k) {
        return val[i + nx * (j + ny * k)];
    }
    const double &operator()(int i, int j, int k) const {
        return val[i + nx * (j + ny * k)];
    }
    const point center(int i, int j, int k) const {
        point p(ll);
        p.x += (i + 0.5) * h.x;
        p.y += (j + 0.5) * h.y;
        p.z += (k + 0.5) * h.z;
        return p;
    }
    const point edge_point(const edge &e, double w) const {
        point p1(ll);
        p1 += point(e.i * h.x, e.j * h.y, e.k * h.z);
        point p2(p1);
        if (e.dir == 0) p2.x += h.x;
        if (e.dir == 1) p2.y += h.y;
        if (e.dir == 2) p2.z += h.z;
        return p1 * (1 - w) + p2 * (w);
    }
};

int main() {
    std::vector<point> pts;
    std::vector<triangle> tri;
    VoxelMesh vm(10, 10, 10, point(0, 0, 0), point(1, 1, 1));
    save("test.vtk", pts, tri);
    return 0;
}
