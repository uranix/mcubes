#include <iostream>

#include "vtk.h"
#include "mesh.h"

struct metaball {
    const point c;
    double charge;
    metaball(const point &c, double charge = 1) : c(c), charge(charge) { }
    double operator()(const point &p) const {
        return charge / (p - c).dot(p - c);
    }
};

double randv() {
    return 1. * rand() / RAND_MAX;
}

int main() {
    std::vector<point> pts;
//    std::vector<triangle> tri;
    std::vector<quad> quads;
    VoxelMesh vm(dim3(23, 21, 18), point(0, 0, 0), point(1, 1, 1));

    std::vector<metaball> mbs;

    for (int i = 0; i < 20; i++) {
        mbs.push_back(metaball(point(randv(), randv(), randv()), 2 * randv() - 1));
    }

    vm.for_each_voxel([&mbs] (VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
            vm[vox] = 0;
            for (const auto &m : mbs)
                vm[vox] -= m(p);
        });

    SurfaceMesh<quad> sm(vm, -20);
    sm.save("test.vtk");
    return 0;
}
