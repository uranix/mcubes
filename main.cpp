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

int main(int argc, char **argv) {
    VoxelMesh vm(dim3(30, 40, 50), point(0, 0, 0), point(1, 1, 1));

    std::vector<metaball> mbs;

    for (int i = 0; i < 50; i++) {
        mbs.push_back(metaball(point(randv(), randv(), randv()), 2 * randv() - 1));
    }

    std::cout << "Setting voxel data" << std::endl;

    vm.for_each_voxel([&mbs] (VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
            vm[vox] = 0;
            for (const auto &m : mbs)
                vm[vox] += m(p);
        });

    bool tri = true;
    if (argc > 1)
        tri = std::string(argv[1])[0] == 't';

    std::cout << "Generating mesh" << std::endl;
    if (tri) {
        SurfaceMesh<triangle> sm(vm, 1.2 * mbs.size());
        std::cout << "Saving mesh" << std::endl;
        sm.save("res_tri.vtk");
    } else {
        SurfaceMesh<quad> sm(vm, mbs.size());
        std::cout << "Saving mesh" << std::endl;
        sm.save("res_quad.vtk");
    }
    return 0;
}
