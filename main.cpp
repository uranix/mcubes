#include <iostream>

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
    VoxelMesh vm(dim3(20, 20, 20), point(0, 0, 0), point(1, 1, 1));

    std::vector<metaball> mbs;

    for (int i = 0; i < 20; i++) {
        mbs.push_back(metaball(point(randv(), randv(), randv()), 3 * randv() - 1));
    }

    std::cout << "Setting voxel data" << std::endl;

    vm.for_each_voxel([&mbs] (VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
/*            vm[vox] = 0;
            for (const auto &m : mbs)
                vm[vox] += m(p); */
            vm[vox] = 1. / (p - point(.5, .5, .5)).norm();
        });

    bool tri = true;
    if (argc > 1)
        tri = std::string(argv[1])[0] == 't';

    std::cout << "Generating mesh" << std::endl;
    double th = 2.1;
    if (tri) {
        SurfaceMesh<triangle> sm(vm, th);
        std::cout << "Saving TRI mesh" << std::endl;
        sm.save_vtk("res_tri.vtk");
        sm.save_txt("res.tri");
    } else {
        SurfaceMesh<quad> sm(vm, th);
        std::cout << "Saving QUAD mesh" << std::endl;
        sm.save_vtk("res_quad.vtk");
        sm.save_txt("res.quad");
    }
    return 0;
}
