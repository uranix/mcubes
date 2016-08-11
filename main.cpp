#include <iostream>

struct point;
std::ostream &operator<<(std::ostream &o, const point &p);

#include "vtk.h"
#include "mesh.h"

std::ostream &operator<<(std::ostream &o, const point &p) {
    return o << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

int main() {
    std::vector<point> pts;
//    std::vector<triangle> tri;
    std::vector<quad> quads;
    VoxelMesh vm(dim3(40, 40, 40), point(0, 0, 0), point(1, 1, 1));

    std::vector<size_t> ptsoffs(vm.val.size());

    vm.for_each_voxel([] (VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
            const point &dp1 = p - point(.7, .3, .3);
            const point &dp2 = p - point(.4, .6, .6);
            vm[vox] = 1. / dp1.dot(dp1) + 2. / dp2.dot(dp2);
        });
    vm.for_each_cube([&pts, &ptsoffs] (VoxelMesh &vm, const dim3 &cube) {
            double level = 30;
            ptsoffs[cube.linear_index(vm.n)] = pts.size();
            vm.gen_vertices(cube, level, pts);
        });
//    vm.for_each_edge([&tri, &ptsoffs] (
    vm.for_each_edge([&quads, &ptsoffs] (
                const VoxelMesh &vm,
                const std::array<char, 4> &edgedata,
                const dim3 cubes[4])
        {
            int ec = edgedata[0] & 12;
            for (int j = 1; j < 4; j++)
                assert((edgedata[j] & 12) == ec);
            ec >>= 2;
            if (ec == 0 || ec == 3) {
                for (int j = 0; j < 4; j++)
                    assert((edgedata[j] & 3) == 3);
                return;
            }
            int ptsid[4];
            for (int j = 0; j < 4; j++)
                ptsid[j] = ptsoffs[cubes[j].linear_index(vm.n)] + (edgedata[j] & 3);
            if (ec == 1) {
                // 0 1 3 2 -> 0 1 3 + 0 3 2
                quads.push_back(quad{ptsid[0], ptsid[1], ptsid[3], ptsid[2]});
//                tri.push_back(triangle{ptsid[0], ptsid[1], ptsid[3]});
//                tri.push_back(triangle{ptsid[0], ptsid[3], ptsid[2]});
            } else {
                // 0 2 3 1 -> 0 2 3 + 0 3 1
                quads.push_back(quad{ptsid[0], ptsid[2], ptsid[3], ptsid[1]});
//                tri.push_back(triangle{ptsid[0], ptsid[2], ptsid[3]});
//                tri.push_back(triangle{ptsid[0], ptsid[3], ptsid[1]});
            }
        });
//    save("test.vtk", pts, tri);
    save("test.vtk", pts, quads);
    return 0;
}
