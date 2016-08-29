#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <iostream>

#include "geom.h"
#include "lut.h"

class EdgeMesh {
public:
    using mark = std::array<char, 4>;
private:
    const dim3 n;
    std::vector<mark> _xedge, _yedge, _zedge;
    mark &xedge(const dim3 &_idx);
    mark &yedge(const dim3 &_idx);
    mark &zedge(const dim3 &_idx);
public:
    EdgeMesh(const dim3 &n);
    void update_edge(const dim3 &cube, const edge &e, int edgeid, bool cb, bool ce, int patchid);
    const mark &xedge(const dim3 &_idx) const;
    const mark &yedge(const dim3 &_idx) const;
    const mark &zedge(const dim3 &_idx) const;
};

class VoxelMesh {
public:
    const dim3 n;
private:
    std::vector<double> val;
    const point ll, ur, h;
    const point edge_point(const dim3 &cube, const edge &e, double w) const;
public:
    VoxelMesh(const dim3 &n, const point &ll, const point &ur);
    double &operator[](const dim3 &idx);
    const double &operator[](const dim3 &idx) const;
    const point center(const dim3 &idx) const;

    int gen_vertices(const dim3 &cube, unsigned int cubeflag, double level, std::vector<point> &pts, EdgeMesh &em) const;

    using VoxelFunctor = std::function<void(VoxelMesh &, const dim3 &)>;
    using CubeFunctor = std::function<void(const VoxelMesh &, const dim3 &, int)>;
    using EdgeFunctor = std::function<void(const VoxelMesh &,
            const EdgeMesh::mark &, const dim3 [4], int, const dim3 &)>;

    void for_each_voxel(VoxelFunctor f);
    void for_each_cube(CubeFunctor f) const;
    void for_each_edge(const EdgeMesh &em, EdgeFunctor f) const;
};

template<typename Element>
class SurfaceMesh {
    std::vector<point> pts;
    std::vector<Element> elems;
    std::vector<int> facetype;

    EdgeMesh em;
public:
    SurfaceMesh(const VoxelMesh &vm, const double level) : em(vm.n) {
        std::unordered_map<dim3, size_t, dim3::hash> ptsoffs;
        vm.for_each_cube([this, &ptsoffs, level]
            (const VoxelMesh &vm, const dim3 &cube, unsigned int cubeflag) {
                ptsoffs[cube] = pts.size();
                int added = vm.gen_vertices(cube, cubeflag, level, pts, em);
                if (added > 1)
                    std::cout << "Patches: " << added << std::endl;
            });
        vm.for_each_edge(em, [this, &ptsoffs] (
                    const VoxelMesh &vm,
                    const std::array<char, 4> &edgedata,
                    const dim3 cubes[4], int sig,
                    const dim3 &dir)
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
                    ptsid[j] = ptsoffs[cubes[j]] + (edgedata[j] & 3);
                if (ec == 1)
                    this->add_face(ptsid[0], ptsid[1], ptsid[3], ptsid[2], sig, dir);
                else
                    this->add_face(ptsid[0], ptsid[2], ptsid[3], ptsid[1], sig, dir);
            });
    }
    void add_face(int p1, int p2, int p3, int p4, int sig, const dim3 &dir);
    void save(const std::string &filename) const;
};
