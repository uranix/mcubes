#pragma once

#include <vector>
#include <array>

#include "geom.h"

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
    const point ll, ur, h;
private:
    std::vector<double> val;
    std::vector<int> _mark;
    const point edge_point(const dim3 &cube, const edge &e, double w) const;
public:
    VoxelMesh(const dim3 &n, const point &ll, const point &ur);
    double &operator[](const dim3 &idx);
    const double &operator[](const dim3 &idx) const;
    int &mark(const dim3 &idx);
    const int &mark(const dim3 &idx) const;
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

template<typename T>
class VoxelData {
public:
    const dim3 n;
private:
    std::vector<T> val;
public:
    VoxelData(const dim3 &n) : n(n), val(n.prod()) { }
    T &operator[](const dim3 &idx) {
        return val[idx.linear_index(n)];
    }
    const T &operator[](const dim3 &idx) const {
        return val[idx.linear_index(n)];
    }
};

template<typename Element>
class SurfaceMesh {
    std::vector<point> pts;
    std::vector<Element> elems;
    std::vector<int> facetype;

    EdgeMesh em;
public:
    const std::vector<point> &points() const { return pts; }
    const std::vector<Element> &elements() const { return elems; }
    SurfaceMesh(VoxelMesh &vm, const double level) : em(vm.n) {
        std::vector<int> ptsoffs((vm.n.i + 1) * (vm.n.j + 1) * (vm.n.k + 1));
        vm.for_each_cube([this, &ptsoffs, level]
            (const VoxelMesh &vm, const dim3 &cube, unsigned int cubeflag) {
                int cid = (cube + dim3(1, 1, 1)).linear_index(vm.n + dim3(1, 1, 1));
                ptsoffs[cid] = pts.size();
                vm.gen_vertices(cube, cubeflag, level, pts, em);
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
                for (int j = 0; j < 4; j++) {
                    int cid = (cubes[j] + dim3(1, 1, 1)).linear_index(vm.n + dim3(1, 1, 1));
                    ptsid[j] = ptsoffs[cid] + (edgedata[j] & 3);
                }
                if ((ec == 1) ^ (dir.j == 1))
                    this->add_face(ptsid[0], ptsid[1], ptsid[3], ptsid[2], sig, dir);
                else
                    this->add_face(ptsid[0], ptsid[2], ptsid[3], ptsid[1], sig, dir);
            });
        const auto &cem = em;
        vm.for_each_voxel([&cem, level] (VoxelMesh &vm, const dim3 &vox) {
                auto lf = cem.xedge(vox - dim3(1, 0, 0))[0] >> 2;
                auto rt = cem.xedge(vox)[0] >> 2;
                auto dn = cem.yedge(vox - dim3(0, 1, 0))[0] >> 2;
                auto up = cem.yedge(vox)[0] >> 2;
                auto fa = cem.zedge(vox - dim3(0, 0, 1))[0] >> 2;
                auto ne = cem.zedge(vox)[0] >> 2;
                int inout = (lf & 1) ? 1 : 0;
                assert(inout == ((dn & 1) ? 1 : 0));
                assert(inout == ((fa & 1) ? 1 : 0));
                assert(inout == ((rt & 2) ? 1 : 0));
                assert(inout == ((up & 2) ? 1 : 0));
                assert(inout == ((ne & 2) ? 1 : 0));
                int boundary = 0;
                if (lf == 2 || lf == 1) boundary = 1;
                if (rt == 2 || rt == 1) boundary = 1;
                if (dn == 2 || dn == 1) boundary = 1;
                if (up == 2 || up == 1) boundary = 1;
                if (fa == 2 || fa == 1) boundary = 1;
                if (ne == 2 || ne == 1) boundary = 1;
                vm.mark(vox) = (boundary << 1) | inout;
            });
    }
    void add_face(int p1, int p2, int p3, int p4, int sig, const dim3 &dir);
    void save(const std::string &filename) const;
};
