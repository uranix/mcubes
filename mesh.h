#pragma once

#include <vector>
#include <array>

#include "geom.h"
#include "lut.h"

struct VoxelMesh;

struct EdgeMesh {
    const dim3 n;
    std::vector<std::array<char, 4>> _xedge, _yedge, _zedge;
    EdgeMesh(const dim3 &n)
        : n(n),
        _xedge((n.i - 1) * n.j * n.k),
        _yedge(n.i * (n.j - 1) * n.k),
        _zedge(n.i * n.j * (n.k - 1))
    {
    }
    std::array<char,4> &xedge(const dim3 &idx) {
        return _xedge[idx.linear_index(n - dim3::X())];
    }
    std::array<char,4> &yedge(const dim3 &idx) {
        return _yedge[idx.linear_index(n - dim3::Y())];
    }
    std::array<char,4> &zedge(const dim3 &idx) {
        return _zedge[idx.linear_index(n - dim3::Z())];
    }
    const std::array<char,4> &xedge(const dim3 &idx) const {
        return _xedge[idx.linear_index(n - dim3::X())];
    }
    const std::array<char,4> &yedge(const dim3 &idx) const {
        return _yedge[idx.linear_index(n - dim3::Y())];
    }
    const std::array<char,4> &zedge(const dim3 &idx) const {
        return _zedge[idx.linear_index(n - dim3::Z())];
    }
    void update_edge(const dim3 &cube, const edge &e, int edgeid, bool cb, bool ce, int patchid) {
        int dir = edgeid >> 2;
        int ab = edgeid & 3;
        int ccb = cb ? 1 : 0;
        int cce = ce ? 1 : 0;
        int sig = (ccb << 3) | (cce << 2) | patchid;
        if (dir == 0) {
            std::array<char,4> &edgeinfo = xedge(cube + e.beg);
            edgeinfo[ab] = sig;
        } else if (dir == 1) {
            std::array<char,4> &edgeinfo = yedge(cube + e.beg);
            edgeinfo[ab] = sig;
        } else {
            std::array<char,4> &edgeinfo = zedge(cube + e.beg);
            edgeinfo[ab] = sig;
        }
    }
};

struct VoxelMesh {
    const dim3 n;
    std::vector<double> val;
    const point ll, ur, h;
    VoxelMesh(const dim3 &n, const point &ll, const point &ur)
        : n(n),
        val(n.i * n.j * n.k),
        ll(ll), ur(ur),
        h((ur.x - ll.x) / n.i, (ur.y - ll.y) / n.j, (ur.z - ll.z) / n.k)
    {
    }
    double &operator[](const dim3 &idx) {
        return val[idx.linear_index(n)];
    }
    const double &operator[](const dim3 &idx) const {
        return val[idx.linear_index(n)];
    }
    const point center(const dim3 &idx) const {
        return point(idx, ll, h);
    }
    const point edge_point(const dim3 &cube, const edge &e, double w) const {
        point p1(cube + e.beg, ll, h);
        point p2(cube + e.end, ll, h);
        return p1 * (1 - w) + p2 * w;
    }
    int gen_vertices(const dim3 &cube, double level, std::vector<point> &pts, EdgeMesh &em) const {
        double v[8];
        point newpts[4];
        int cnt[4] = {0, 0, 0, 0};
        for (int i = 0; i < 4; i++)
            newpts[i] = point(0, 0, 0);
        int cases = 0;
        for (int i = 0; i < 8; i++) {
            double val = (*this)[cube + dim3::cube_vertex(i)];
            v[i] = val;
            if (val > level)
                cases |= (1 << i);
        }
        int code = edgeGroup[cases];
        for (int i = 0; i < 12; i++) {
            edge e(i);
            int patchid = (code >> (2 * i)) & 3;
            double v1 = v[e.beg.vertex_id()];
            double v2 = v[e.end.vertex_id()];
            bool c1 = v1 > level;
            bool c2 = v2 > level;
            em.update_edge(cube, e, i, c1, c2, patchid);
            if (c1 == c2) {
                assert(patchid == 3);
                continue;
            }
            double w = (level - v1) / (v2 - v1);
            cnt[patchid]++;
            newpts[patchid] += edge_point(cube, e, w);
        }
        int patches = 0;
        for (int i = 0; i < 4; i++)
            if (cnt[i] > 0) {
                pts.push_back(newpts[i] * (1. / cnt[i]));
                patches++;
            }
        for (int i = patches; i < 4; i++)
            assert(cnt[i] == 0);
        return patches;
    }
    template<typename VoxelFunctor>
    void for_each_voxel(VoxelFunctor &&f) {
        dim3 vox;
        for (vox.k = 0; vox.k < n.k; vox.k++)
            for (vox.j = 0; vox.j < n.j; vox.j++)
                for (vox.i = 0; vox.i < n.i; vox.i++)
                    f(*this, vox);
    }
    template<typename CubeFunctor>
    void for_each_cube(CubeFunctor &&f) {
        dim3 vox;
        for (vox.k = 0; vox.k < n.k - 1; vox.k++)
            for (vox.j = 0; vox.j < n.j - 1; vox.j++)
                for (vox.i = 0; vox.i < n.i - 1; vox.i++)
                    f(*this, vox);
    }
    template<typename EdgeFunctor>
    void for_each_edge(const EdgeMesh &em, EdgeFunctor &&f) {
        dim3 vox;
        for (vox.k = 1; vox.k < n.k - 1; vox.k++)
            for (vox.j = 1; vox.j < n.j - 1; vox.j++)
                for (vox.i = 0; vox.i < n.i - 1; vox.i++) {
                    dim3 cubes[4] = {vox, vox - dim3::Y(), vox - dim3::Z(), vox - dim3::Y() - dim3::Z()};
                    f(*this, em.xedge(vox), cubes);
                }
        for (vox.k = 1; vox.k < n.k - 1; vox.k++)
            for (vox.j = 0; vox.j < n.j - 1; vox.j++)
                for (vox.i = 1; vox.i < n.i - 1; vox.i++) {
                    dim3 cubes[4] = {vox, vox - dim3::X(), vox - dim3::Z(), vox - dim3::X() - dim3::Z()};
                    f(*this, em.yedge(vox), cubes);
                }
        for (vox.k = 0; vox.k < n.k - 1; vox.k++)
            for (vox.j = 1; vox.j < n.j - 1; vox.j++)
                for (vox.i = 1; vox.i < n.i - 1; vox.i++) {
                    dim3 cubes[4] = {vox, vox - dim3::X(), vox - dim3::Y(), vox - dim3::X() - dim3::Y()};
                    f(*this, em.zedge(vox), cubes);
                }
    }
};

template<typename Element>
struct SurfaceMesh {
    std::vector<point> pts;
    std::vector<Element> elems;
    EdgeMesh em;
    SurfaceMesh(const VoxelMesh &vm, const double level) : em(vm.n) {
        std::vector<size_t> ptsoffs(vm.val.size());
        vm.for_each_cube([this, &ptsoffs, level]
            (const VoxelMesh &vm, const dim3 &cube) {
                ptsoffs[cube.linear_index(vm.n)] = pts.size();
                vm.gen_vertices(cube, level, pts, em);
            });
    }
};
