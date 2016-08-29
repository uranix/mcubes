#include "mesh.h"

#include <iostream>

EdgeMesh::EdgeMesh(const dim3 &n)
    : n(n),
    _xedge((n.i + 1) * (n.j + 2) * (n.k + 2)),
    _yedge((n.i + 2) * (n.j + 1) * (n.k + 2)),
    _zedge((n.i + 2) * (n.j + 2) * (n.k + 1))
{
}

EdgeMesh::mark &EdgeMesh::xedge(const dim3 &_idx) {
    const dim3 &idx(_idx + dim3(1, 1, 1));
    return _xedge[idx.linear_index(n + dim3(1, 2, 2))];
}
EdgeMesh::mark &EdgeMesh::yedge(const dim3 &_idx) {
    const dim3 &idx(_idx + dim3(1, 1, 1));
    return _yedge[idx.linear_index(n + dim3(2, 1, 2))];
}
EdgeMesh::mark &EdgeMesh::zedge(const dim3 &_idx) {
    const dim3 &idx(_idx + dim3(1, 1, 1));
    return _zedge[idx.linear_index(n + dim3(2, 2, 1))];
}
const EdgeMesh::mark &EdgeMesh::xedge(const dim3 &_idx) const {
    const dim3 &idx(_idx + dim3(1, 1, 1));
    return _xedge[idx.linear_index(n + dim3(1, 2, 2))];
}
const EdgeMesh::mark &EdgeMesh::yedge(const dim3 &_idx) const {
    const dim3 &idx(_idx + dim3(1, 1, 1));
    return _yedge[idx.linear_index(n + dim3(2, 1, 2))];
}
const EdgeMesh::mark &EdgeMesh::zedge(const dim3 &_idx) const {
    const dim3 &idx(_idx + dim3(1, 1, 1));
    return _zedge[idx.linear_index(n + dim3(2, 2, 1))];
}

void EdgeMesh::update_edge(const dim3 &cube, const edge &e, int edgeid,
        bool cb, bool ce, int patchid)
{
    int dir = edgeid >> 2;
    int ab = edgeid & 3;
    int ccb = cb ? 1 : 0;
    int cce = ce ? 1 : 0;
    int sig = (ccb << 3) | (cce << 2) | patchid;
    if (dir == 0)
        xedge(cube + e.beg)[ab] = sig;
    else if (dir == 1)
        yedge(cube + e.beg)[ab] = sig;
    else
        zedge(cube + e.beg)[ab] = sig;
}

void VoxelMesh::for_each_edge(const EdgeMesh &em, EdgeFunctor f) const {
    dim3 vox;
    for (vox.k = 0; vox.k < n.k; vox.k++)
        for (vox.j = 0; vox.j < n.j; vox.j++)
            for (vox.i = -1; vox.i < n.i; vox.i++) {
                int sig = 0;
                if (vox.i == -1)
                    sig = 1;
                if (vox.i == n.i - 1)
                    sig = 2;
                dim3 cubes[4] = {vox, vox - dim3::Y(), vox - dim3::Z(), vox - dim3::Y() - dim3::Z()};
                f(*this, em.xedge(vox), cubes, sig, dim3::X());
            }
    for (vox.k = 0; vox.k < n.k; vox.k++)
        for (vox.j = -1; vox.j < n.j; vox.j++)
            for (vox.i = 0; vox.i < n.i; vox.i++) {
                int sig = 0;
                if (vox.j == -1)
                    sig = 3;
                if (vox.j == n.j - 1)
                    sig = 4;
                dim3 cubes[4] = {vox, vox - dim3::X(), vox - dim3::Z(), vox - dim3::X() - dim3::Z()};
                f(*this, em.yedge(vox), cubes, sig, dim3::Y());
            }
    for (vox.k = -1; vox.k < n.k; vox.k++)
        for (vox.j = 0; vox.j < n.j; vox.j++)
            for (vox.i = 0; vox.i < n.i; vox.i++) {
                int sig = 0;
                if (vox.k == -1)
                    sig = 5;
                if (vox.k == n.k - 1)
                    sig = 6;
                dim3 cubes[4] = {vox, vox - dim3::X(), vox - dim3::Y(), vox - dim3::X() - dim3::Y()};
                f(*this, em.zedge(vox), cubes, sig, dim3::Z());
            }
}

VoxelMesh::VoxelMesh(const dim3 &n, const point &ll, const point &ur)
    : n(n),
    val(n.i * n.j * n.k),
    ll(ll), ur(ur),
    h((ur.x - ll.x) / n.i, (ur.y - ll.y) / n.j, (ur.z - ll.z) / n.k)
{
}

double &VoxelMesh::operator[](const dim3 &idx) {
    return val[idx.linear_index(n)];
}
const double &VoxelMesh::operator[](const dim3 &idx) const {
    return val[idx.linear_index(n)];
}
const point VoxelMesh::center(const dim3 &idx) const {
    return point(idx, ll, h);
}
const point VoxelMesh::edge_point(const dim3 &cube, const edge &e, double w) const {
    point p1(cube + e.beg, ll, h);
    point p2(cube + e.end, ll, h);
    return p1 * (1 - w) + p2 * w;
}
int VoxelMesh::gen_vertices(const dim3 &cube, unsigned int cubeflag, double level, std::vector<point> &pts, EdgeMesh &em) const {
    double v[8];
    point newpts[4];
    int cnt[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++)
        newpts[i] = point(0, 0, 0);
    int cases = 0;
//        std::cout << "cube, (" << cube.i << ", " << cube.j << ", " << cube.k << ") -> " << cubeflag << ", case ";
    for (int i = 0; i < 8; i++) {
        dim3 cv = dim3::cube_vertex(i);
        double val = level;
        if (!(  (cv.i == 0 && (cubeflag & (1 << 0))) ||
                (cv.i == 1 && (cubeflag & (1 << 1))) ||
                (cv.j == 0 && (cubeflag & (1 << 2))) ||
                (cv.j == 1 && (cubeflag & (1 << 3))) ||
                (cv.k == 0 && (cubeflag & (1 << 4))) ||
                (cv.k == 1 && (cubeflag & (1 << 5)))))
            val = (*this)[cube + cv];
        v[i] = val;
        if (val > level)
            cases |= (1 << i);
    }
//        std::cout << cases << std::endl;
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
            point p = newpts[i] * (1. / cnt[i]);
//                std::cout << "yield point " << p.x << ", " << p.y << ", " << p.z << std::endl;

            if (cubeflag & (1 << 0)) p.x = ll.x;
            if (cubeflag & (1 << 1)) p.x = ur.x;
            if (cubeflag & (1 << 2)) p.y = ll.y;
            if (cubeflag & (1 << 3)) p.y = ur.y;
            if (cubeflag & (1 << 4)) p.z = ll.z;
            if (cubeflag & (1 << 5)) p.z = ur.z;

            pts.push_back(p);
            patches++;
        }
    for (int i = patches; i < 4; i++)
        assert(cnt[i] == 0);
    return patches;
}

void VoxelMesh::for_each_voxel(VoxelFunctor f) {
    dim3 vox;
    for (vox.k = 0; vox.k < n.k; vox.k++)
        for (vox.j = 0; vox.j < n.j; vox.j++)
            for (vox.i = 0; vox.i < n.i; vox.i++)
                f(*this, vox);
}
void VoxelMesh::for_each_cube(CubeFunctor f) const {
    dim3 cube;
    for (cube.k = -1; cube.k < n.k; cube.k++)
        for (cube.j = -1; cube.j < n.j; cube.j++)
            for (cube.i = -1; cube.i < n.i; cube.i++) {
                unsigned int cubeflag = 0;
                if (cube.i == -1)
                    cubeflag |= (1 << 0);
                if (cube.i == n.i - 1)
                    cubeflag |= (1 << 1);
                if (cube.j == -1)
                    cubeflag |= (1 << 2);
                if (cube.j == n.j - 1)
                    cubeflag |= (1 << 3);
                if (cube.k == -1)
                    cubeflag |= (1 << 4);
                if (cube.k == n.k - 1)
                    cubeflag |= (1 << 5);
                f(*this, cube, cubeflag);
            }
}

template<>
void SurfaceMesh<quad>::add_face(int p1, int p2, int p3, int p4, int sig, const dim3 &) {
    elems.push_back(quad{p1, p2, p3, p4});
    facetype.push_back(sig);
}

double quality(const point &p1, const point &p2, const point &p3) {
    double a = (p2 - p3).norm();
    double b = (p1 - p3).norm();
    double c = (p1 - p2).norm();
    double r = 0.5 * std::sqrt((b + c - a) * (c + a - b) * (a + b - c) / (a + b + c));
    double s = 0.5 * (a + b + c);
    double R = a * b * c / (4 * r * s);
    return r / R;
}

template<>
void SurfaceMesh<triangle>::add_face(int p1, int p2, int p3, int p4, int sig, const dim3 &dir) {
    point proj[4];
    proj[0] = pts[p1];
    proj[1] = pts[p2];
    proj[2] = pts[p3];
    proj[3] = pts[p4];
    for (int i = 0; i < 4; i++) {
        proj[i].x *= (1 - dir.i);
        proj[i].y *= (1 - dir.j);
        proj[i].z *= (1 - dir.k);
    }

    double q11 = quality(proj[0], proj[1], proj[2]);
    double q12 = quality(proj[0], proj[2], proj[3]);

    double q21 = quality(proj[0], proj[1], proj[3]);
    double q22 = quality(proj[1], proj[2], proj[3]);

    if (std::min(q11, q12) > std::min(q21, q22)) {
        elems.push_back(triangle{p1, p2, p3});
        elems.push_back(triangle{p1, p3, p4});
    } else {
        elems.push_back(triangle{p1, p2, p4});
        elems.push_back(triangle{p2, p3, p4});
    }

    facetype.push_back(sig);
    facetype.push_back(sig);
}
