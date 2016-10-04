#include <iostream>

#include "vtk.h"
#include "mesh.h"

#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>

struct metaball {
    const point c;
    double charge;
    metaball(const point &c, double charge = 1) : c(c), charge(charge) { }
    double operator()(const point &p) const {
        return charge / (p - c).dot(p - c);
    }
};

template<int n>
struct simplex;

template<>
struct simplex<1> {
    static constexpr int nver = 1;
    point p[nver];
    simplex(const point &, const double scale = 0) {
        p[0] = point(0, 0, 0);
    }
    double operator()(const double v[nver], const point &pr) const {
        return v[0];
    }
};

template<>
struct simplex<4> {
    static constexpr int nver = 4;
    point p[nver];
    simplex(const point &h, const double scale = .7) {
        p[0] = point( h.x/2, -h.y/2, -h.z/2) * scale;
        p[1] = point(-h.x/2,  h.y/2, -h.z/2) * scale;
        p[2] = point(-h.x/2, -h.y/2,  h.z/2) * scale;
        p[3] = point( h.x/2,  h.y/2,  h.z/2) * scale;
    }
    double operator()(const double v[nver], const point &pr) const {
        double vc = 0;
        double vp = 0;
        double x = pr.x;
        double y = pr.y;
        double z = pr.z;
        for (int i = 0; i < nver; i++) {
            vc += v[i];
            vp += v[i] * (x / p[i].x + y / p[i].y + z / p[i].z);
        }
        return 0.25 * (vc + vp);
    }
};

template<>
struct simplex<10> {
    static constexpr int nver = 10;
    point p[nver];
    Eigen::FullPivLU<Eigen::Matrix<double, 10, 10>> LU;
    simplex(const point &h, const double scale = .7) {
        p[0] = point( h.x/2, -h.y/2, -h.z/2) * scale;
        p[1] = point(-h.x/2,  h.y/2, -h.z/2) * scale;
        p[2] = point(-h.x/2, -h.y/2,  h.z/2) * scale;
        p[3] = point( h.x/2,  h.y/2,  h.z/2) * scale;
        int k = 4;
        for (int i = 0; i < 4; i++)
            for (int j = i+1; j < 4; j++) {
                p[k] = (p[i] + p[j]) * 0.5;
                k++;
            }
        Eigen::Matrix<double, 10, 10> A;
        for (int i = 0; i < 10; i++) {
            double x = p[i].x;
            double y = p[i].y;
            double z = p[i].z;
            A(i, 0) = 1;
            A(i, 1) = x;
            A(i, 2) = y;
            A(i, 3) = z;
            A(i, 4) = x*x;
            A(i, 5) = x*y;
            A(i, 6) = x*z;
            A(i, 7) = y*y;
            A(i, 8) = y*z;
            A(i, 9) = z*z;
        }
        LU = A.fullPivLu();
    }
    double operator()(const double v[nver], const point &pr) const {
        Eigen::Matrix<double, 10, 1> f(v);
        Eigen::Matrix<double, 10, 1> s = LU.solve(f);

        double x = pr.x;
        double y = pr.y;
        double z = pr.z;

        return s[0] + s[1] * x + s[2] * y + s[3] * z +
            s[4] * x*x + s[5] * x*y + s[6] * x*z +
            s[7] * y*y + s[8] * y*z + s[9] * z*z;
    }
};

double randv() {
    return 1. * rand() / RAND_MAX;
}

double solid(const point &a, const point &b, const point &c) {
    const double an = a.norm();
    const double bn = b.norm();
    const double cn = c.norm();
    double num = a.cross(b).dot(c);
    double denom = an * bn * cn + a.dot(b) * cn + a.dot(c) * bn + b.dot(c) * an;
    return 2 * std::atan2(num, denom);
}

double density(const point &p) {
    return 2 + sin(3 * p.x) * cos(5 * p.y * p.z) + 2 * p.y - p.x*p.x;
}

int main(int argc, char **argv) {
    VoxelMesh vm(dim3(30, 30, 30), point(0, 0, 0), point(1, 1, 1));

    std::cout << "Setting voxel data" << std::endl;

    vm.for_each_voxel([] (VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
            vm[vox] = -(p - point(.5, .5, .5)).norm();
        });

    std::cout << "Generating mesh" << std::endl;

    double th = -0.3;

    SurfaceMesh<triangle> sm(vm, th);
    std::cout << "Saving tri mesh" << std::endl;
    sm.save("res_tri.vtk");

    struct data {
        double exact;
        double interp;
        double v[10];
        dim3 closest;
        int distance;
    };

    VoxelData<data> vd(vm.n);

    simplex<10> simpl(vm.h);

    vm.for_each_voxel([&vd, &sm, &simpl] (const VoxelMesh &vm, const dim3 &vox) {
            double val = 0;
            double v[10] = {0, 0, 0, 0};
            vd[vox].distance = -1;
            if (vm.mark(vox) & 1) {
                const std::vector<point> &pts = sm.points();
                const std::vector<triangle> &tri = sm.elements();
                const point &rv = vm.center(vox);
                for (const auto &t : tri) {
                    const point &r1 = pts[t.p1];
                    const point &r2 = pts[t.p2];
                    const point &r3 = pts[t.p3];
                    const point &a = r1 - rv;
                    const point &b = r2 - rv;
                    const point &c = r3 - rv;
                    val += solid(a, b, c) * density((r1 + r2 + r3) * (1. / 3));
                    if (vm.mark(vox) == 1) {
                        for (int j = 0; j < simpl.nver; j++) {
                            const point &a = r1 - rv - simpl.p[j];
                            const point &b = r2 - rv - simpl.p[j];
                            const point &c = r3 - rv - simpl.p[j];
                            v[j] += solid(a, b, c) * density((r1 + r2 + r3) * (1. / 3));
                        }
                    }
                }
                if (vm.mark(vox) == 3) {
                    dim3 neib;
                    int range = 0;
                    int maxdist = vm.n.i + vm.n.j + vm.n.k;
                    int dist = maxdist;
                    while (dist == maxdist) {
                        range++;
                        dim3 nb;
                        for (nb.i = vox.i - range; nb.i <= vox.i + range; nb.i++)
                            for (nb.j = vox.j - range; nb.j <= vox.j + range; nb.j++)
                                for (nb.k = vox.k - range; nb.k <= vox.k + range; nb.k++) {
                                    if (nb.i < 0 || nb.j < 0 || nb.k < 0)
                                        continue;
                                    if (nb.i >= vm.n.i || nb.j >= vm.n.j || nb.k >= vm.n.k)
                                        continue;
                                    if (vm.mark(nb) == 1) {
                                        if (vox.dist(nb) < dist) {
                                            neib = nb;
                                            dist = vox.dist(nb);
                                        }
                                    }
                                }
                    }
                    vd[vox].closest = neib;
                    vd[vox].distance = dist;
                }
            }
            vd[vox].exact = val;
            for (int j = 0; j < simpl.nver; j++)
                vd[vox].v[j] = v[j];
        });

    std::vector<int> hist(10, 0);
    vm.for_each_voxel([&vd, &hist] (const VoxelMesh &vm, const dim3 &vox) {
            int d = vd[vox].distance;
            if (d + 1 > hist.size())
                hist[hist.size()-1]++;
            else if (d >= 0)
                hist[d]++;
        });
/*
    std::cout << "Distance hist: \n";
    for (int k = 0; k < hist.size(); k++)
        std::cout << k << " " << hist[k] << std::endl;
*/
    double sumdiff = 0, sumsq = 0;
    int num = 0;
    double minv = 100;
    double maxv = 0;
    vm.for_each_voxel([&simpl, &vd, &num, &sumdiff, &sumsq, &minv, &maxv] (const VoxelMesh &vm, const dim3 &vox) {
            vd[vox].interp = 0;
            if (vm.mark(vox) == 1) {
                vd[vox].interp = simpl(vd[vox].v, point(0, 0, 0));
                num++;
                double err = 100 * std::abs((vd[vox].exact - vd[vox].interp) / vd[vox].exact);
                if (minv > err)
                    minv = err;
                if (maxv < err)
                    maxv = err;
                sumdiff += err;
                sumsq += err * err;
            }
        });

    double mean = sumdiff / num;
    double stddev = std::sqrt(sumsq / num - mean * mean);

    std::cout << "Interpolation error: mean: " << mean << "%, std dev: " << stddev << "%, min: " <<
        minv << "%, max: " << maxv << "%" << std::endl;

    num = 0;
    sumdiff = sumsq = 0;
    vm.for_each_voxel([&vd, &num, &sumdiff, &sumsq, &simpl, &minv, &maxv] (const VoxelMesh &vm, const dim3 &vox) {
            if (vm.mark(vox) == 3) {
                const auto &p = vm.center(vox);
                const auto &nb = vd[vox].closest;
                const auto &pc = vm.center(nb);
                vd[vox].interp = simpl(vd[nb].v, p - pc);
                num++;
                double err = 100 * std::abs((vd[vox].exact - vd[vox].interp) / vd[vox].exact);
                if (minv > err)
                    minv = err;
                if (maxv < err)
                    maxv = err;
                sumdiff += err;
                sumsq += err * err;
            }
        });
    mean = sumdiff / num;
    stddev = std::sqrt(sumsq / num - mean * mean);

    std::cout << "Extrapolation error: mean: " << mean << "%, std dev: " << stddev << "%, min: " <<
        minv << "%, max: " << maxv << "%" << std::endl;

    std::ofstream fp("voxels.vtk");
    fp << "# vtk DataFile Version 3.0\n";
    fp << "Trinagulated surface\n";
    fp << "ASCII\n";
    fp << "DATASET STRUCTURED_POINTS\n";
    fp << "DIMENSIONS " << vm.n.i+1 << " " << vm.n.j+1 << " " << vm.n.k+1 << "\n";
    fp << "ORIGIN " << vm.ll.x << " " << vm.ll.y << " " << vm.ll.z << "\n";
    fp << "SPACING " << vm.h.x << " " << vm.h.y << " " << vm.h.z << "\n";
    /*
    vm.for_each_voxel([&fp] (const VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
            fp << p.x << " " << p.y << " " << p.z << "\n";
        });
    fp << "CELLS 0 0\n";
    fp << "CELL_TYPES 0\n"; */
    fp << "CELL_DATA " << vm.n.prod() << "\n";
    fp << "SCALARS type int\nLOOKUP_TABLE default\n";
    vm.for_each_voxel([&fp] (const VoxelMesh &vm, const dim3 &vox) {
            fp << vm.mark(vox) << "\n";
        });
    fp << "SCALARS ex double\nLOOKUP_TABLE default\n";
    vm.for_each_voxel([&fp, &vd] (const VoxelMesh &vm, const dim3 &vox) {
            fp << vd[vox].exact << "\n";
        });
    fp << "SCALARS in double\nLOOKUP_TABLE default\n";
    vm.for_each_voxel([&fp, &vd] (const VoxelMesh &vm, const dim3 &vox) {
            fp << vd[vox].interp << "\n";
        });
    fp << "SCALARS diff double\nLOOKUP_TABLE default\n";
    vm.for_each_voxel([&fp, &vd] (const VoxelMesh &vm, const dim3 &vox) {
            fp << std::abs(vd[vox].exact - vd[vox].interp) / (std::abs(vd[vox].exact) + 1e-300) << "\n";
        });
    fp << "SCALARS dist int\nLOOKUP_TABLE default\n";
    vm.for_each_voxel([&fp, &vd] (const VoxelMesh &vm, const dim3 &vox) {
            fp << vd[vox].distance << "\n";
        });
    fp << "VECTORS ijk int\n";
    vm.for_each_voxel([&fp, &vd] (const VoxelMesh &vm, const dim3 &vox) {
            fp << vox.i << " " << vox.j << " " << vox.k << "\n";
        });

    std::ofstream ftr("train.txt");
    vm.for_each_voxel([&ftr, &vd] (const VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
            if (vm.mark(vox) != 1)
                return;
            ftr << p.x << " " << p.y << " " << p.z << " " << vd[vox].exact << "\n";
        });
    std::ofstream fte("test.txt");
    vm.for_each_voxel([&fte, &vd] (const VoxelMesh &vm, const dim3 &vox) {
            const point &p = vm.center(vox);
            if (vm.mark(vox) != 3)
                return;
            fte << p.x << " " << p.y << " " << p.z << " " << vd[vox].exact << " " << vd[vox].interp << "\n";
        });

    return 0;
}
