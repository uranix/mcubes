#include "vtk.h"

#include <fstream>
#include <algorithm>

#include "mesh.h"

template<>
void SurfaceMesh<quad>::save(const std::string &filename) const {
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    f << "# vtk DataFile Version 3.0\n";
    f << "Trinagulated surface\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << pts.size() << " float\n";
    for (const auto &p : pts)
        f << p.x << " " << p.y << " " << p.z << "\n";
    f << "CELLS " << elems.size() << " " << 5 * elems.size() << "\n";
    for (const auto &q : elems)
        f << "4 " << q.p1 << " " << q.p2 << " " << q.p3 << " " << q.p4 << "\n";
    f << "CELL_TYPES " << elems.size() << "\n";
    std::for_each(elems.begin(), elems.end(), [&f] (const quad &) { f << "9\n"; });
    f << "CELL_DATA " << elems.size() << "\n";
    f << "SCALARS type int\nLOOKUP_TABLE default\n";
    for (const auto &v : facetype)
        f << v << "\n";
}

void save(const std::string &filename, const std::vector<point> &pts, const std::vector<triangle> &tri) {
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    f << "# vtk DataFile Version 3.0\n";
    f << "Trinagulated surface\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << pts.size() << " float\n";
    for (const auto &p : pts)
        f << p.x << " " << p.y << " " << p.z << "\n";
    f << "CELLS " << tri.size() << " " << 4 * tri.size() << "\n";
    for (const auto &t : tri)
        f << "3 " << t.p1 << " " << t.p2 << " " << t.p3 << "\n";
    f << "CELL_TYPES " << tri.size() << "\n";
    std::for_each(tri.begin(), tri.end(), [&f] (const triangle &) { f << "5\n"; });
}

void save(const std::string &filename, const std::vector<point> &pts, const std::vector<quad> &quads) {
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    f << "# vtk DataFile Version 3.0\n";
    f << "Trinagulated surface\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << pts.size() << " float\n";
    for (const auto &p : pts)
        f << p.x << " " << p.y << " " << p.z << "\n";
    f << "CELLS " << quads.size() << " " << 5 * quads.size() << "\n";
    for (const auto &q : quads)
        f << "4 " << q.p1 << " " << q.p2 << " " << q.p3 << " " << q.p4 << "\n";
    f << "CELL_TYPES " << quads.size() << "\n";
    std::for_each(quads.begin(), quads.end(), [&f] (const quad &) { f << "9\n"; });
}
