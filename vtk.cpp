#include <fstream>
#include <algorithm>

#include "mesh.h"

template<>
void SurfaceMesh<quad>::save_vtk(const std::string &filename) const {
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

template<>
void SurfaceMesh<triangle>::save_vtk(const std::string &filename) const {
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    f << "# vtk DataFile Version 3.0\n";
    f << "Trinagulated surface\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << pts.size() << " float\n";
    for (const auto &p : pts)
        f << p.x << " " << p.y << " " << p.z << "\n";
    f << "CELLS " << elems.size() << " " << 4 * elems.size() << "\n";
    for (const auto &t : elems)
        f << "3 " << t.p1 << " " << t.p2 << " " << t.p3 << "\n";
    f << "CELL_TYPES " << elems.size() << "\n";
    std::for_each(elems.begin(), elems.end(), [&f] (const triangle &) { f << "5\n"; });
    f << "CELL_DATA " << elems.size() << "\n";
    f << "SCALARS type int\nLOOKUP_TABLE default\n";
    for (const auto &v : facetype)
        f << v << "\n";
}

template<>
void SurfaceMesh<quad>::save_txt(const std::string &filename) const {
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    f << pts.size() << "\n";
    for (const auto &p : pts)
        f << p.x << " " << p.y << " " << p.z << "\n";
    f << elems.size() << "\n";
    for (size_t i = 0; i < elems.size(); i++) {
        const auto &q = elems[i];
        const auto &v = facetype[i];
        f << q.p1 << " " << q.p2 << " " << q.p3 << " "
            << q.p4 << " " << v << "\n";
    }
}

template<>
void SurfaceMesh<triangle>::save_txt(const std::string &filename) const {
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    f << pts.size() << "\n";
    for (const auto &p : pts)
        f << p.x << " " << p.y << " " << p.z << "\n";
    f << elems.size() << "\n";
    for (size_t i = 0; i < elems.size(); i++) {
        const auto &q = elems[i];
        const auto &v = facetype[i];
        f << q.p1 << " " << q.p2 << " " << q.p3 << " " << v << "\n";
    }
}
