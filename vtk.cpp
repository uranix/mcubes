#include "vtk.h"

#include <fstream>
#include <algorithm>

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
