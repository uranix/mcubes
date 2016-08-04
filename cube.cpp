#include "geom.h"

struct tet {
    int a, b, c, d;
};

int cubefaces[6][4] = {
    {0, 2, 6, 4}, {1, 5, 7, 3},
    {0, 4, 5, 1}, {2, 3, 7, 6},
    {0, 1, 3, 2}, {4, 6, 7, 5}
};


#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

void save(const std::string &filename, const std::vector<point> &pts, const std::vector<tet> &tets, const std::vector<double> &vals) {
    std::ofstream f(filename, std::ios::binary);
    f << "# vtk DataFile Version 3.0\n";
    f << "Cube\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << pts.size() << " float\n";
    for (const auto &p : pts)
        f << p.x << " " << p.y << " " << p.z << "\n";
    f << "CELLS " << tets.size() << " " << 5 * tets.size() << "\n";
    for (const auto &t : tets)
        f << "4 " << t.a << " " << t.b << " " << t.c << " " << t.d << "\n";
    f << "CELL_TYPES " << tets.size() << "\n";
    std::for_each(tets.begin(), tets.end(), [&f] (const tet &) { f << "10\n"; });
    f << "POINT_DATA " << pts.size() << "\n";
    f << "SCALARS u float\nLOOKUP_TABLE default\n";
    for (const auto &v: vals)
        f << v << "\n";
}

int edgeid(int u, int v) {
    return std::min(u, v) * 100 + std::max(u, v);
}

int main() {
    std::vector<point> pts;

    for (int z = 0; z < 2; z++)
        for (int y = 0; y < 2; y++)
            for (int x = 0; x < 2; x++)
                pts.emplace_back(point(x, y, z));

    for (int f = 0; f < 6; f++) {
        point sum(0, 0, 0);
        for (int v = 0; v < 4; v++)
            sum += pts[cubefaces[f][v]];
        pts.push_back(sum * 0.25);
    }

    pts.emplace_back(point(.5, .5, .5));

    std::vector<tet> tets;

    for (int f = 0; f < 6; f++)
        for (int j = 0; j < 4; j++)
            tets.emplace_back(tet{cubefaces[f][(j+1)&3], cubefaces[f][j], 8+f, 14});

    std::unordered_map<int, int> edgemap;
    std::vector<int> iedgemap;
    for (const auto &t : tets) {
        int p[4] = {t.a, t.b, t.c, t.d};
        int i[6] = {0, 0, 0, 1, 1, 2};
        int j[6] = {1, 2, 3, 2, 3, 3};
        for (int k = 0; k < 6; k++) {
            int id = edgemap.size();
            int eid = edgeid(p[i[k]], p[j[k]]);
            if (edgemap.find(eid) == edgemap.end()) {
                edgemap[eid] = id;
                iedgemap.push_back(eid);
            }
        }
    }

    std::vector<std::unordered_set<int>> graph(edgemap.size());

    for (const auto &t : tets) {
        int p[4] = {t.a, t.b, t.c, t.d};
        int i[6] = {0, 0, 0, 1, 1, 2};
        int j[6] = {1, 2, 3, 2, 3, 3};
        for (int k = 0; k < 6; k++) {
            int eid = edgemap[edgeid(p[i[k]], p[j[k]])];
            for (int ks = k + 1; ks < 6; ks++) {
                int eids = edgemap[edgeid(p[i[ks]], p[j[ks]])];
                graph[eid].insert(eids);
                graph[eids].insert(eid);
            }
        }
    }

    std::vector<double> vals(pts.size());

    double theta = .7;

    std::vector<int> marks(edgemap.size());

    for (int cases = 0; cases < (1 << 8); cases++) {
        int lo, hi;

        for (int k = 0; k < 8; k++)
            vals[k] = (cases & (1 << k)) ? 1 : 0;

        vals[14] = 0;
        for (int k = 0; k < 8; k++)
            vals[14] += vals[k];
        vals[14] *= 0.125;

        for (int f = 0; f < 6; f++) {
            vals[8 + f] = 0;
            for (int j = 0; j < 4; j++)
                vals[8 + f] += vals[cubefaces[f][j]];
            vals[8 + f] *= 0.25;
        }

        for (int e = 0; e < edgemap.size(); e++) {
            int eid = iedgemap[e];
            int u = eid % 100;
            int v = eid / 100;
            marks[e] = (vals[u] > theta) ^ (vals[v] > theta);
        }

        std::ofstream dot("cube." + std::to_string(cases) + ".dot", std::ios::binary);

        dot << "strict graph cube {\n";
        for (int u = 0; u < graph.size(); u++) {
            if (marks[u] == 0)
                continue;
            int eu = iedgemap[u];
            if ((eu / 100) < 8 && (eu % 100) < 8) {
                dot << "e" << eu / 100 << "_" << eu % 100;
                dot << " [style = bold];\n";
            }
            dot << "e" << eu / 100 << "_" << eu % 100 << " -- {";
            for (int v : graph[u]) {
                int ev = iedgemap[v];
                if (marks[v] == 0)
                    continue;
                dot << "e" << ev / 100 << "_" << ev % 100 << " ";
            }
            dot << "}\n";
        }
        dot << "}\n";

        save("cube." + std::to_string(cases) + ".vtk", pts, tets, vals);
    }

    return 0;
}
