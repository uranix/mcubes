#pragma once

#include "point.h"

struct edge {
    const int i, j, k, dir;
    edge(int i, int j, int k, int dir)
        : i(i), j(j), k(k), dir(dir)
    { }
    static const edge fromid(int id) {
        switch (id) {
            case 0: 
        }
    }
    const edge shift(int i, int j, int k) const {
        return edge(i + this->i, j + this->j, k + this->k, dir);
    }
};

struct triangle {
    int p1, p2, p3;
};

struct triedge {
    edge p1, p2, p3;
};
