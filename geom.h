#pragma once

#include "point.h"

struct triangle {
    int p1, p2, p3;
};

struct quad {
    int p1, p2, p3, p4;
};

struct mesh_vertex {
    dim3 cube; // cube coordinates
    int patchid; // 0-3
};

struct edge {
    dim3 beg;
    dim3 end;

    edge(int num) {
        int idir = num / 4;
        int a = num & 1;
        int b = (num & 2) >> 1;
        if (idir == 0) {
            beg.i = 0;
            end.i = 1;
            beg.j = end.j = a;
            beg.k = end.k = b;
        } else if (idir == 1) {
            beg.i = end.i = a;
            beg.j = 0;
            end.j = 1;
            beg.k = end.k = b;
        } else {
            beg.i = end.i = a;
            beg.j = end.j = b;
            beg.k = 0;
            end.k = 1;
        }
    }
};
