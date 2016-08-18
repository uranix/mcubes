#pragma once

#include <cmath>
#include <cstddef>
#include <cassert>
#include <functional>

struct dim3 {
    std::ptrdiff_t i, j, k;

    struct hash {
        size_t operator()(const dim3 &o) const {
            size_t i = o.i;
            size_t j = o.j;
            size_t k = o.k;
            std::hash<size_t> hasher;
            return hasher(hasher(hasher(i) ^ j) ^ k);
        }
    };
    dim3(const std::ptrdiff_t i, const std::ptrdiff_t j, const std::ptrdiff_t k) : i(i), j(j), k(k) { }
    dim3() : dim3(0, 0, 0) { }
    dim3(const dim3 &o) : dim3(o.i, o.j, o.k) { }
    bool operator==(const dim3 &o) const {
        return (i == o.i) && (j == o.j) && (k == o.k);
    }
    dim3 &operator=(const dim3 &o) {
        i = o.i;
        j = o.j;
        k = o.k;
        return *this;
    }
    std::ptrdiff_t linear_index(const dim3 &n) const {
        assert(i < n.i && j < n.j && k < n.k && i >= 0 && j >= 0 && k >= 0);
        return i + n.i * (j + n.j * k);
    }
    const dim3 operator+(const dim3 &o) const {
        return dim3(i + o.i, j + o.j, k + o.k);
    }
    const dim3 operator-(const dim3 &o) const {
        return dim3(i - o.i, j - o.j, k - o.k);
    }
    dim3 &operator+=(const dim3 &o) {
        *this = *this + o;
        return *this;
    }
    dim3 &operator-=(const dim3 &o) {
        *this = *this - o;
        return *this;
    }
    static const dim3 X() {
        return dim3(1, 0, 0);
    }
    static const dim3 Y() {
        return dim3(0, 1, 0);
    }
    static const dim3 Z() {
        return dim3(0, 0, 1);
    }
    static const dim3 cube_vertex(int v) {
        assert(v < 8 && v >= 0);
        return dim3(v & 1, (v & 2) >> 1, (v & 4) >> 2);
    }
    int vertex_id() const {
        assert(i < 2 && j < 2 && k < 2 && i >= 0 && j >= 0 && k >= 0);
        return i + 2 * j + 4 * k;
    }
};

struct point {
    double x, y, z;
    point() { }
    point(const double x, const double y, const double z) : x(x), y(y), z(z) { }
    point(const point &p) : point(p.x, p.y, p.z) { }
    explicit point(const dim3 &p, const point &ll = point(0, 0, 0), const point &h = point(1, 1, 1))
        : point(ll.x + (p.i + 0.5) * h.x, ll.y + (p.j + 0.5) * h.y, ll.z + (p.k + 0.5) * h.z)
    { }
    double dot(const point &p) const {
        return x * p.x + y * p.y + z * p.z;
    }
    const point cross(const point &p) const {
        return point(
                y * p.z - p.y * z,
                z * p.x - p.z * x,
                x * p.y - p.x * y
            );
    }
    const point operator-(const point &p) const {
        return point(
                x - p.x,
                y - p.y,
                z - p.z
            );
    }
    const point operator+(const point &p) const {
        return point(
                x + p.x,
                y + p.y,
                z + p.z
            );
    }
    point &operator-=(const point &p) {
        x -= p.x;
        y -= p.y;
        z -= p.z;
        return *this;
    }
    point &operator+=(const point &p) {
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }
    point &operator*=(const double m) {
        x *= m;
        y *= m;
        z *= m;
        return *this;
    }
    const double norm() const {
        return std::sqrt(dot(*this));
    }
    const point operator*(const double m) const {
        return point(m * x, m * y, m * z);
    }
};
