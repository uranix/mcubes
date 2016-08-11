#pragma once

#include "geom.h"

#include <string>
#include <vector>

void save(const std::string &filename, const std::vector<point> &pts, const std::vector<triangle> &tri);
void save(const std::string &filename, const std::vector<point> &pts, const std::vector<quad> &quads);
