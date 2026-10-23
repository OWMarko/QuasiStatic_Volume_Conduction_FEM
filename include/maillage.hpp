#pragma once

#include <Eigen/Dense>
#include <array>
#include <vector>

struct Maillage
{
    std::vector<Eigen::Vector3d> noeuds;
    std::vector<std::array<int, 4>> tetraedres;
};

Maillage creer_maillage(int n, double L);
