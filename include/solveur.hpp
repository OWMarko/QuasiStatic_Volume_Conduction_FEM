#pragma once

#include "maillage.hpp"

#include <Eigen/Dense>

struct Element
{
    double volume;
    std::array<Eigen::Vector3d, 4> gradient;
};

Element geometrie(const Maillage& maillage, const std::array<int, 4>& tetraedre);
Eigen::Vector3d courant_primaire(const Eigen::Vector3d& x, double L, double intensite);

struct Resultat
{
    Eigen::VectorXd potentiel;
    int iterations = 0;
    double residu = 0.0;
    double erreur_nodale = 0.0;
};

// Mode manufacture : -div(sigma grad(phi)) = f
// Mode primaire : div(Jp - sigma grad(phi)) = 0
Resultat resoudre(const Maillage& maillage, double L, double sigma,
                 double intensite, bool manufacture);
