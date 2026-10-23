#pragma once

#include "solveur.hpp"

void exporter_resultats(const Maillage& maillage, const Resultat& resultat,
                       double L, double sigma, double intensite,
                       bool manufacture);