#include "output.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

void exporter_resultats(const Maillage& maillage, const Resultat& resultat,
                        double L, double sigma, double intensite, bool manufacture)
{
    std::filesystem::create_directories("data");

    // Potentiel aux noeuds.
    std::ofstream potentiel("data/potentiel.csv");
    if(!potentiel)
    {
        throw std::runtime_error("Impossible de creer data/potentiel.csv");
    }
    potentiel << std::setprecision(16);
    potentiel << "x_m,y_m,z_m,phi_V\n";

    for(std::size_t i = 0; i < maillage.noeuds.size(); ++i)
    {
        const auto& x = maillage.noeuds[i];
        potentiel << x.x() << ',' << x.y() << ',' << x.z() << ','
                  << resultat.potentiel[static_cast<Eigen::Index>(i)] << '\n';
    }
    potentiel.close();
    if(!potentiel)
    {
        throw std::runtime_error("Echec de l'ecriture de data/potentiel.csv");
    }

    // Champs au barycentre de chaque tetraedre : E = -grad(phi), J = Jp + sigma E.
    std::ofstream champs("data/champs.csv");
    if(!champs)
    {
        throw std::runtime_error("Impossible de creer data/champs.csv");
    }
    champs << std::setprecision(16);
    champs << "x_m,y_m,z_m,E_x_V_m,E_y_V_m,E_z_V_m,J_x_A_m2,J_y_A_m2,J_z_A_m2\n";

    for(const auto& tetraedre : maillage.tetraedres)
    {
        const Element element = geometrie(maillage, tetraedre);
        Eigen::Vector3d centre = Eigen::Vector3d::Zero();
        Eigen::Vector3d gradient_potentiel = Eigen::Vector3d::Zero();
        for(int i = 0; i < 4; ++i)
        {
            centre += maillage.noeuds[tetraedre[i]] / 4.0;
            gradient_potentiel += resultat.potentiel[tetraedre[i]] * element.gradient[i];
        }

        const Eigen::Vector3d champ_electrique = -gradient_potentiel;
        Eigen::Vector3d courant = sigma * champ_electrique;
        if(!manufacture)
        {
            courant += courant_primaire(centre, L, intensite);
        }

        champs << centre.x() << ',' << centre.y() << ',' << centre.z();
        for(int i = 0; i < 3; ++i)
        {
            champs << ',' << champ_electrique[i];
        }
        for(int i = 0; i < 3; ++i)
        {
            champs << ',' << courant[i];
        }
        champs << '\n';
    }
    champs.close();
    if(!champs)
    {
        throw std::runtime_error("Echec de l'ecriture de data/champs.csv");
    }
}
