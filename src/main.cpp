#include "maillage.hpp"
#include "solveur.hpp"
#include "output.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
    try
    {
        const double L = 0.01;          // longueur du cube (m)
        const double sigma = 0.3;       // conductivite (S/m)
        const double intensite = 0.1;   // amplitude du courant primaire (A/m2)
        int n = 8;                      // subdivisions dans chaque direction
        std::string mode = "primaire";  // "primaire" ou "manufacture"

        if(argc >= 2){mode = argv[1];}
        if(argc >= 3){n = std::stoi(argv[2]);}
        if((mode != "primaire" && mode != "manufacture") || argc > 3)
        {
            std::cerr << "Usage : ./conduction [primaire|manufacture] [n]\n";
            return 1;
        }

        const bool manufacture = (mode == "manufacture");

        const Maillage maillage = creer_maillage(n, L);
        const Resultat resultat = resoudre(maillage, L, sigma, intensite, manufacture);
        exporter_resultats(maillage, resultat, L, sigma, intensite, manufacture);

        std::cout << "Maillage : " << maillage.noeuds.size() << " noeuds, "
                  << maillage.tetraedres.size() << " tetraedres\n";
        std::cout << "CG : " << resultat.iterations << " iterations\n";
        std::cout << "Residu relatif : " << std::scientific << resultat.residu << '\n';
        if(manufacture){std::cout << "Erreur RMS nodale : " << resultat.erreur_nodale << " V\n";}
        std::cout << "Resultats : data/potentiel.csv, data/champs.csv\n";
    }

    catch(const std::exception& erreur)
    {
        std::cerr << "Erreur : " << erreur.what() << '\n';
        return 1;
    }
    
    return 0;
}
