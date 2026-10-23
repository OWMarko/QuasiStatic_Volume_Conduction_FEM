#include "solveur.hpp"

#include <Eigen/Sparse>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

    constexpr double pi = 3.14159265358979323846;

    Element geometrie(const Maillage& maillage, const std::array<int, 4>& tetraedre)
    {
        Eigen::Matrix3d B;
        for(int j = 0; j < 3; ++j)
        {
            B.col(j) = maillage.noeuds[tetraedre[j + 1]] - maillage.noeuds[tetraedre[0]];
        }

        Element element;
        element.volume = B.determinant() / 6.0;
        if(!std::isfinite(element.volume) || element.volume <= 0.0)
        {
            throw std::runtime_error("Tetraedre de volume nul, negatif ou non fini");
        }

        // grad(Nj) = B^(-T) ej pour j = 1, 2, 3 la somme des gradients est nulle
        const auto decomposition = B.transpose().partialPivLu();

        for(int j = 1; j < 4; ++j)
        {
            Eigen::Vector3d direction = Eigen::Vector3d::Zero();
            direction[j - 1] = 1.0;
            element.gradient[j] = decomposition.solve(direction);
        }

        element.gradient[0] = -element.gradient[1] - element.gradient[2] - element.gradient[3];

        return element;
    }

    static double potentiel_exact(const Eigen::Vector3d& x, double L)
    {
        const double q = pi / L;
        return 1e-3 * std::cos(q*x.x()) * std::cos(q*x.y()) * std::cos(q*x.z());
    }

    Eigen::Vector3d courant_primaire(const Eigen::Vector3d& x, double L, double intensite)
    {
        const Eigen::Vector3d position = x/L;
        const Eigen::Vector3d centre(0.42, 0.52, 0.48);
        double enveloppe = 1.0;

        for(int i = 0; i < 3; ++i)
        {
            const double a = 4.0 * position[i] * (1.0-position[i]);
            enveloppe *= a*a;
        }

        const double largeur = 0.18;
        const double gaussienne = std::exp(-(position - centre).squaredNorm() / (2.0*largeur*largeur));
        return {intensite*enveloppe*gaussienne, 0.0, 0.0};
    }

Resultat resoudre(const Maillage& maillage, double L, double sigma,
                 double intensite, bool manufacture)
{
    if(!std::isfinite(sigma) || sigma <= 0.0)
    {
        throw std::invalid_argument("La conductivite doit etre positive et finie");
    }

    const int N = static_cast<int>(maillage.noeuds.size());
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
    std::vector<Eigen::Triplet<double>> termes;
    termes.reserve(maillage.tetraedres.size() * 16);

    // Quadrature a quatre points barycentriques
    const double a = 0.5854101966249685;
    const double c = 0.1381966011250105;
    const std::array<std::array<double, 4>, 4> points =
    {{
        {a, c, c, c},
        {c, a, c, c},
        {c, c, a, c},
        {c, c, c, a}
    }};

    for(const auto& tetraedre : maillage.tetraedres)
    {
        const Element element = geometrie(maillage, tetraedre);

        // K_ij^T = |T| sigma grad(Ni).grad(Nj)
        // La jauge phi(0) = 0 permet d'omettre la ligne et la colonne du noeud 0
        for(int i = 0; i < 4; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                if(tetraedre[i] > 0 && tetraedre[j] > 0)
                {
                    const double rigidite = element.volume * sigma * element.gradient[i].dot(element.gradient[j]);
                    termes.emplace_back(tetraedre[i] - 1, tetraedre[j] - 1, rigidite);
                }
            }
        }

        // Second membre : integrale de f Ni (manufacture) ou de Jp.grad(Ni)
        for(const auto& lambda : points)
        {
            Eigen::Vector3d x = Eigen::Vector3d::Zero();
            for(int j = 0; j < 4; ++j)
            {
                x += lambda[j] * maillage.noeuds[tetraedre[j]];
            }

            if(manufacture)
            {
                const double q = pi / L;
                const double f = 3.0 * sigma * q * q * potentiel_exact(x, L);
                for(int i = 0; i < 4; ++i)
                {
                    b[tetraedre[i]] += element.volume * 0.25 * f * lambda[i];
                }
            }
            else
            {
                const Eigen::Vector3d Jp = courant_primaire(x, L, intensite);
                for(int i = 0; i < 4; ++i)
                {
                    b[tetraedre[i]] += element.volume * 0.25 * Jp.dot(element.gradient[i]);
                }
            }
        }
    }

    // Aucun terme de bord.
    // La somme du second membre complet doit etre nulle
    if(std::abs(b.sum()) > 1e-9 * b.lpNorm<1>())
    {
        throw std::runtime_error("Second membre de Neumann incompatible");
    }

    Eigen::SparseMatrix<double> A(N - 1, N - 1);
    A.setFromTriplets(termes.begin(), termes.end());

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
    cg.setTolerance(1e-10);
    cg.setMaxIterations(3000);
    cg.compute(A);
    const Eigen::VectorXd potentiel_libre = cg.solve(b.tail(N - 1));
    
    if(cg.info() != Eigen::Success || !potentiel_libre.allFinite())
    {
        throw std::runtime_error("Le gradient conjugue n'a pas converge vers une solution finie");
    }

    Resultat resultat;
    resultat.potentiel = Eigen::VectorXd::Zero(N);
    resultat.potentiel.tail(N - 1) = potentiel_libre;
    resultat.iterations = cg.iterations();
    const double norme_b = b.tail(N - 1).norm();
    resultat.residu = (A * potentiel_libre - b.tail(N - 1)).norm() / (norme_b > 0.0 ? norme_b : 1.0);

    if(manufacture)
    {
        // Erreur RMS aux noeuds
        Eigen::VectorXd ecart(N);
        for(int i = 0; i < N; ++i)
        {
            ecart[i] = resultat.potentiel[i] - potentiel_exact(maillage.noeuds[i], L);
        }

        ecart.array() -= ecart.mean();
        resultat.erreur_nodale = ecart.norm() / std::sqrt(static_cast<double>(N));
    }

    return resultat;
}

