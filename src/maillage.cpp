#include "maillage.hpp"

#include <cmath>
#include <stdexcept>

Maillage creer_maillage(int n, double L)
{
    if(n < 2 || n > 40 || !std::isfinite(L) || L <= 0.0)
    {
        throw std::invalid_argument("Maillage : choisir 2 <= n <= 40 et L > 0 fini");
    }

    Maillage maillage;
    const double h = L / n;
    const int N = n + 1;

    // x, puis y, puis z.
    for(int k = 0; k <= n; ++k)
    {
        for(int j = 0; j <= n; ++j)
        {
            for(int i = 0; i <= n; ++i)
            {
                maillage.noeuds.emplace_back(i * h, j * h, k * h);
            }
        }
    }

    for(int k = 0; k < n; ++k)
    {
        for(int j = 0; j < n; ++j)
        {
            for(int i = 0; i < n; ++i)
            {
                // decalages suivant x, y, z
                const int n000 = i + N * (j + N * k);
                const int n100 = n000 + 1;
                const int n010 = n000 + N;
                const int n110 = n010 + 1;
                const int n001 = n000 + N * N;
                const int n101 = n001 + 1;
                const int n011 = n001 + N;
                const int n111 = n011 + 1;

                // Six tetraedres orientes positivement autour de n000 n111
                maillage.tetraedres.push_back({n000, n100, n110, n111});
                maillage.tetraedres.push_back({n000, n100, n111, n101});
                maillage.tetraedres.push_back({n000, n010, n111, n110});
                maillage.tetraedres.push_back({n000, n010, n011, n111});
                maillage.tetraedres.push_back({n000, n001, n101, n111});
                maillage.tetraedres.push_back({n000, n001, n111, n011});
            }
        }
    }

    return maillage;
}
