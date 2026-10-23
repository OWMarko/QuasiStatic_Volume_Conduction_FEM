# Éléments finis pour la conduction volumique en 3D

Petit solveur C++ de la conduction électrique quasi-statique dans un cube, avec éléments finis tétraédriques P1. L'objectif est de traduire directement les équations en un calcul numérique simple, puis de visualiser le potentiel.

## Modèle

Le courant total et sa conservation s'écrivent

\[
J=J_p-\sigma\nabla\phi, \qquad \nabla\cdot J=0.
\]

Pour un domaine isolé, \(n\cdot J=0\). La formulation faible correspondante est

\[
\int_\Omega\sigma\nabla\phi\cdot\nabla v\,dx
=\int_\Omega J_p\cdot\nabla v\,dx.
\]

Sur chaque tétraèdre, les fonctions de base sont linéaires et l'on assemble

\[
K^T_{ij}=|T|\sigma\nabla N_i\cdot\nabla N_j.
\]

La matrice de Neumann est singulière : on fixe le potentiel du premier nœud à zéro, puis on résout le système réduit par gradient conjugué (Eigen). Les calculs géométriques, l'assemblage et la source sont écrits dans le projet.

## Organisation

```text
include/maillage.hpp    maillage du cube
include/solveur.hpp    declarations du probleme FEM
src/maillage.cpp        subdivision en tetraedres
src/solveur.cpp        rigidite, second membre, CG, champs et export CSV
src/main.cpp           parametres et lancement du calcul
scripts/plot.py        figure d'une coupe 2D du potentiel
```

Les paramètres de base sont visibles au début de `src/main.cpp`.

## Compilation

Dépendances : compilateur C++17, CMake, Eigen3. Pour les figures : Python 3, NumPy et Matplotlib.

Sous Ubuntu/Debian :

```bash
sudo apt install g++ cmake libeigen3-dev python3-numpy python3-matplotlib
```

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Depuis **la racine du projet**, exécuter :

```bash
./build/conduction primaire 8
python3 scripts/plot.py
```

Le programme écrit `data/potentiel.csv`, `data/champs.csv` et le script produit `data/potentiel.png`.

Sous Windows, avec un générateur Visual Studio, lancer depuis la racine :

```powershell
.\build\Release\conduction.exe primaire 8
python scripts/plot.py
```

Eigen3 doit être installé et accessible à CMake (par exemple via `-DEIGEN3_INCLUDE_DIR=...`).

## Ce que l'on peut modifier

Dans `src/main.cpp` :

- `L` : côté du cube, en mètres ;
- `sigma` : conductivité constante, en S/m ;
- `intensite` : amplitude du courant primaire, en A/m² ;
- `n` : subdivisions par direction (également modifiable en argument, de 2 à 40) ;
- `mode` : `primaire` ou `manufacture`.

Exemples :

```bash
./build/conduction primaire 12
./build/conduction manufacture 4
./build/conduction manufacture 8
./build/conduction manufacture 12
```

Le mode `primaire` impose un courant localisé, lisse et nul sur les faces du cube. Le mode `manufacture` utilise la solution connue

\[
\phi(x,y,z)=10^{-3}\cos(\pi x/L)\cos(\pi y/L)\cos(\pi z/L),
\]

avec \(-\nabla\cdot(\sigma\nabla\phi)=f\), \(f=3\sigma(\pi/L)^2\phi\) et \(\nabla\phi\cdot n=0\), afin d'observer comment diminue **l'erreur RMS aux nœuds** lorsque le maillage est raffiné. Ce n'est pas une norme \(L^2\) intégrée sur le domaine. La quadrature du second membre comporte quatre points par tétraèdre ; elle n'est pas exacte pour la fonction trigonométrique.

Dans `champs.csv`, le mode `primaire` exporte le courant total \(J_p-\sigma\nabla\phi_h\). Le mode `manufacture` exporte uniquement le courant de conduction \(-\sigma\nabla\phi_h\), associé à la source volumique \(f\).

Lancer de nouvelles expériences écrase les deux CSV de `data/` : conserver ailleurs les fichiers intéressants avant une nouvelle exécution.

## Limites

Le domaine et les courants sont synthétiques : il ne s'agit pas d'un modèle anatomique, d'une simulation clinique ou d'une prédiction physiologique. Seule la conductivité scalaire constante est implémentée. Les champs électriques exportés sont constants sur chaque tétraèdre P1 et le courant reconstruit n'est pas localement conservatif au sens fort. Cette version est volontairement pédagogique ; elle ne reprend pas la validation, la robustesse numérique et les variantes de calcul du dépôt scientifique complet.

## Licence

MIT, voir `LICENSE`.
