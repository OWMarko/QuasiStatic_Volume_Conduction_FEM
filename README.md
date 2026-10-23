# Finite Elements for 3D Volume Conduction

A small C++ solver for quasi-static electrical conduction in a cube, using P1 tetrahedral finite elements.

## Model

The total current and its conservation are given by

\[
J=J_p-\sigma\nabla\phi, \qquad \nabla\cdot J=0
\]

For an insulated domain, \(n\cdot J=0\). The corresponding weak formulation is

\[
\int_\Omega\sigma\nabla\phi\cdot\nabla v\,dx
=\int_\Omega J_p\cdot\nabla v\,dx
\]

On each tetrahedron, the basis functions are linear, and we assemble

\[
K^T_{ij}=|T|\sigma\nabla N_i\cdot\nabla N_j
\]

The Neumann matrix is singular: we set the potential at the first node to zero, then solve the reduced system using the conjugate gradient method.

## Project Structure

```text
include/maillage.hpp    cube mesh
include/solveur.hpp     FEM problem declarations
src/maillage.cpp        subdivision into tetrahedra
src/solveur.cpp         stiffness matrix, right-hand side, CG, fields and CSV export
src/main.cpp            parameters and program execution
scripts/plot.py         plot of a 2D slice of the potential
```

The main parameters are defined at the beginning of `src/main.cpp`.

## Build

Dependencies: C++17 compiler, CMake, Eigen3.  
For plotting: Python 3, NumPy and Matplotlib.

On Ubuntu/Debian:

```bash
sudo apt install g++ cmake libeigen3-dev python3-numpy python3-matplotlib
```

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

From **the project root directory**, run:

```bash
./build/conduction primaire 8
python3 scripts/plot.py
```

The program writes `data/potentiel.csv` and `data/champs.csv`, and the script generates `data/potentiel.png`.

On Windows, run from the project root directory:

```powershell
.\build\Release\conduction.exe primaire 8
python scripts/plot.py
```

Eigen3 must be installed and accessible to CMake (for example, using `-DEIGEN3_INCLUDE_DIR=...`).

## What You Can Modify

In `src/main.cpp`:

- `L`: cube side length, in metres;
- `sigma`: constant conductivity, in S/m;
- `intensite`: primary current amplitude, in A/m²;
- `n`: subdivisions per direction (also configurable as an argument, from 2 to 40);
- `mode`: `primaire` or `manufacture`.

Examples:

```bash
./build/conduction primaire 12
./build/conduction manufacture 4
./build/conduction manufacture 8
./build/conduction manufacture 12
```

The `primaire` mode uses a smooth, localized current that vanishes on the faces of the cube. The `manufacture` mode uses the known solution

\[
\phi(x,y,z)=10^{-3}\cos(\pi x/L)\cos(\pi y/L)\cos(\pi z/L)
\]

with \(-\nabla\cdot(\sigma\nabla\phi)=f\), \(f=3\sigma(\pi/L)^2\phi\) and \(\nabla\phi\cdot n=0\), to observe how the nodal RMS error decreases as the mesh is refined. The right-hand side is integrated using four quadrature points per tetrahedron; this quadrature is not exact for the trigonometric function.

In `champs.csv`, the `primaire` mode exports the total current \(J_p-\sigma\nabla\phi_h\). The `manufacture` mode only exports the conduction current \(-\sigma\nabla\phi_h\), associated with the volume source \(f\).

The exported electric fields are constant on each P1 tetrahedron, and the reconstructed current is not locally conservative in the strong sense.

## License

MIT, see `LICENSE`.