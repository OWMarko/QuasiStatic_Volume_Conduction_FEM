from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    fichier = Path('data/potentiel.csv')
    valeurs = np.genfromtxt(fichier, delimiter=',', names=True)
    x, y, z = valeurs['x_m'], valeurs['y_m'], valeurs['z_m']
    phi = valeurs['phi_V']

    niveau = np.unique(z)[len(np.unique(z)) // 2]
    coupe = np.isclose(z, niveau)
    xs, ys = np.unique(x[coupe]), np.unique(y[coupe])
    image = np.full((len(ys), len(xs)), np.nan)
    for xi, yi, pi in zip(x[coupe], y[coupe], phi[coupe]):
        i = np.argmin(np.abs(xs - xi))
        j = np.argmin(np.abs(ys - yi))
        image[j, i] = pi * 1e3

    plt.figure(figsize=(6, 5))
    plt.imshow(image, origin='lower', extent=[xs[0], xs[-1], ys[0], ys[-1]],
               aspect='equal', cmap='coolwarm')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title(f'Potentiel a z = {niveau:.4g} m')
    plt.colorbar(label='Potentiel (mV)')
    plt.tight_layout()
    sortie = fichier.with_suffix('.png')
    plt.savefig(sortie, dpi=180)
    print(f'Figure : {sortie}')


if __name__ == '__main__':
    main()
