import matplotlib.pyplot as plt
import numpy as np

# assumes scores are in [sil, vrc] form
def main(points, cell_class, clustering_mthd, scores = None):
    fig, ax = plt.subplots()
    x = points[:, 0]
    y = points[:, 1]
    for label in np.unique(cell_class):
        mask = np.array(cell_class) == label
        ax.scatter(x[mask], y[mask], label=label)

    ax.legend(title="Cell Types")
    ax.set_xlabel(f'{clustering_mthd} 1')
    ax.set_ylabel(f'{clustering_mthd} 2')
    if scores is not None:
        ax.set_title(f'SS: {scores[0]}, VRC: {scores[1]}')
    plt.show()

