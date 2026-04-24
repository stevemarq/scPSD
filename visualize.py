import matplotlib.pyplot as plt
import numpy as np

# assumes scores are in [sil, vrc] form
def main(points, cell_class, clustering_mthd, scores = None):
    fig, ax = plt.subplots()
    x = points[:, 0]
    y = points[:, 1]

    colors = plt.cm.rainbow(np.linspace(0, 1, len(np.unique(cell_class))))
    for color, label in zip(colors, np.unique(cell_class)):
        mask = cell_class == label
        ax.scatter(x[mask], y[mask], color = color, label=label)


    ax.legend(title="Cell Types", loc="best", bbox_to_anchor=(1, 1), borderaxespad=0.0)
    ax.set_xlabel(f'{clustering_mthd} 1')
    ax.set_ylabel(f'{clustering_mthd} 2')
    if scores is not None:
        ax.set_title(f'SS: {scores[0]:.3f}, VRC: {scores[1]:.3f}')
    plt.show()

