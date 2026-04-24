import matplotlib.pyplot as plt
import numpy as np

def main(points, cell_class):
    fig, ax = plt.subplots()
    x = points[:, 0]
    y = points[:, 1]
    for label in np.unique(cell_class):
        mask = np.array(cell_class) == label
        ax.scatter(x[mask], y[mask], label=label)

    ax.legend(title="Cell Types")
    plt.show()

