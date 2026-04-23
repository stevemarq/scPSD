import scPSD
import visualize


import pandas as pd

import argparse
from scipy.io import mmread
#pip install fast-matrix-market

import fast_matrix_market as fmm
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, HDBSCAN
import umap #pip install umap-learn


import random
import numpy as np
import os

random.seed(42)
def main():
    parser = argparse.ArgumentParser(
                        prog='scPSD clustering',
                        description='Performing power spectral density analysis on scRNA data')

    parser.add_argument("-c",
                        '--clustering_mthd',
                        default='SVM',
                        type=str,
                        help="Type of clustering algorithm used. DEFAULT: PCA")

    parser.add_argument("-d",
                        "--dataset_loc",
                        type=str,
                        help="Location of the datasets")

    parser.add_argument("-l",
                        "--labels",
                        type=str,
                        help"Location of dataset labels for classification tasks")

    parser.add_argument("-s",
                        "--sparse",
                        type=bool,
                        default=True,
                        help="Whether dataset file used is sparse matrix or not")

    args = parser.parse_args()
    dataset_loc = args.dataset_loc
    classes = args.labels
    clustering_mthd = args.clustering
    sparse = args.sparse

    if clustering_mthd not in ["TSNE", "PCA", "UMAP"]:
        raise ValueError("Invalid clustering algorithm. Options: SVM, KNN")

    if sparse:
        try:
            data = fmm.mmread(dataset_loc)
            data = data.toarray()
        except FileNotFoundError:
            print("Unable to locate" + dataset_loc)
        except ValueError:
            print("Not a Matrix Market file. Missing banner.")
    else:
        raise NotImplementedError

    psd = scPSD()
    transformed = psd.transform(data)
    cell_classes = pd.read_csv(classes)
    cell_classes = np.asarray(cell_classes)
    n_clusters = np.unique(cell_classes).size


    if clustering_mthd == "TSNE":
        clusters = TSNE(n_components=2).fit_transform(transformed)
        centers = KMeans(n_clusters=n_clusters).fit(clusters)
    elif clustering_mthd == "PCA":
        clusters = PCA(n_components=2).fit_transform(transformed)
        centers = KMeans(n_clusters=n_clusters).fit(clusters)
    elif clustering_mthd == "UMAP":
        clusters = umap.UMAP(n_components=2,
                             n_neighbors=30,
                             min_dist = 0.0).fit_transform(transformed)
        centers = HDBSCAN(min_samples=n_clusters,
                         min_cluster_size=500).fit(clusters)
    labels = centers.labels_

    visualize.main(clusters, labels, cell_classes)


if __name__ == "__main__":
    main()



