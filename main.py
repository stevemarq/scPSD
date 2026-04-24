import scPSD as scpsd
import visualize

import pandas as pd
import argparse
from scipy.io import mmread
#pip install fast-matrix-market

import fast_matrix_market as fmm
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, HDBSCAN
from sklearn.metrics import silhouette_score, calinski_harabasz_score as VRC

import umap #pip install umap-learn
import random
import numpy as np

random.seed(42)

def evaluate_clusters(x, labels):
    sil_score = silhouette_score(x, labels, metric='euclidean')
    vrc_score = VRC(x, labels)
    return sil_score, vrc_score


def transform(x, cell_classes, clustering_mthd):
    psd = scpsd.scPSD()
    print("Starting scPSD transformation")
    transformed = psd.transform(x)
    # transformed = x
    # transformed = (x - np.min(x)) / (np.max(x) - np.min(x))
    n_clusters = np.unique(cell_classes).size

    print("Starting downstream analysis ")
    clusters = None
    if clustering_mthd == "TSNE":
        clusters = TSNE(n_components=2).fit_transform(transformed.T)
    elif clustering_mthd == "PCA":
        clusters = PCA(n_components=2).fit_transform(transformed.T)
    elif clustering_mthd == "UMAP":
        clusters = umap.UMAP(n_components=2,
                             n_neighbors=30,
                             min_dist=0.0).fit_transform(transformed.T)
        # centers = HDBSCAN(min_samples=n_clusters, min_cluster_size=500).fit(clusters)
    kmeans = KMeans(n_clusters=n_clusters).fit(clusters)
    labels = kmeans.labels_
    sil, vrc = evaluate_clusters(clusters, labels)

    results = {"Clusters": clusters, "Scores": [sil, vrc]}
    return results

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
                        help="Location of dataset labels for classification tasks")

    parser.add_argument("-s",
                        "--sparse",
                        type=bool,
                        default=True,
                        help="Whether dataset file used is sparse matrix or not")

    args = parser.parse_args()
    dataset_loc = args.dataset_loc
    classes = args.labels
    clustering_mthd = (args.clustering_mthd).upper()
    sparse = args.sparse

    if clustering_mthd not in ["TSNE", "PCA", "UMAP"]:
        raise ValueError("Invalid clustering algorithm. Options: SVM, KNN")


    # fix this ! looks ugly
    data = None
    if sparse:
        try:
            data = fmm.mmread(dataset_loc)
        except FileNotFoundError:
            print("Unable to locate" + dataset_loc)
        except ValueError:
            print("Not a Matrix Market file. Missing banner.")
    else:
        raise NotImplementedError

    data = data.toarray()
    cell_classes = pd.read_csv(classes)
    cell_classes = np.asarray(cell_classes).flatten()
    results = transform(data, cell_classes, clustering_mthd)

    # visualize.main(clusters, labels, cell_classes)
    visualize.main(results["Clusters"], cell_classes, clustering_mthd, results["Scores"])
    print("Analysis complete")


if __name__ == "__main__":
    main()



# python main.py -c "TSNE" -d "../scPSD/processed_data/prep_out.mtx" -l "../scPSD/processed_data/cells_out.csv" -s True