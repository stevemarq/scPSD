import scPSD
import argparse
from scipy.io import mmread
import fast_matrix_market as fmm

import numpy as np
import os


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
    parser.add_argument("-s",
                        "--sparce",
                        type=bool,
                        default=False,
                        help="Whether dataset file used is sparse matrix or not")

    args = parser.parse_args()
    dataset_loc = args.dataset_loc
    clustering_mthd = args.clustering
    sparce = args.sparce

    if clustering_mthd not in ["TSNE", "PCA", "UMAP"]:
        raise ValueError("Invalid clustering algorithm. Options: SVM, KNN")

    if sparce:
        try:
            data = fmm.mmread(dataset_loc)
            data = data.toarray()
        except FileNotFoundError:
            print("Unable to locate" + dataset_loc)
        except ValueError:
            print("Not a Matrix Market file. Missing banner.")
    else:
        raise NotImplementedError

    algo = scPSD(data)
    transformed = algo.transform()



if __name__ == "__main__":
    main()



