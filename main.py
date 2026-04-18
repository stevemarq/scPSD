import scPSD
import scPSD_preproc as preproc
import argparse

def main():
    parser = argparse.ArgumentParser(
                        prog='scPSD clustering',
                        description='What the program does')

    parser.add_argument('filename',
                        help="scRNA-seq file must be in the ... format",
                        type=str)
    parser.add_argument('clustering_mthd',
                        default='SVM',
                        type=str,
                        help="Type of clustering algorithm used. DEFAULT: SVM")

    args = parser.parse_args()

    filename = args.filename
    norm_type = args.preprocess
    clustering_mthd = args.clustering

    if clustering_mthd not in ['SVM', 'KNN', "TSNE", "PCA", "UMAP"]:
        raise ValueError("Invalid clustering algorithm. Options: SVM, KNN")


if __name__ == "__main__":
    main()



