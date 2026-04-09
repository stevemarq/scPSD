import scPSD
import scPSD_preproc as preproc


parser = argparse.ArgumentParser(
                    prog='scPSD clustering',
                    description='What the program does')

parser.add_argument('filename',
                    help="scRNA-seq file must be in the ... format",
                    type=str)
parser.add_argument('preprocesss_mthd',
                    default="Raw",
                    # choices=['Raw', 'TMM', 'CPM', 'Scone', 'Linnorm', 'Scran', 'Seurat'],
                    type=str,
                    help="Type of scRNA-seq normalization to be done before clustering. DEFAULT: Raw")
parser.add_argument('clustering_mthd',
                    default='SVM',
                    type=str,
                    help="Type of clustering algorithm used. DEFAULT: SVM")

args = parser.parse_args()

filename = args.filename
norm_type = args.preprocess
clustering_mthd = args.clustering

if norm_type not in ['Raw', 'TMM', 'CPM', 'Scone', 'Linnorm', 'Scran', 'Seurat']:
    raise ValueError("Invalid normalization type. Options: Raw, TMM, CPM, Scone, and Linnorm")

if clustering_mthd not in ['SVM', 'KNN']:
    raise ValueError("Invalid clustering algorithm. Options: SVM, KNN")



