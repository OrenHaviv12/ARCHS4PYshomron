import archs4py as a4
import pandas as pd
import numpy as np  # Import numpy for mathematical operations
import argparse


def load_scrnaseq_data(file_path, num_samples=300):
    """
    Load random samples for scRNA-seq data from the dataset.
    """
    print(f"Loading {num_samples} random samples from: {file_path}")
    try:
        random_data = a4.data.rand(file_path, num_samples, remove_sc=False)
        print("Data successfully loaded.")
        return pd.DataFrame(random_data)
    except Exception as e:
        print(f"Error loading data: {e}")
        return None


def load_gene_names(file_path):
    """
    Load gene names from a file.
    """
    print(f"Loading gene names from: {file_path}")
    try:
        with open(file_path, 'r') as file:
            genes = [line.strip() for line in file if line.strip()]
        print(f"Loaded {len(genes)} genes.")
        return genes
    except Exception as e:
        print(f"Error loading gene names: {e}")
        return None


def normalize_data(df, method="log"):
    """
    Normalize scRNA-seq data.
    """
    print(f"Normalizing data using {method} normalization...")
    if method == "log":
        return df.apply(lambda x: np.log1p(x))  # Use numpy's log1p directly
    elif method == "z-score":
        return df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    else:
        print("Unsupported normalization method. Skipping normalization.")
        return df


def filter_low_expression(df, threshold=1.0):
    """
    Filter out genes with low expression.
    """
    print(f"Filtering genes with mean expression below {threshold}...")
    return df[df.mean(axis=1) > threshold]


def deduplicate_genes(df):
    """
    Deduplicate genes by summing their expression values.
    """
    print("Deduplicating genes...")
    return df.groupby(df.index).sum()


def filter_genes_for_csv(df, genes_of_interest):
    """
    Filter scRNA-seq data for specific genes to save as a CSV file.
    """
    print(f"Filtering data for genes: {genes_of_interest}...")
    found_genes = df.index.intersection(genes_of_interest)
    if not found_genes.empty:
        print(f"Data for genes found: {list(found_genes)}")
        return df.loc[found_genes]
    else:
        print("None of the specified genes were found in the dataset.")
        return None


def main():
    parser = argparse.ArgumentParser(description="Process scRNA-seq data and save filtered genes to CSV.")
    parser.add_argument(
        "--data_file_path",
        type=str,
        default=r"C:\Users\orenh\OneDrive\מסמכים\ביולוגיה חישובית על גווניה\humanDNA.h5",
        help="Path to the scRNA-seq data file",
    )
    parser.add_argument(
        "--genes_file_path",
        type=str,
        default=r"C:\Users\orenh\OneDrive\מסמכים\ARCHS4PY\GENESc.text",
        help="Path to the file containing gene names",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        default="filtered_expression_matrix.csv",
        help="Output CSV file path",
    )
    parser.add_argument(
        "--num_samples",
        type=int,
        default=300,
        help="Number of random samples to load from the dataset",
    )
    parser.add_argument(
        "--normalize_method",
        type=str,
        default="log",
        choices=["log", "z-score"],
        help="Normalization method to apply",
    )
    parser.add_argument(
        "--expression_threshold",
        type=float,
        default=1.0,
        help="Threshold for low-expression filtering",
    )
    parser.add_argument(
        "--skip_deduplication",
            action="store_true",
        help="Skip gene deduplication step",
    )

    args = parser.parse_args()

    # Step 1: Load the Data
    data = load_scrnaseq_data(args.data_file_path, num_samples=args.num_samples)
    if data is None:
        return

    # Step 2: Normalize Data
    data = normalize_data(data, method=args.normalize_method)

    # Step 3: Filter Low Expression
    data = filter_low_expression(data, threshold=args.expression_threshold)

    # Step 4: Deduplicate Genes
    if not args.skip_deduplication:
        data = deduplicate_genes(data)

    # Step 5: Load Genes of Interest
    genes_of_interest = load_gene_names(args.genes_file_path)
    if genes_of_interest is None:
        return

    # Step 6: Filter Genes and Save to CSV
    gene_data = filter_genes_for_csv(data, genes_of_interest)
    if gene_data is None:
        return

    print(f"Saving filtered data to {args.output_csv}...")
    gene_data.to_csv(args.output_csv)
    print("Data saved successfully.")


if __name__ == "__main__":
    main()