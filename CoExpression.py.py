import os
import archs4py as a4
import pandas as pd
import numpy as np
import argparse


def check_file_existence(file_path):
    """
    Check if a file exists and is accessible.
    """
    if not os.path.isfile(file_path):
        print(f"Error: File {file_path} does not exist.")
        return False
    return True


def explore_hdf5_structure(file_path):
    """
    Explore and list the structure of the ARCHS4 HDF5 file using h5py.
    """
    print(f"Exploring the HDF5 file structure for: {file_path}")
    try:
        if not check_file_existence(file_path):
            raise FileNotFoundError(f"HDF5 file {file_path} not found.")

        # Use h5py to explore the structure
        import h5py
        with h5py.File(file_path, 'r') as f:
            print("HDF5 File Structure:")
            structure = {}
            for key in f.keys():
                print(f"{key}: {list(f[key])}")
                structure[key] = list(f[key])
            return structure
    except Exception as e:
        print(f"Error exploring HDF5 structure: {e}")
        return None


def load_scrnaseq_data(file_path, num_samples=300):
    """
    Load random samples for scRNA-seq data from the dataset.
    """
    print(f"Loading {num_samples} random samples from: {file_path}")
    try:
        if not check_file_existence(file_path):
            raise FileNotFoundError(f"Data file {file_path} not found.")
        random_data = a4.data.rand(file_path, num_samples, remove_sc=False)
        if random_data is None:
            raise ValueError("Loaded data is None.")
        print("Data successfully loaded.")
        return pd.DataFrame(random_data)
    except Exception as e:
        print(f"Error loading data: {e}")
        return None


def normalize_data(data, method="log"):
    """
    Normalize scRNA-seq data.
    """
    print(f"Normalizing data using {method} normalization...")
    try:
        if data is None:
            raise ValueError("Input data for normalization is None.")
        if method == "log":
            normalized_data = data.apply(lambda x: np.log1p(x), axis=0)
        elif method == "z-score":
            normalized_data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
        else:
            raise ValueError(f"Unsupported normalization method: {method}")
        print("Normalization successful.")
        return normalized_data
    except Exception as e:
        print(f"Error during normalization: {e}")
        return None


def filter_low_expression(data, threshold=1.0):
    """
    Filter genes with mean expression below the threshold.
    """
    print(f"Filtering genes with mean expression below {threshold}...")
    try:
        initial_genes = data.shape[0]
        filtered_data = data[data.mean(axis=1) > threshold]
        remaining_genes = filtered_data.shape[0]
        print(f"Filtered genes: {initial_genes} -> {remaining_genes}")
        return filtered_data
    except Exception as e:
        print(f"Error during filtering: {e}")
        return data


def aggregate_genes(data):
    """
    Aggregate duplicate genes by summing their expression values.
    """
    print("Aggregating duplicate genes...")
    try:
        initial_rows = data.shape[0]
        aggregated_data = data.groupby(data.index).sum()
        final_rows = aggregated_data.shape[0]
        print(f"Aggregation complete: {initial_rows} -> {final_rows}")
        return aggregated_data
    except Exception as e:
        print(f"Error during aggregation: {e}")
        return data


def extract_metadata(file_path, term=None, output_file="metadata.xlsx"):
    """
    Extract metadata and save it to an Excel file. If no term is provided, extract all metadata.
    """
    print(f"Extracting metadata for term: {term if term else 'all terms'}")
    try:
        import h5py
        with h5py.File(file_path, 'r') as f:
            metadata = {}
            if term:
                if term in f["meta"]:
                    metadata = f["meta"][term][:]
                else:
                    raise ValueError(f"Metadata term '{term}' not found.")
            else:
                for key in f["meta"].keys():
                    metadata[key] = f["meta"][key][:]

            metadata_df = pd.DataFrame(metadata)
            metadata_df.to_excel(output_file, index=False)
            print(f"Metadata saved to {output_file}")
            return metadata_df
    except Exception as e:
        print(f"Error extracting metadata: {e}")
        return None


def evaluate_coexpression(data, primary_gene, target_genes, method="pearson"):
    """
    Evaluate co-expression between a primary gene and target genes.

    Parameters:
        data (pd.DataFrame): The gene expression dataset.
        primary_gene (str): The name of the primary gene.
        target_genes (list of str): A list of target gene names.
        method (str): Correlation method, either "pearson" or "spearman".

    Returns:
        pd.DataFrame: Co-expression values for the target genes with the primary gene.
    """
    print(f"Evaluating co-expression using {method} correlation...")
    try:
        if data is None or data.empty:
            raise ValueError("Expression data is empty or None.")

        # Convert all gene names to lowercase for case-insensitive matching
        data.index = data.index.str.lower()
        primary_gene = primary_gene.lower()
        target_genes = [gene.lower() for gene in target_genes]

        if primary_gene not in data.index:
            raise ValueError(f"Primary gene '{primary_gene}' not found in the dataset.")

        missing_genes = [gene for gene in target_genes if gene not in data.index]
        if missing_genes:
            print(f"Warning: The following target genes were not found: {missing_genes}")
            target_genes = [gene for gene in target_genes if gene in data.index]

        if not target_genes:
            raise ValueError("No valid target genes found in the dataset.")

        # Retrieve expression values
        primary_expression = data.loc[primary_gene]
        target_expression = data.loc[target_genes]

        # Calculate correlation
        if method == "pearson":
            coexpression = target_expression.apply(lambda x: np.corrcoef(primary_expression, x)[0, 1], axis=1)
        elif method == "spearman":
            from scipy.stats import spearmanr
            coexpression = target_expression.apply(lambda x: spearmanr(primary_expression, x)[0], axis=1)
        else:
            raise ValueError(f"Unsupported correlation method: {method}")

        # Create results DataFrame
        results = pd.DataFrame({
            "Target Gene": target_genes,
            "Coexpression": coexpression
        })
        print("Co-expression evaluation complete.")
        return results
    except Exception as e:
        print(f"Error evaluating co-expression: {e}")
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
        "--num_samples", type=int, default=300, help="Number of random samples to load from the dataset"
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
        "--metadata_term",
        type=str,
        help="Metadata term to extract and save",
    )

    args = parser.parse_args()

    # Validate file paths
    if not check_file_existence(args.data_file_path):
        print("Data file path is invalid. Exiting.")
        return
    if not check_file_existence(args.genes_file_path):
        print("Genes file path is invalid. Exiting.")
        return

    # Step 1: Explore HDF5 Structure
    explore_hdf5_structure(args.data_file_path)

    # Step 2: Load the Data
    data = load_scrnaseq_data(args.data_file_path, num_samples=args.num_samples)
    if data is None:
        return

    # Step 3: Normalize Data
    data = normalize_data(data, method=args.normalize_method)
    if data is None:
        return

    # Step 4: Filter Low Expression
    data = filter_low_expression(data, threshold=args.expression_threshold)

    # Step 5: Aggregate Genes
    data = aggregate_genes(data)

    # Step 6: Extract Metadata (if provided)
    if args.metadata_term:
        extract_metadata(args.data_file_path, args.metadata_term)

    # Step 7: Co-expression Evaluation
    primary_gene = input("Enter the primary gene name: ").strip()
    target_genes = input("Enter target gene names (comma-separated): ").split(",")
    target_genes = [gene.strip() for gene in target_genes]

    coexpression_results = evaluate_coexpression(data, primary_gene, target_genes)
    if coexpression_results is not None:
        print(coexpression_results)

        # Updated save path
        save_path = r"C:\Users\orenh\PycharmProjects\pythonProject\coexpression_results.csv"
        coexpression_results.to_csv(save_path, index=False)
        print(f"Co-expression results saved to '{save_path}'.")


if __name__ == "__main__":
    main()