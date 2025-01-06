import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import argparse


def visualize_scrnaseq_heatmap(csv_file, num_cells=30):
    """
    Visualize a heatmap for scRNA-seq data from a CSV file.
    """
    print(f"Loading data from {csv_file}...")
    try:
        data = pd.read_csv(csv_file, index_col=0)
        print("Data loaded successfully.")
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return

    print(f"Creating heatmap for {num_cells} random cells...")
    try:
        cell_subset = random.sample(list(data.columns), min(num_cells, data.shape[1]))
        subset_data = data.loc[:, cell_subset]

        plt.figure(figsize=(12, 10))
        sns.heatmap(
            subset_data,
            cmap="RdBu_r",
            linewidths=0.5,
            linecolor="black",
            cbar_kws={"label": "Expression Level"},
        )
        plt.title("scRNA-seq Heatmap for Multiple Genes")
        plt.xlabel("Cells")
        plt.ylabel("Genes")
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"Error visualizing heatmap: {e}")


def main():
    parser = argparse.ArgumentParser(description="Visualize a heatmap from scRNA-seq data.")
    parser.add_argument(
        "--csv_file",
        type=str,
        default="filtered_expression_matrix.csv",
        help="Path to the input CSV file",
    )
    parser.add_argument(
        "--num_cells",
        type=int,
        default=30,
        help="Number of random cells to visualize",
    )

    args = parser.parse_args()
    visualize_scrnaseq_heatmap(args.csv_file, num_cells=args.num_cells)


if __name__ == "__main__":
    main()
