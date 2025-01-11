import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


def generate_clustermap_and_matrix(csv_file, clustermap_save_path, matrix_save_path, cmap="Spectral"):
    """
    Generate a clustermap from the co-expression results CSV file and save the correlation matrix.

    Parameters:
        csv_file (str): Path to the CSV file containing co-expression data.
        clustermap_save_path (str): Path to save the clustermap image.
        matrix_save_path (str): Path to save the correlation matrix CSV.
        cmap (str): Colormap for the clustermap (default: "Spectral").
    """
    try:
        # Check if file exists
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"The file {csv_file} does not exist. Please verify the path.")

        # Load the co-expression data
        coexpression_data = pd.read_csv(csv_file)

        # Debug: Print the first few rows and column names for verification
        print("Loaded Co-expression Data:")
        print(coexpression_data.head())
        print(f"Columns: {list(coexpression_data.columns)}")

        # Ensure required columns exist
        required_columns = {'Primary Gene', 'Target Gene', 'Coexpression'}
        if not required_columns.issubset(coexpression_data.columns):
            raise ValueError(f"CSV file must contain columns: {required_columns}")

        # Pivot the data to prepare for heatmap/clustermap
        heatmap_data = coexpression_data.pivot_table(
            values='Coexpression',
            index='Target Gene',
            columns='Primary Gene'
        )

        # Save the correlation matrix as a CSV
        heatmap_data.to_csv(matrix_save_path)
        print(f"Correlation matrix successfully saved to: {matrix_save_path}")

        # Generate the clustermap with hierarchical clustering
        sns.clustermap(
            heatmap_data,
            cmap=cmap,
            figsize=(14, 12),
            annot=False,
            cbar_kws={"label": "Coexpression Coefficient"},
            center=0,
            xticklabels=True,
            yticklabels=True,
            linewidths=0.5,
            dendrogram_ratio=(0.2, 0.2),  # Adjust the size of dendrograms
        )
        plt.title("Gene Co-expression Clustermap", fontsize=16)
        plt.savefig(clustermap_save_path, dpi=300, bbox_inches="tight")
        plt.show()
        print(f"Clustermap successfully saved to: {clustermap_save_path}")

    except Exception as e:
        print(f"Error generating clustermap: {e}")


if __name__ == "__main__":
    # Specify the path to the co-expression CSV file
    csv_file = r"C:\Users\orenh\PycharmProjects\pythonProject5\coexpression_results.csv"  # Adjust the path as needed

    # Specify the path to save the clustermap image
    clustermap_save_path = r"C:\Users\orenh\PycharmProjects\pythonProject5\clustermap_coexpression_spectral.png"

    # Specify the path to save the correlation matrix
    matrix_save_path = r"C:\Users\orenh\PycharmProjects\pythonProject5\correlation_matrix_spectral.csv"

    # Generate the clustermap and save the matrix
    generate_clustermap_and_matrix(csv_file, clustermap_save_path, matrix_save_path, cmap="Spectral")
