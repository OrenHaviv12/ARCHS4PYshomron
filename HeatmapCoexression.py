
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def generate_heatmap(csv_file, save_path):
    """
    Generate a heatmap from the co-expression results CSV file.

    Parameters:
        csv_file (str): Path to the CSV file containing co-expression data.
        save_path (str): Path to save the heatmap image.
    """
    try:
        # Load co-expression data
        coexpression_data = pd.read_csv(csv_file)

        # Ensure the required columns exist
        if 'Target Gene' not in coexpression_data.columns or 'Coexpression' not in coexpression_data.columns:
            raise ValueError("CSV file must contain 'Target Gene' and 'Coexpression' columns.")

        # Pivot data to prepare for heatmap
        heatmap_data = coexpression_data.pivot_table(values='Coexpression', index='Target Gene')

        # Generate the heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, annot=True, cmap="RdBu_r", center=0, fmt=".2f",
                    cbar_kws={"label": "Coexpression Coeff."})
        plt.title("Gene Co-expression Heatmap")
        plt.xlabel("Genes")
        plt.ylabel("Target Gene")

        # Save the heatmap as an image
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.show()
        print(f"Heatmap successfully saved to: {save_path}")
    except Exception as e:
        print(f"Error generating heatmap: {e}")


if __name__ == "__main__":
    # Specify paths
    csv_file = r"C:\Users\orenh\PycharmProjects\pythonProject\coexpression_results.csv"
    heatmap_save_path = r"C:\Users\orenh\PycharmProjects\pythonProject\heatmapcoexpression.png"

    # Generate the heatmap
    generate_heatmap(csv_file, heatmap_save_path)
