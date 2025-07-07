import numpy as np
import pandas as pd
import polars as pl
import anndata as ad
from numpy.typing import NDArray
import argparse

def random_predictor(
    pert_names: NDArray[np.str_],
    cell_counts: NDArray[np.int64],
    gene_names: NDArray[np.str_],
    max_count: int | float = 1e4,
    log1p: bool = True,
) -> ad.AnnData:
    """Generate a random AnnData with the expected number of cells per perturbation."""
    total_cells = cell_counts.sum()
    matrix = np.random.randint(0, int(max_count), size=(total_cells, gene_names.size))
    if log1p:
        matrix = np.log1p(matrix)
    return ad.AnnData(
        X=matrix,
        obs=pd.DataFrame({"target_gene": np.repeat(pert_names, cell_counts)},
                         index=np.arange(total_cells).astype(str)),
        var=pd.DataFrame(index=gene_names)
    )

def main(args):
    # Load expected perturbation counts
    pert_counts = pl.read_csv(args.pert_counts).to_pandas()
    gene_names = pl.read_csv(args.gene_names, has_header=False).to_numpy().flatten()

    # Generate predictions
    adata = random_predictor(
        pert_names=pert_counts["target_gene"].to_numpy(),
        cell_counts=pert_counts["n_cells"].to_numpy(),
        gene_names=gene_names,
    )

    # Optionally add non-targeting controls
    if args.training_adata:
        tr_adata = ad.read_h5ad(args.training_adata)
        ntc_adata = tr_adata[tr_adata.obs["target_gene"] == "non-targeting"]
        if "non-targeting" not in adata.obs["target_gene"].unique():
            assert np.all(adata.var_names.values == ntc_adata.var_names.values), (
                "Gene names are out of order or unequal"
            )
            adata = ad.concat([adata, ntc_adata])

    # Save output
    adata.write_h5ad(args.output_h5ad)
    print(f"Saved prediction AnnData to: {args.output_h5ad}")
    print("Run the following command to prepare for scoring:")
    print(f"cell-eval prep -i {args.output_h5ad} --genes {args.gene_names}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare VCC submission AnnData")
    parser.add_argument("--pert_counts", type=str, required=True,
                        help="Path to pert_counts_Validation.csv")
    parser.add_argument("--gene_names", type=str, required=True,
                        help="Path to gene_names.csv")
    parser.add_argument("--training_adata", type=str, default=None,
                        help="Path to training AnnData (optional, for non-targeting control)")
    parser.add_argument("--output_h5ad", type=str, default="example.h5ad",
                        help="Output path for generated .h5ad file")

    args = parser.parse_args()
    main(args)
