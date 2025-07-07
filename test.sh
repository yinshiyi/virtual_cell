python prepare_vcc_submission.py \
  --pert_counts ./pert_counts_Validation.csv \
  --gene_names ./gene_names.csv \
  --training_adata ./adata_Training.h5ad \
  --output_h5ad ./example.h5ad
cell-eval prep -i ./example.h5ad --genes ./gene_names.csv
