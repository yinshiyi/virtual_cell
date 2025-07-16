data=/Users/yins/Desktop/vcc_data/
output=20250710
python main.py \
  --pert_counts custom.csv \
  --gene_names $data/gene_names.csv \
  --output_h5ad $data/$output.h5ad