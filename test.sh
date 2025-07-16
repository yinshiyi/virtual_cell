# pyenv virtualenv 3.10.13 vcc-env
# pyenv activate vcc-env
# pip install -r requirements.txt

data1=/Users/yins/Desktop/vcc_data/
data=vcc_data/
output=20250710
python main_og.py \
  --pert_counts $data/pert_counts_Validation.csv \
  --gene_names $data/gene_names.csv \
  --training_adata $data/adata_Training.h5ad \
  --output_h5ad $data/$output.h5ad
cell-eval prep -i $data/$output.h5ad --genes $data/gene_names.csv
