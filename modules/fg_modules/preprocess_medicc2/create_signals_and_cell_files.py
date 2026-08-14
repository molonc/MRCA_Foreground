import sys
import os
import argparse
import scgenome
import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description='Import HMMcopy data using scgenome')
parser.add_argument('hscn', help='hscn input data # single-cell DNA seq data with allele')
parser.add_argument('metrics', help='metrics data # must contain cell_id,sample_type,multiplier')
parser.add_argument('--cell_name_output', required=False, default='cell_list.txt', help='name of cell output file ".txt" (default: cell_list.txt)')
parser.add_argument('--ann_data', required=False, default='anndata.h5ad', help='name of ann data output file ".h5ad" (default: anndata.h5ad)')
parser.add_argument('--number_of_cells_returned', required=False, type=int, help='number of cells to include in the cell list')

## help
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

# Parse arguments, reads in terminal to make list called args
args = parser.parse_args()
# Print the arguments for confirmation
print(f"hscn: {args.hscn}") 
print(f"metrics: {args.metrics}")
print(f"--cell_name_output: {args.cell_name_output}")
print(f"--ann_data: {args.ann_data}")
print(f"--number_of_cells_returned: {args.number_of_cells_returned}")

metrics = pd.read_csv(args.metrics)
hscn = pd.read_csv(args.hscn)
hscn['chr'] = hscn['chr'].astype(str)
hscn['start'] = hscn['start'].astype(int)
hscn['end'] = hscn['end'].astype(int)
if (hscn.loc[:,'state_AS_phased'] == "NA|NA").any():
    print(f"{(hscn.loc[:,'state_AS_phased']=='NA|NA').sum()} Allele Phase Assignments did not work out of {hscn.shape[0]}")
    print(f"These alleles are from cells: {hscn[hscn.loc[:,'state_AS_phased']=='NA|NA'].cell_id.unique()}")
    print(hscn[hscn.loc[:,'state_AS_phased']=="NA|NA"])
    print("Filtering out NA Phased Alleles")
    hscn = hscn[hscn.loc[:,'state_AS_phased'] != "NA|NA"]
    metrics = metrics[metrics.loc[:,'cell_id'].isin(hscn['cell_id'].unique())]
    
    
    
# list the column names
# list(hscn.columns.values)

print(hscn['state_AS_phased'].unique())

# hscn = hscn[hscn['chr']!="X"]

# Rename from "A" to "Maj" and from "B" to "Min"
hscn = hscn.rename(columns={'A': 'Maj', 'B': 'Min'})
hscn = hscn.astype({"start": int, "end": int})



adata = scgenome.pp.convert_dlp_signals(hscn, metrics)
# I have to remove sample_type column because it contains NAs and can't write
adata.obs = adata.obs.drop('sample_type', axis=1)
# Adding brief_cell_id column, needed at the next step
adata.obs["brief_cell_id"] = adata.obs.axes[0].to_list()
# Coerce some columns to strings
adata.obs['is_contaminated'] = adata.obs['is_contaminated'].astype(str)
adata.obs['is_s_phase'] = adata.obs['is_s_phase'].astype(bool)
adata.obs['brief_cell_id'] = adata.obs['brief_cell_id'].astype(str)
adata.obs['cell_call'] = adata.obs['cell_call'].astype(str)
adata.obs['experimental_condition'] = adata.obs['experimental_condition'].astype(str)
adata.obs['index_i5'] = adata.obs['index_i5'].astype(str)
adata.obs['index_i7'] = adata.obs['index_i7'].astype(str)
adata.obs['index_sequence'] = adata.obs['index_sequence'].astype(str)
adata.obs['is_control'] = adata.obs['is_control'].astype(str)
adata.obs['library_id'] = adata.obs['library_id'].astype(str)
adata.obs['primer_i5'] = adata.obs['primer_i5'].astype(str)
adata.obs['primer_i7'] = adata.obs['primer_i7'].astype(str)
adata.obs['sample_id'] = adata.obs['sample_id'].astype(str)

# Make sure pandas prints everything
pd.set_option('display.max_columns', None)     # Unlimited columns
print(adata.obs.dtypes)

adata.write(args.ann_data, compression="gzip")

## get the cell ids - a subset or all
# can select a subset of cells
if args.number_of_cells_returned:
    cell_name_output = adata.obs.axes[0].to_list()[0:args.number_of_cells_returned]
else:
    cell_name_output = adata.obs.axes[0].to_list()  # Default to all

with open(args.cell_name_output, "w") as outfile:
    outfile.write("\n".join(cell_name_output))
