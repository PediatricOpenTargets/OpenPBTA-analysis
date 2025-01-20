#!/bin/bash
# Module author: Aryan Neupane
# Date: 2025-01
# Purpose: Run medulloblastoma immune composition analysis for both deconvolution methods

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Run medulloblastoma analysis for xCell
echo "Running medulloblastoma analysis for xCell..."
Rscript --vanilla 02-medulloblastoma-analysis.R \
--deconv_method 'xcell' \
--clin_file '../../data/histologies.tsv' \
--results_dir 'results' \
--output_dir 'results'

## Run medulloblastoma analysis for quanTIseq
#echo "Running medulloblastoma analysis for quanTIseq..."
#Rscript --vanilla 02-medulloblastoma-analysis.R \
#--deconv_method 'quantiseq' \
#--clin_file '../../data/histologies.tsv' \
#--results_dir 'results' \
#--output_dir 'results'