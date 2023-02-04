#!/bin/sh

set -e
set -o pipefail

printf "Start generating pre-release files...\n\n"

# Set locations for s3 bucket that contains release files
URL="s3://d3b-openaccess-us-east-1-prd-pbta/open-targets"
RELEASE="analysis-files-pre-release"

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# If RUN_LOCAL is used, the time-intensive steps are skipped because they cannot
# be run on a local computer -- the idea is that setting RUN_LOCAL=1 will allow for
# local testing running/testing of all other steps
RUN_LOCAL=${RUN_LOCAL:-0}

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -
  
analyses_dir="$BASEDIR/analyses"
data_dir="$BASEDIR/data"
scratch_dir="$BASEDIR/scratch"

# Compile all the files that need to be included in the release in one place
# in the scratch directory
release_dir="${scratch_dir}/analysis-files-pre-release"
mkdir -p ${release_dir}

# Create the independent sample list using the *base* histology file (i.e. - histologies-base.tsv)
echo "Create independent sample list for fusion filtering module"
cd ${analyses_dir}/independent-samples
OPENPBTA_BASE_SUBTYPING=1 bash run-independent-samples.sh

# Run fusion filtering
echo "Create fusion filtered list"
cd ${analyses_dir}/fusion_filtering
OPENPBTA_BASE_SUBTYPING=1 bash run_fusion_merged.sh

# Copy over fusions lists
cp ${analyses_dir}/fusion_filtering/results/fusion-putative-oncogenic.tsv ${release_dir}

# Run modules that cannot be run locally due to memory requirements
if [ "$RUN_LOCAL" -lt "1" ]; then
  
  # Run GISTIC step -- only the part that generates ZIP file
  echo "Run GISTIC"
  # Run a step that subs ploidy for NA to allow GISTIC to run
  Rscript ${analyses_dir}/run-gistic/scripts/prepare_seg_for_gistic.R \
  --in_consensus ${data_dir}/cnv-consensus.seg.gz \
  --out_consensus ${analyses_dir}/run-gistic/results/cnv-consensus-gistic-only.seg.gz \
  --histology ${data_dir}/histologies-base.tsv
  
  # This will use the file that just got generated above
  bash ${analyses_dir}/run-gistic/scripts/run-gistic-opentargets.sh
  
  # Copy over GISTIC
  cp ${analyses_dir}/run-gistic/results/cnv-consensus-gistic.zip ${release_dir}
  cp ${analyses_dir}/run-gistic/results/cnv-consensus-gistic-only.seg.gz ${release_dir}
  
  # Run step that generates "most focal CN" files (annotation) using the *BASE* histology file
  echo "Run focal CN file preparation"
  cd ${analyses_dir}/focal-cn-file-preparation
  RUN_ORIGINAL=1 OPENPBTA_BASE_SUBTYPING=1 bash run-prepare-cn.sh
  
  ## Copy over focal CN
  cp ${analyses_dir}/focal-cn-file-preparation/results/consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz ${release_dir}
  cp ${analyses_dir}/focal-cn-file-preparation/results/consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz ${release_dir}
  cp ${analyses_dir}/focal-cn-file-preparation/results/consensus_wgs_plus_cnvkit_wxs.tsv.gz ${release_dir}

  # Copy over the consensus with status file
  cp ${analyses_dir}/focal-cn-file-preparation/results/consensus_seg_with_status.tsv ${release_dir}
  cp ${analyses_dir}/focal-cn-file-preparation/results/cnvkit_with_status.tsv ${release_dir}

fi

# Run fusion summary
echo "Run fusion summary for subtypes"
cd ${analyses_dir}/fusion-summary
bash run-new-analysis.sh

# Copy over fusion summary
cp ${analyses_dir}/fusion-summary/results/* ${release_dir}

# Run TMB
echo "Create TMB results"
cd ${analyses_dir}/tmb-calculation
bash run_tmb_calculation.sh

# Copy over TMB results
cp ${analyses_dir}/tmb-calculation/results/snv-mutation-tmb-coding.tsv ${release_dir}
cp ${analyses_dir}/tmb-calculation/results/snv-mutation-tmb-all.tsv ${release_dir}

# Create an md5sum file for all the files in the directories where the analysis
# files are compiled
cd ${release_dir}
# Remove old md5sum release file if it exists
rm -f analysis_files_release_md5sum.txt
# Create a new md5sum release file
md5sum * > analysis_files_release_md5sum.txt

# Upload all release files s3 bucket in their respective folders
aws s3 cp ${release_dir}/ $URL/$RELEASE/ --recursive

printf "\nDone generating pre-release files...\n\n"
