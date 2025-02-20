#!/bin/bash
# PediatricOpenTargets 2021
# Yuanchao Zhang
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
# copied from the run_in_ci.sh file at
# <https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/scripts/>
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

mkdir -p results

# if results directory already exists, remove gzip files
# if not, will print message saying files do not exist and then continue
rm -f results/*.gz || true

Rscript --vanilla '01-snv-frequencies.R'

# Convert JSON to JSON Lines (JSONL) format with jq
#
# jq: https://stedolan.github.io/jq/ JSONL: https://jsonlines.org/
#
# Each 01-tpm-summary-stats.R output json file is an array of objects, i.e.
# [{k11:v11,k12:v12,...}, {k21:v21,k22:v22,...}, ...].
#
# jq docs:
#
# --compact-output / -c: By default, jq pretty-prints JSON output. Using this
# option will result in more compact output by instead **putting each JSON
# object on a single line**.
#
# Note putting each JSON object on a single line follows JSONL
#
# Array/Object Value Iterator: .[]
#
# If you use the .[index] syntax, but omit the index entirely, it will **return
# all of the elements of an array**. Running .[] with the input [1,2,3] will
# produce the numbers as three separate results, rather than as a single array.
#
# You can also use this on an object, and it will return all the values of the
# object.
#
# Commands adapted from
#
# - https://stackoverflow.com/a/48711608/4638182
# - https://stackoverflow.com/a/66709708/4638182
echo 'Convert JSON files to JSONL files...'

jq --compact-output '.[]' \
  results/variant-level-snv-consensus-annotated-mut-freq.json \
  > results/variant-level-snv-consensus-annotated-mut-freq.jsonl

jq --compact-output '.[]' \
  results/gene-level-snv-consensus-annotated-mut-freq.json \
  > results/gene-level-snv-consensus-annotated-mut-freq.jsonl

echo 'Remove JSON files...'

rm results/variant-level-snv-consensus-annotated-mut-freq.json
rm results/gene-level-snv-consensus-annotated-mut-freq.json

# --no-name option stops the filename and timestamp from being stored in the
# output file. So rerun will have the same file.

echo 'gzip JSONL files...'

gzip --no-name results/variant-level-snv-consensus-annotated-mut-freq.jsonl
gzip --no-name results/gene-level-snv-consensus-annotated-mut-freq.jsonl

echo 'gzip results/variant-level-snv-consensus-annotated-mut-freq.tsv...'

gzip --no-name results/variant-level-snv-consensus-annotated-mut-freq.tsv

echo 'gzip results/gene-level-snv-consensus-annotated-mut-freq.tsv...'

gzip --no-name results/gene-level-snv-consensus-annotated-mut-freq.tsv

echo 'Done running run-snv-frequencies.sh'
