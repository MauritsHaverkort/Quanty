#!/bin/bash
set -e

echo "Generating matrix-jobs.yml..."

# Define test folders and their corresponding CI stages
declare -A folders
folders["TestScripts/Math"]="test_math"
folders["TestScripts/Matrix"]="test_matrix"
folders["TestScripts/Responsefunction"]="test_responsefunction"
folders["TestScripts/TightBindingDef"]="test_tightbinding"
folders["TestScripts/Plot"]="test_plot"
folders["TestScripts"]="test"

echo "" > matrix-jobs.yml

for folder in "${!folders[@]}"; do
  stage="${folders[$folder]}"
  for filepath in "$folder"/*.Quanty; do
    [ -e "$filepath" ] || continue
    testname=$(basename "$filepath" .Quanty)
    # Normalize folder name for job naming
    foldername=$(basename "$folder")
    jobname="test-${foldername,,}-$testname"
    jobname="${jobname//_/-}"
    displayname="$(basename "$folder"): $testname"

    cat >> matrix-jobs.yml <<EOF
$jobname:
  stage: $stage
  script:
    - export OMP_NUM_THREADS=8
    - ./Quanty.exe "$filepath"
  dependencies:
    - build-job

EOF
  done
done

