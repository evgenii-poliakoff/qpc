#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

echo INFO Cleaning up build artifacts in the repository

rm -rf "${script_dir}"/../conda-env

eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda clean -f --packages
