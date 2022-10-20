# AlphaFold-Multimer Peptide-Receptor ranking
Identifying peptide-receptor interactions using AlphaFold-Multimer.

![alphafold_receptors_ranking](https://user-images.githubusercontent.com/56223326/197050884-daa0c071-42b1-4d1c-bbb3-4a2088fbf1cb.png)

## Prerequisites

- Installation of AlphaFold 2.2.0 - we used the docker-free version provided in https://github.com/kalininalab/alphafold_non_docker
- As we split MSA generation from prediction, copy `af_scripts/run_alphafold_msaonly.py` into the root directory of alphafold (that contains `run_alphafold.py`). This script only runs the data generation pipeline and omits the neural network execution.


## Run AlphaFold
- Execute `af_scripts/precompute_msas.py` to make all MSAs. The working directory needs to be the alphafold root dir. To change the data or run parameters, modify the variables on lines 14 to 21.
- Execute `af_scripts/predict_from_precomputed.py` to predict all pairwise complexes. Modify the variables on lines 18 to 26 if you changed the data or msa directories. The script is meant to be executed on a GPU node and spawns multiple AlphaFold processes in parallel. Modify `GPU_AVAILABLE` starting from line 34 to match your GPU setup (default assumes 8 GPUs available)


## Rank receptors

- The function to extract the metrics from a single alphafold result is defined in `qc_metrics.py`. In `benchmark.ipynb`, we apply this function to all results, aggregate the metrics and rank the receptors. The notebook produces the results presented in the manuscript.
