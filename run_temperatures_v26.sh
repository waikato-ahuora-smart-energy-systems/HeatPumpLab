#!/bin/bash
# Four v26 mixtures (compressor model 5, motor 0.93, 1 wt% minimum) at 150-190 C, 1,000 evaluations each.
cd "$(dirname "$0")"
for T in 150 160 170 180 190; do
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python main_RHPfluid_sweep.py \
    --temperature $T --output results/RHPfixed${T}_v26_comp5_min1 --slots 4 --top 0 \
    --comp-model 5 --motor-eff 0.93 \
    --mixtures Propane/Pentane Propane/Isopentane Propane/Butane/Pentane Propane/Butane/Isopentane
done
echo ALL TEMPERATURES DONE
