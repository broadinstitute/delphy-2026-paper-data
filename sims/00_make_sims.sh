#!/usr/bin/bash

for n in 100 1000 10000 100000
do
    dir=exp_${n}
    mkdir -p ${dir}/ground_truth
    mkdir -p ${dir}/inputs
    mkdir -p ${dir}/ml_outputs
    
    if (( n <= 1000 )); then
        fasta_out="--out-fasta ${dir}/inputs/exp_${n}.fasta"
    else
        fasta_out=
    fi
    ../sapling \
        --seed $n \
        -n $n \
        --exp-pop-n0 6.0 \
        --exp-pop-g 10 \
        --t0 2024-07-31 \
        --min-tip-t 2024-01-01 \
        --max-tip-t 2024-07-31 \
        --hky-kappa 5.0 \
        --hky-pi-A 0.3 \
        --hky-pi-C 0.18 \
        --hky-pi-G 0.20 \
        --hky-pi-T 0.32 \
        --mu 0.001 \
        -L 30000 \
        --out-info ${dir}/ground_truth/exp_${n}_info.json \
        --out-newick ${dir}/ground_truth/exp_${n}_tree.nwk \
        --out-nexus ${dir}/ground_truth/exp_${n}_tree.nexus \
        ${fasta_out} \
        --out-maple ${dir}/inputs/exp_${n}.maple
done

for n in 100 1000 10000 100000
do
    dir=const_${n}
    mkdir -p ${dir}/ground_truth
    mkdir -p ${dir}/inputs
    mkdir -p ${dir}/ml_outputs

    if (( n <= 1000 )); then
        fasta_out="--out-fasta ${dir}/inputs/const_${n}.fasta"
    else
        fasta_out=
    fi
    ../sapling \
        --seed $n \
        -n $n \
        --const-pop-n0 2.0 \
        --t0 2024-07-31 \
        --min-tip-t 2024-01-01 \
        --max-tip-t 2024-07-31 \
        --hky-kappa 5.0 \
        --hky-pi-A 0.3 \
        --hky-pi-C 0.18 \
        --hky-pi-G 0.20 \
        --hky-pi-T 0.32 \
        --mu 0.001 \
        -L 30000 \
        --out-info ${dir}/ground_truth/const_${n}_info.json \
        --out-newick ${dir}/ground_truth/const_${n}_tree.nwk \
        --out-nexus ${dir}/ground_truth/const_${n}_tree.nexus \
        ${fasta_out} \
        --out-maple ${dir}/inputs/const_${n}.maple
done
