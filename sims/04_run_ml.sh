#!/bin/bash

for n in 100 1000  # Only for runs where we output a FASTA file
do
    for sim in exp_${n} const_${n}
    do
        time ../iqtree2 -s ${sim}/inputs/${sim}.fasta -m HKY+FO --prefix ${sim}/ml_outputs/${sim}
        time treetime --tree ${sim}/ml_outputs/${sim}.treefile --dates ${sim}/ml_outputs/${sim}_dates.csv --aln ${sim}/inputs/${sim}.fasta --outdir ${sim}/ml_outputs --coalescent skyline --n-skyline 2 --stochastic-resolve
        
        # Tweak TreeTime output to not confuse baltic
        sed -i -e 's/ Tree /tree /g' ${sim}/ml_outputs/timetree.nexus
    done
done

