#!/bin/bash

for d in 0[2-9]* 1* 2*
do
    echo $d
    ./plot_compact_summary.py $d
done
