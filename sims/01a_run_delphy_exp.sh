#!/usr/bin/bash

./run.py --rep a --sim exp_100
./run.py --rep a --sim exp_1000
./run.py --rep a --sim exp_10000
./run.py --rep a --sim exp_100000 --coal-cells 10000  # Mitigate discretization error in parallelized coal. prior
