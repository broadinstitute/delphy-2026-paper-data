#!/usr/bin/bash

./run.py --rep b --sim const_100
./run.py --rep b --sim const_1000
./run.py --rep b --sim const_10000
./run.py --rep b --sim const_100000 --coal-cells 10000  # Mitigate discretization error in parallelized coal. prior
