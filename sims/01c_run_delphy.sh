#!/usr/bin/bash

# Calibrate effect of discretization error in parallelized coal. prior
# Not really timing these, so it's ok for them to run in parallel to better use the 96-vCPU AWS instances
./run.py --rep c_5000bins --sim exp_100000 --coal-cells 5000 &
./run.py --rep c_2500bins --sim exp_100000 --coal-cells 2500 &
./run.py --rep c_1250bins --sim exp_100000 --coal-cells 1250 &
./run.py --rep c_625bins  --sim exp_100000 --coal-cells 625 &

