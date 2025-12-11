#!/usr/bin/bash

# delphy_mcc has a hard-coded 30% burn-in

for SIM in exp_100 \
           exp_1000 \
           exp_10000 \
           exp_100000 \
           const_100 \
           const_1000 \
           const_10000 \
           const_100000
do
    ../delphy_mcc ${SIM}/delphy_outputs_a/${SIM}.trees ${SIM}/delphy_outputs_a/${SIM}.mcc
    ../delphy_mcc ${SIM}/delphy_outputs_b/${SIM}.trees ${SIM}/delphy_outputs_b/${SIM}.mcc
done

for COAL_BINS in 625 1250 2500 5000
do
    SIM=exp_100000
    ../delphy_mcc ${SIM}/delphy_outputs_c_${COAL_BINS}bins/${SIM}.trees ${SIM}/delphy_outputs_c_${COAL_BINS}bins/${SIM}.mcc
done
