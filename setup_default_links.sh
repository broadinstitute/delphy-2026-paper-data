#!/usr/bin/bash

[ -f mafft ] || ln -s "`which mafft`"
[ -f delphy ] || ln -s "${HOME}/now/delphy/build/release/delphy"
[ -f delphy_mcc ] || ln -s "${HOME}/now/delphy/build/release/delphy_mcc"
[ -f loganalyser277 ] || ln -s "${HOME}/github/CompEvol/beast2.7.7/bin/loganalyser" loganalyser277
[ -f sapling ] || ln -s "${HOME}/now/sapling/build/release/sapling"
[ -f iqtree2 ] || ln -s "${HOME}/tools/iqtree-2.3.6-Linux-intel/bin/iqtree2"
[ -f beastX1050 ] || ln -s "${HOME}/tools/BEASTv10.5.0/bin/beast" beastX1050
[ -f beast277 ] || ln -s "${HOME}/github/CompEvol/beast2.7.7/bin/beast" beast277
[ -f calc-tree-ess ] || ln -s "`pwd`/tree_ess/target/release/calc-tree-ess" calc-tree-ess
[ -f createMapleFile.py ] || ln -s "${HOME}/github/NicolaDM/MAPLE/createMapleFile.py" createMapleFile.py

#[ -f h5n1-daily-updated-tree/mafft ] || ln -s "`pwd`/mafft" h5n1-daily-updated-tree/mafft
#[ -f h5n1-daily-updated-tree/delphy ] || ln -s "`pwd`/delphy" h5n1-daily-updated-tree/delphy
