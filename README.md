# Delphy paper data

Datasets, scripts for analyzing them, and scripts to generate figures.

Although these scripts were last successfully run in late 2025 to prepare, execute and
analyze the Delphy and BEAST2 runs, there is no guarantee that they will run perfectly in
the future and/or on your machine.  We do not intend to maintain these scripts moving
forward.  Instead, the aim of this repo is to be living documentation of all the details
of the runs and plots.

Some files in this repo (notably `*.trees` and `*.dphy` run results file) were to big too
push to GitHub.  Except for clones of public GitHub repos, the full contents are available
at Zenodo
[![here](https://zenodo.org/badge/DOI/10.5281/zenodo.17899935.svg)](https://doi.org/10.5281/zenodo.17899935).
See below for URLs of these public GitHub repos.

This repository supercedes the [analogous
repository](https://github.com/broadinstitute/delphy-2025-paper-data) that accompanied
[our preprint](https://www.biorxiv.org/content/10.1101/2025.03.25.645253v1).

## Usage

Create a Python virtualenv and install the associated requirements:
```
python3 -m venv delphy-env
source ./delphy-env/bin/activate
pip install -r requirements.txt
```

Compile the associated `calc-tree-ess` analysis program (you need a Rust installation):
```
cd tree_ess
cargo build --release
cd ..
```

Make symbolic links to the helper binaries you want to use here (you can run `setup_default_links.sh` to look up binaries via `which`, but you'll most likely need to adapt these steps manually):
```
beast277 (bin/beast inside the BEAST2 2.7.7 distribution; see below)
beastX1050 (bin/beast inside the BEAST X 10.5.0 distribution; see below)
delphy
delphy_mcc
iqtree2
loganalyser277 (bin/loganalyser in BEAST2 2.7.7 distribution; see below)
mafft
sapling (https://github.com/broadinstitute/sapling)
calc-tree-ess
```

Then enter each of the dataset directory and run the numbered scripts in order.  Warning: the BEAST runs can take a very long time to complete.

In the present repo, all the files created by these runs and analyses have been uploaded.  If you want to regenerate them from scratch, you should delete every *directory* inside each of the dataset directory (e.g., in `sars-cov-2-lemieux`, delete things like `delphy_outputs_a` and `raw`, but not `00_prepare_runs.py` nor `sample_ids.csv`).

The final plots that were composed into the paper figures are in the `plots` directory of each dataset.


## Software versions used

- Delphy Version 1.1.4 (build 2044, commit a50e378)
- MAFFT Version 7.505 (2022/Apr/10)
- BEAST2 v2.7.7
- BEAGLE commit `6480ad3` (Mon Sep 15 2025)
- BEAST X 10.5.0
- Sapling Version 0.1.1 (build 2, commit a0b9da1)
- IQ-TREE 2.3.6
- TreeTime 0.11.4
- calc-tree-ess 0.1.0 (in this repository)

# Preparing the AWS machine for the benchmarks

- Launch an Ubuntu 24.04 LTS x86-64 instance of type `c7a.2xlarge` (8 vCPUs & 16GB memory) with 8GB gp3 storage
  * For the `ebola-dudas-2017`, `h5n1-andersen-2025` and `sars-cov-2-gisaid-week-by-week` dataset, launch an instance of type `c7a.4xlarge` (16 vCPUs & 32GB memory)
  * For the `sims` dataset, launch an instance of type `c7a.24xlarge` (96 vCPUs & 192GB memory)
- SSH into the machine, e.g.
```
  ssh -i "~/.ssh/2023-01-29-aws-vs.pem" ubuntu@ec2-3-78-245-33.eu-central-1.compute.amazonaws.com
```
- Upgrade Ubuntu packages and restart
```
    sudo apt update -y
    sudo apt upgrade -y
    sudo shutdown -r now  # will disconnect; log back in when restart complete
```
- Install latest available Java LTS release (25 as of this writing) - only used by BEAST X; BEAST 2 releases includes its own JRE (17.0.3 as of this writing)
```
    sudo apt install -y openjdk-25-jdk
```
- Check Java works and print version:
```
    java -version

    > openjdk version "25" 2025-09-16
    > OpenJDK Runtime Environment (build 25+36-Ubuntu-124.04.2)
    > OpenJDK 64-Bit Server VM (build 25+36-Ubuntu-124.04.2, mixed mode, sharing)
```
- Download, unpack BEAST X and make a symbolic link to it
```
    wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0/BEAST_X_v10.5.0.tgz
    tar -xzvf BEAST_X_v10.5.0.tgz
    ln -s BEASTv10.5.0/bin/beast beastX1050
```
- Test it
```
    ./beastX1050 -version
    >
    >  BEAST v10.5.0, 2002-2025
    >  ...
```
- Download, unpack BEAST 2 and make a symbolic link to it
```
    wget https://github.com/CompEvol/beast2/releases/download/v2.7.7/BEAST.v2.7.7.Linux.x86.tgz
    tar -xvzf BEAST.v2.7.7.Linux.x86.tgz
    mv beast beast2.7.7
    ln -s beast2.7.7/bin/beast beast277
```
- Test it
```
    ./beast277 -version
    > BEAST v2.7.7
    > ---
    > BEAST.base v2.7.7
    > BEAST.app v2.7.7
    > ---
    > Java version 17.0.3
```
- Build and install BEAGLE from source (following these instructions: https://github.com/beagle-dev/beagle-lib/wiki/LinuxInstallInstructions)
```
    # Don't download JDK 11, already got JDK 25 above
    sudo apt-get install -y cmake build-essential autoconf automake libtool git pkg-config # openjdk-11-jdk
    export JAVA_HOME=/usr/lib/jvm/java-25-openjdk-amd64/  # Need this for CMake to find JDK libs below
    echo 'export JAVA_HOME=/usr/lib/jvm/java-25-openjdk-amd64/' >> ~/.bashrc # Also add that same line to ~/.bashrc

    git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    git checkout 6480ad3   # Pin the version of Beagle we use for benchmarking
    mkdir build
    cd build
    cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DCMAKE_INSTALL_PREFIX:PATH=$HOME ..
    make -j install

    export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH  # So that BEAST finds BEAGLE
    echo 'LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH' >> ~/.bashrc # Also add that same line to ~/.bashrc
    
    make test  # Should work
    cd ../..  # Back to home
```
- Ensure that BEASTX finds BEAGLE
```
    ./beastX1050 -beagle_info
    > ...
    > --- BEAGLE RESOURCES ---
    > 
    >0 : CPU (x86_64)
    > Flags: PRECISION_SINGLE PRECISION_DOUBLE COMPUTATION_SYNCH EIGEN_REAL EIGEN_COMPLEX SCALING_MANUAL SCALING_AUTO SCALING_ALWAYS SCALERS_RAW SCALERS_LOG VECTOR_SSE VECTOR_NONE THREADING_CPP THREADING_NONE PROCESSOR_CPU FRAMEWORK_CPU PREORDER_TRANSPOSE_MANUAL PREORDER_TRANSPOSE_AUTO PREORDER_TRANSPOSE_LOW_MEMORY
```
- Ensure that BEAST2 finds BEAGLE
```
    ./beast277 -beagle_info
    > ...
    > --- BEAGLE RESOURCES ---
    > 
    >0 : CPU (x86_64)
    > Flags: PRECISION_SINGLE PRECISION_DOUBLE COMPUTATION_SYNCH EIGEN_REAL EIGEN_COMPLEX SCALING_MANUAL SCALING_AUTO SCALING_ALWAYS SCALERS_RAW SCALERS_LOG VECTOR_SSE VECTOR_NONE THREADING_CPP THREADING_NONE PROCESSOR_CPU FRAMEWORK_CPU

```
- Now upload any BEAST2.7.7 XML file and run it with the following command (`-threads -1` uses as many threads as there are CPUs, `-beagle` enforces the use of BEAGLE)
```
    cd path/containing/beast2/xml
    time ~/beast277 -threads -1 -beagle <beast2_input.xml>
``` 
- And any BEAST X v10.5.0 XML file is run like this (BEAST X uses as many threads as there are CPUs by default, `-beagle` enforces the use of BEAGLE)
```
    cd path/containing/beastX/xml
    time ~/beastX1050 -beagle <beastX_input.xml>
``` 

# Benchmark datasets

## ebola-dudas-2017

- 1610 Ebola sequences (18,996 sites) covering the 2014-16 West Africa outbreak
- Paper: Dudas et al 2017 - https://dx.doi.org/10.7488/ds/1711
- Data: https://github.com/ebov/space-time
- Run with Skygrid.  Cutoff date = 22-Oct-2013.  24 intervals.
- For metadata, use 'Geo' by default

## ebola-gire-2014

- 81 sequences Ebola sequences (18,959 sites) covering the first few months of the 2014-16 West Africa outbreak
- Paper: Gire et al 2014 - https://dx.doi.org/10.1126/science.1259657
- Data: (SI of paper)
- Run with all default parameters
- For metadata, use 'Country' by default

## h3n2-rambaut-2008

- 165 H3N2 HA sequences (1698 sites) for 3 waves of influenza A in New York covering 2000-2003
- Paper: Rambaut et al 2008, https://dx.doi.org/10.1038/nature06945
- Data: https://beast.community/workshop_influenza_phylodynamics  (BEAST X Skygrid tutorial)
- Run with Skygrid.  Cutoff date = 23-Dec-1997.  72 intervals.
- No metadata

## h5n1-andersen-2025

- 3339 H5N1 genotype B3.13 sequences (all segments concatenated - 13,136 sites) covering the 2024-25 US _cattle_ outbreak (excludes sequences found in non-cattle hosts).  Stops at 3 Aug 2025 !
- No specific paper
- Data: https://github.com/andersen-lab/avian-influenza at frozen at commit e756a15
- Run with Skygrid.  Cutoff date = 3-Oct-2023.  22 intervals.
- For metadata, use 'geo' by default

## mpox-otoole-2023

- 41 mpox sequences (197,209 sites) around the 2017-2020 Nigerian mpox outbreak that turned into the 2021-2022 worldwide mpox outbreak
- Paper: O'Toole et al 2023, https://dx.doi.org/10.1126/science.adg8116
- Data: https://github.com/hmpxv/apobec3
- Run with "Model APOBEC action"
- For metadata, use 'geo' by default

## mpox-parker-2025

- 177 mpox sequences (197,209 sites), almost all from Nigeria, mostly preceding the 2021-2022 worlwide mpox outbreak
- Paper: Parker et al 2025, https://doi.org/10.1038/s41586-025-09128-2
- Data: https://github.com/andersen-lab/Mpox_West_Africa
- Run with "Model APOBEC action"
- For metadata, use 'Region' by default

## sars-cov-2-lemieux

- 757 SARS-CoV-2 sequences (29,903 sites) from Massachussetts in the first months of the SARS-CoV-2 pandemic (2020-)
- Paper: Lemieux et al 2021, https://doi.org/10.1126/science.abe3261
- Data: GenBank, accession IDs obtained from authors
- Run with all-default parameters
- For metadata, use 'clade' by default

## zika-metsky-2017

- 174 Zika virus sequences (10,807 sites) from the Americas spanning the 2015-16 Zika outbreak
- Paper: Metsky et al 2017, https://dx.doi.org/10.1038/nature22402
- Data: SI from the paper
- Run with all-default parameters
- For metadata, use 'Geo' by default


# Omitted public GitHub repos

To keep the size of this data repository manageable, we do not include the contents of the
following GitHub repos.  However, running the associated benchmarks requires you to have
checked out the repo (or run `00_prepare_runs.py` in each directory to clone them automatically):

- `ebola-dudas-2017/raw/space-time`:
```
git clone https://github.com/ebov/space-time.git
cd space-time
git checkout 9db59a4
```

- `h5n1-andersen-2025/raw/avian-influenza`:
```
git clone https://github.com/andersen-lab/avian-influenza.git
cd avian-influenza
git checkout e756a15
```

- `h5n1-daily-updated-tree/avian-influenza`:
```
git clone https://github.com/andersen-lab/avian-influenza.git
```

- `mpox-otoole-2023`
```
git clone https://github.com/hmpxv/apobec3.git
```

- `mpox-parker-2025`
```
git clone https://github.com/andersen-lab/Mpox_West_Africa.git
```

Additionally, the raw inputs and outputs in `sars-cov-2-gisaid-week-by-week` cannot be shared publicly owing to [GISAID's Data Access Agreement](https://gisaid.org/terms-of-use/).  See [sars-cov-2-gisaid-week-by-week/README.md](sars-cov-2-gisaid-week-by-week/README.md).