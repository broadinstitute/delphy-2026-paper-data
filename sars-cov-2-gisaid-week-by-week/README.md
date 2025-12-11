This directory would include sequences and metadata from GISAID, which have been removed to avoid making them public
under [GISAID's Data Access Agreement](https://gisaid.org/terms-of-use/).  Original files may be obtained by researchers
with GISAID access by contacting `pvarilly@broadinstitute.org`.

Specifically, the file `20200331.fasta.xz` contains all (unaligned) GISAID sequences collected on or before 2020-03-31,
while `metadata_20200331.tsv` contains metadata for them with the following columns (after the first header row):

- Virus name
- Accession ID
- Collection date (YYYY-MM-DD or YYYY-MM)
- Location
- Additional location information
- Sequence length
- Host
- Submission date (YYYY-MM-DD)

With a suitably populated FASTA and metadata file, running the scripts in this folder (00_prepare_runs.sh, 01_run_by_collection_date.sh, 02_run_by_submission_date.sh) should produce the results reported in the Delphy paper.

# Preparing the AWS machine for these runs

- Launch an Ubuntu 22.04 LTS x86-64 instance of type `c5a.4xlarge` (16 vCPUs & 32 GB memory) with 8GB gp2 storage
- SSH into the machine, e.g.
```
  ssh -i "~/.ssh/2023-01-29-aws-vs.pem" ubuntu@ec2-3-78-245-33.eu-central-1.compute.amazonaws.com
```
- Upgrade Ubuntu packages
```
    sudo apt update
    sudo apt upgrade  # instance may need to be restarted; do it (`sudo shutdown -r now`) and log back in
```
- Upload `delphy`, inputs and scripts, e.g.
```
  ssh ubuntu@AWS-INSTANCE-IP mkdir sars-cov-2-gisaid-week-by-week
  scp -r ../delphy ubuntu@AWS-INSTANCE-IP:.
  scp -r inputs_by_submission_date/ 01_run_by_submission_date.sh run.py ubuntu@AWS-INSTANCE-IP:sars-cov-2-gisaid-week-by-week
```
Now run the driver script and copy back the output to analyze it

The runs were performed twice, and the two independent outputs were saved to the _a and _b folders