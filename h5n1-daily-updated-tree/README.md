# Background

The `run.py` script in this folder is meant to run as a daily cron job, as detailed below.
It should run in a dedicated Python venv, created as below using the sibling `requirements.txt`.

# Running

Make a Python venv with the requirements:
```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Link to `mafft` and `delphy`, wherever they are:
```
ln -s `which mafft` ./mafft
ln -s `which delphy` ./delphy
```

Then run this script often (say, daily):
```
./run.py my@email.com    # The email address is for attaching to GenBank queries
```
If there is an update to the Andersen lab's `avian-influenza` repo, after under an hour,
there should be a complete Delphy run ready to publish under:
```
results/{commit_date}-{commit_hash}/us-h5n1-{commit_date}-{commit_hash}.dphy
```
(e.g., `results/2025-03-20-c4ea76da2/us-h5n1-2025-03-20-c4ea76da2.dphy`).
