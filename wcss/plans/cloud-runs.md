# Plan: Offload WCSS runs to Hetzner Cloud instances

## Goal

Three Python scripts to manage WCSS runs on Hetzner Cloud:

1. **`cloud_setup.py`** — Commission a Hetzner instance and upload
   everything needed to run a WCSS study.
2. **`cloud_run.py`** — Kick off a configurable subset of runs on a
   particular instance, then periodically report progress and rsync
   results back.
3. **`cloud_teardown.py`** — Delete a Hetzner instance and clean up
   local metadata.

## Prerequisites

- **`hcloud` CLI** installed locally (`apt install hcloud` or from
  GitHub releases).  Configured with a Hetzner API token
  (`hcloud context create delphy`).
- **SSH key** registered with Hetzner (`hcloud ssh-key create`),
  with the corresponding private key available locally.
- Delphy binary built locally (release build).
- Sim inputs already generated locally via `01_generate.py`.

## What gets transferred

Per study, the Makefile drives runs from `sims/`.  Each `sim_NNN/`
needs only `sim.maple` (64K) as input.  The Makefile and the `delphy`
binary (2.7M) are also needed.

Files **not** transferred to the instance:
- `sim.fasta`, `sim-COMPLETE.*`, `sim.nexus`, `sim.nwk` (only needed
  for analysis/plotting, not for Delphy runs)
- `sim_info.json` (only needed for analysis)

Files **rsynced back** from the instance:
- `delphy.log`, `delphy.trees`, `.done` — all outputs in each
  `sim_NNN/` directory.
- `delphy.dphy` — optionally excluded with `--no-dphy-files` (these
  can be ~100MB each).

## Script 1: `cloud_setup.py`

### Usage

```bash
# First study on a new instance:
./cloud_setup.py \
  --study 06_missing_data \
  --name wcss-06-a \
  --type cpx52 \
  --location nbg1 \
  --ssh-key hetzner \
  --ssh-key-file ~/.ssh/id_hetzner \
  --delphy ~/now/delphy/build/release/delphy

# Add another study to the same instance:
./cloud_setup.py \
  --study 07_missing_data_no_alpha \
  --name wcss-06-a
```

### What it does

Each instance gets its own `/root/wcss/<study>/sims/` directory, so
different studies can coexist on the same instance without cross-talk.
The Delphy binary is shared at `/root/wcss/delphy`.  Lock files are
per-study (`/root/wcss/<study>/make.lock`) to allow concurrent runs of
different studies.

To add multiple studies to one instance, call `cloud_setup.py`
multiple times with the same `--name`.  If the instance already exists
(detected by finding `.wcss-cloud/<name>.json`), server creation is
skipped and only the file upload for the new study is performed.  The
metadata file is updated to record the additional study.

1. **Create the instance** via `hcloud server create` (skipped if
   `.wcss-cloud/<name>.json` already exists — the IP is read from
   it instead):
   - Image: `ubuntu-24.04`
   - Server type from `--type` (default: `cpx52` — 12 vCPUs, 24GB
     RAM, ~EUR 0.09/hr)
   - Location from `--location` (default: `nbg1` — Nuremberg)
   - SSH key from `--ssh-key` (name registered with Hetzner)
   - Name from `--name`
   - `--without-ipv4` (IPv6 only — avoids the IPv4 surcharge)
   - Local SSH private key from `--ssh-key-file` (stored in metadata
     for use by `cloud_run.py`)

2. **Remove stale host keys** for the IP from `~/.ssh/known_hosts`
   via `ssh-keygen -R` (only on new instance creation).  Both bare
   and bracketed forms are removed for IPv6 addresses (e.g.,
   `ssh-keygen -R 2a01:4f8:...::1` and
   `ssh-keygen -R [2a01:4f8:...::1]`).  This prevents
   `REMOTE HOST IDENTIFICATION HAS CHANGED` errors when Hetzner
   reuses an IP from a previously deleted instance.

3. **Wait for SSH** (poll until `ssh root@IP true` succeeds, with
   `StrictHostKeyChecking=accept-new` to auto-accept the new host
   key).

4. **Install make** on the instance (only on new instance creation):
   ```
   apt-get update -qq && apt-get install -y -qq make
   ```
   Ubuntu 24.04 minimal doesn't include make by default.

5. **Create study directory** on the instance:
   ```
   ssh root@IP "mkdir -p /root/wcss/<study>/sims"
   ```

6. **Upload files** via rsync:
   ```
   rsync -avz --include='*/' \
     --include='Makefile' \
     --include='sim_*/sim.maple' \
     --exclude='*' \
     <study>/sims/ root@[IP]:/root/wcss/<study>/sims/
   ```
   Plus (only on first setup for this instance), upload the binary
   specified by `--delphy` (default: `~/now/delphy/build/release/delphy`):
   ```
   scp <delphy-binary> root@[IP]:/root/wcss/delphy
   ssh root@IP "chmod +x /root/wcss/delphy"
   ```

7. **Fix up Makefile DELPHY path** on the instance:
   ```
   ssh root@IP "sed -i 's|DELPHY = .*|DELPHY = /root/wcss/delphy|' /root/wcss/<study>/sims/Makefile"
   ```

8. **Print summary**: instance name, IP, sim count, and a suggested
   `cloud_run.py` command with the correct `--sims` range based on
   the actual number of sim directories found locally.

### Output

Writes instance metadata (name, IP, ssh_key_file, studies) to
`.wcss-cloud/<name>.json` in the current directory for use by
`cloud_run.py`.  On subsequent calls with the same `--name`, the
new study is appended to the metadata.  The `ssh_key_file` is read
from the metadata by `cloud_run.py` so it doesn't need to be
specified again.

## Script 2: `cloud_run.py`

### Usage

```bash
# Kick off runs and monitor:
./cloud_run.py \
  --name wcss-06-a \
  --study 06_missing_data \
  --sims 0-249 \
  --jobs 10 \
  --sync-every 120 \
  --no-dphy-files

# Resume monitoring only (no new make invocation):
./cloud_run.py \
  --name wcss-06-a \
  --study 06_missing_data \
  --sims 0-249 \
  --resume \
  --sync-every 120 \
  --no-dphy-files
```

### What it does

1. **Read instance metadata** from `.wcss-cloud/<name>.json` to
   get IP and SSH key file path.

2. **Start the runs** on the instance via SSH (skipped with
   `--resume`):

   First, check for a lock file `/root/wcss/<study>/make.lock` on the
   instance.  If it exists, abort with an error message explaining
   that another `make` is already running and suggesting `--resume`.
   If not, write a launcher script and start it:

   The make command is written to a launcher script on the instance
   (to avoid SSH channel FD inheritance issues with inline commands),
   with a unique filename including timestamp and PID:

   ```bash
   # /root/wcss/<study>/run_2026-03-18-14-30-00_12345.sh
   #!/bin/bash
   touch /root/wcss/<study>/make.lock
   cd /root/wcss/<study>/sims
   make -j10 sim_000/.done sim_001/.done ... sim_249/.done
   rm -f /root/wcss/<study>/make.lock
   ```

   Then executed via:
   ```bash
   ssh root@IP "/root/wcss/<study>/run_2026-03-18-14-30-00_12345.sh </dev/null >/root/wcss/<study>/make.log 2>&1 & disown"
   ```

   The lock file is removed automatically when `make` finishes
   (whether it succeeds or fails).  The `make` targets are generated
   from the `--sims` range.

   **Why a launcher script + `disown`**: passing `make -jN` inline
   to SSH causes SSH to hang because `make`'s child processes inherit
   the SSH channel's file descriptors.  Writing the command to a
   script avoids this.  `disown` removes the job from bash's job
   table so the shell exits immediately.  Redirecting all three
   stdio streams (`</dev/null >log 2>&1`) fully detaches the process.
   The `make` process survives SSH disconnection and `cloud_run.py`
   termination.

3. **Monitor loop** (runs locally, every `--sync-every` seconds):

   a. **Report progress**: SSH in and gather status for the
      requested sim range:
      - Count `.done` markers (completed sims).
      - For running sims (no `.done`, but `delphy.log` exists),
        read the last log line to extract the current step count
        (first column of the tab-separated log).
      - The total steps per run is extracted from the Makefile's
        `--v0-steps` argument (parsed once at startup from the
        local copy of the Makefile).
      - Display an in-place terminal update (using ANSI escape
        codes to overwrite only the previous progress block):
        ```
        [14:32:05] 47/250 done (18.8%)  ETA: 2h 15m
        Running:  sim_048 ████████░░ 16.2M/20M
                  sim_049 ██████░░░░ 12.1M/20M
                  sim_050 ██░░░░░░░░  4.3M/20M
                  sim_051 █████████░ 18.7M/20M
                  ... (up to ~jobs active sims)
        ```
        The ETA is computed from the rate of progress since
        monitoring started (not since the runs started).  On the
        first progress check, the current overall progress is
        recorded as a baseline; subsequent ETA calculations use
        the delta from that baseline, so ETA works correctly with
        `--resume`.

   b. **Check if all done**: if done count == total targets, perform
      a final rsync, print summary, and exit.

   c. **Rsync results back** with a rolling display of the last 10
      rsync output lines (stdout and stderr interleaved), shown
      indented below a "Syncing results..." status message.  The
      status message and rsync output are cleared after rsync
      completes:
      ```bash
      rsync -avz --include='*/' \
        --include='sim_*/delphy.log' \
        --include='sim_*/delphy.trees' \
        --include='sim_*/.done' \
        --exclude='*' \
        root@[IP]:/root/wcss/<study>/sims/ <local-study>/sims/
      ```
      (With `--no-dphy-files`, `delphy.dphy` is excluded; otherwise
      `--include='sim_*/delphy.dphy'` is added.)

   d. **Show next update time**: after rsync completes, display
      "Next update at HH:MM:SS" (computed from current time +
      `--sync-every`).

4. On Ctrl-C, print a message noting the runs continue on the
   instance (since `make` runs detached via `disown`) and can be
   re-attached with `--resume`.

### The `--resume` flag

`--resume` skips step 2 (launching `make`) and goes straight to the
monitor loop.  Use cases:
- Reconnecting after a local disconnect or Ctrl-C
- Checking progress without launching duplicate `make` instances
- Monitoring runs that were started by a previous `cloud_run.py`
  invocation

### Splitting across multiple instances

To split 1000 replicates across 4 instances:

```bash
./cloud_setup.py --study 06_missing_data --name wcss-06-a --ssh-key hetzner --ssh-key-file ~/.ssh/id_hetzner
./cloud_setup.py --study 06_missing_data --name wcss-06-b --ssh-key hetzner --ssh-key-file ~/.ssh/id_hetzner
./cloud_setup.py --study 06_missing_data --name wcss-06-c --ssh-key hetzner --ssh-key-file ~/.ssh/id_hetzner
./cloud_setup.py --study 06_missing_data --name wcss-06-d --ssh-key hetzner --ssh-key-file ~/.ssh/id_hetzner

# In 4 terminals:
./cloud_run.py --name wcss-06-a --study 06_missing_data --sims 0-249
./cloud_run.py --name wcss-06-b --study 06_missing_data --sims 250-499
./cloud_run.py --name wcss-06-c --study 06_missing_data --sims 500-749
./cloud_run.py --name wcss-06-d --study 06_missing_data --sims 750-999
```

## Script 3: `cloud_teardown.py`

### Usage

```bash
./cloud_teardown.py --name wcss-06-a
```

### What it does

1. Delete the Hetzner server via `hcloud server delete <name>`.
2. Remove the local metadata file `.wcss-cloud/<name>.json`.
   If `.wcss-cloud/` is now empty, remove the directory too.
3. Print confirmation.

## File layout on the instance

```
/root/wcss/
  delphy              # binary (shared across studies)
  <study>/
    make.log          # stdout/stderr from make
    make.lock         # lock file (present while make is running)
    run_*.sh          # launcher script(s) (one per cloud_run.py invocation)
    sims/
      Makefile
      sim_000/
        sim.maple        # input (uploaded)
        delphy.log       # output (rsynced back)
        delphy.trees     # output (rsynced back)
        delphy.dphy      # output (optionally rsynced back)
        .done            # output (rsynced back)
      sim_001/
        ...
```

## Implementation notes

- All three scripts are Python 3, using `subprocess` to call
  `hcloud`, `ssh`, `rsync`, `scp`.  No special Python packages
  required.
- The instance needs `make` installed (`apt-get install make`),
  which `cloud_setup.py` does automatically on new instances.
  `delphy` itself is a fully self-contained Linux x86_64 binary.
- Rsync uses `--compress` for the transfer; the `.maple` inputs are
  small (~64K each) so even 1000 of them is only ~64MB.
- A per-study lock file (`/root/wcss/<study>/make.lock`) prevents
  concurrent `make` instances for the same study.  It is created
  before `make` starts and removed when `make` exits (via the
  launcher script).  `--resume` skips the lock check since it
  doesn't launch `make`.
- SSH commands use `-o StrictHostKeyChecking=accept-new` to
  auto-accept new host keys without prompting (safe for freshly
  created instances).  All SSH/rsync/scp calls pass `-i <key_file>`
  when a `--ssh-key-file` was provided (stored in the metadata and
  read automatically by `cloud_run.py`).
- The total steps per run is parsed from the local Makefile's
  `--v0-steps` argument at startup (e.g., `--v0-steps 20000000`
  → 20M).  The current step is read from the first column of the
  last line of each sim's `delphy.log` on the instance.
- Instances are created with `--without-ipv4` (IPv6 only) to avoid
  the IPv4 surcharge.  The server's IPv6 address is derived from its
  `/64` network (e.g., `2a01:4f8:1c1e:da48::/64` →
  `2a01:4f8:1c1e:da48::1`).  SSH uses the bare address; rsync and
  scp use brackets (`[addr]`) in the `host:path` syntax to
  disambiguate from the path separator.
- Metadata is stored in `.wcss-cloud/` in the current directory
  (typically `wcss/`), not in `~/.wcss-cloud/`.

## Scripts location

All three scripts go in the `wcss/` directory alongside the studies:

```
wcss/
  cloud_setup.py
  cloud_run.py
  cloud_teardown.py
  .wcss-cloud/
    wcss-06-a.json
    ...
  02_simple/
  03_free_mu/
  ...
```
