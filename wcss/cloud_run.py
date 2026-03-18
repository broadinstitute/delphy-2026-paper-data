#!/usr/bin/env python3
"""Kick off WCSS runs on a Hetzner instance and monitor progress."""

import argparse
import json
import os
import re
import signal
import subprocess
import sys
import time

METADATA_DIR = ".wcss-cloud"
BASE_SSH_OPTS = ["-o", "StrictHostKeyChecking=accept-new", "-o", "BatchMode=yes"]

# Set by main() after loading metadata
_ssh_opts = list(BASE_SSH_OPTS)


def ssh_opts():
    return _ssh_opts


def rsync_host(ip):
    """Format IP for rsync/scp host:path syntax (IPv6 needs brackets)."""
    if ":" in ip:
        return f"[{ip}]"
    return ip


def ssh_run(ip, cmd, **kwargs):
    return subprocess.run(
        ["ssh", *ssh_opts(), f"root@{ip}", cmd],
        capture_output=True, text=True, **kwargs,
    )


def parse_sims_range(sims_str):
    """Parse '0-249' or '0-249,500-749' into a list of sim indices."""
    indices = []
    for part in sims_str.split(","):
        if "-" in part:
            start, end = part.split("-", 1)
            indices.extend(range(int(start), int(end) + 1))
        else:
            indices.append(int(part))
    return sorted(indices)


def sim_dir_name(idx):
    return f"sim_{idx:03d}"


def parse_total_steps(study):
    """Extract --v0-steps value from the local Makefile."""
    makefile = os.path.join(study, "sims", "Makefile")
    with open(makefile) as f:
        content = f.read()
    m = re.search(r"--v0-steps\s+(\d+)", content)
    if not m:
        print(f"ERROR: Could not find --v0-steps in {makefile}", file=sys.stderr)
        sys.exit(1)
    return int(m.group(1))


def format_steps(n):
    """Format step count as e.g. '16.2M' or '500K'."""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f}M"
    if n >= 1_000:
        return f"{n / 1_000:.0f}K"
    return str(n)


def make_progress_bar(current, total, width=10):
    filled = int(width * current / total) if total > 0 else 0
    filled = min(filled, width)
    return "\u2588" * filled + "\u2591" * (width - filled)


def gather_progress(ip, study, sim_indices):
    """SSH in and get status for all sims in one command."""
    # Build a shell script that checks each sim
    checks = []
    for idx in sim_indices:
        d = sim_dir_name(idx)
        checks.append(
            f'if [ -f /root/wcss/{study}/sims/{d}/.done ]; then echo "{idx} done"; '
            f'elif [ -f /root/wcss/{study}/sims/{d}/delphy.log ]; then '
            f'echo "{idx} $(tail -1 /root/wcss/{study}/sims/{d}/delphy.log | cut -f1)"; '
            f'else echo "{idx} pending"; fi'
        )
    script = "; ".join(checks)
    result = ssh_run(ip, script, timeout=30)
    if result.returncode != 0:
        return None

    status = {}  # idx -> "done" | step_count (int) | "pending"
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split(None, 1)
        idx = int(parts[0])
        val = parts[1] if len(parts) > 1 else "pending"
        if val == "done":
            status[idx] = "done"
        elif val == "pending":
            status[idx] = "pending"
        else:
            try:
                status[idx] = int(val)
            except ValueError:
                status[idx] = "pending"
    return status


_prev_progress_lines = 0
_start_time = None
_start_overall = None


def format_duration(seconds):
    """Format seconds as e.g. '2h 15m' or '45m' or '3m'."""
    seconds = int(seconds)
    if seconds < 60:
        return f"{seconds}s"
    minutes = seconds // 60
    if minutes < 60:
        return f"{minutes}m"
    hours = minutes // 60
    mins = minutes % 60
    return f"{hours}h {mins:02d}m"


def display_progress(status, total_steps, sim_indices):
    """Print progress display, overwriting only the previous progress block."""
    global _prev_progress_lines, _start_time, _start_overall

    if _start_time is None:
        _start_time = time.time()

    done_count = sum(1 for v in status.values() if v == "done")
    total_sims = len(sim_indices)
    pct = 100.0 * done_count / total_sims if total_sims > 0 else 0

    # Compute overall progress including partial progress of running sims
    running_progress = sum(
        v / total_steps for v in status.values() if isinstance(v, int)
    )
    overall = (done_count + running_progress) / total_sims if total_sims > 0 else 0

    # Record initial progress on first check (important for --resume)
    if _start_overall is None:
        _start_overall = overall

    # Estimate time remaining based on progress delta since monitoring started
    elapsed = time.time() - _start_time
    eta_str = ""
    progress_delta = overall - _start_overall
    if progress_delta > 0:
        remaining = elapsed * (1.0 - overall) / progress_delta
        eta_str = f"  ETA: {format_duration(remaining)}"

    now = time.strftime("%H:%M:%S")
    lines = [f"[{now}] {done_count}/{total_sims} done ({pct:.1f}%){eta_str}"]

    # Collect running sims
    running = []
    for idx in sim_indices:
        s = status.get(idx)
        if isinstance(s, int):
            running.append((idx, s))

    if running:
        row_items = []
        for idx, step in running:
            bar = make_progress_bar(step, total_steps)
            label = f"{sim_dir_name(idx)} {bar} {format_steps(step)}/{format_steps(total_steps)}"
            row_items.append(label)

        for i, item in enumerate(row_items):
            prefix = "Running:  " if i == 0 else "          "
            lines.append(prefix + item)

    # Move cursor up to overwrite previous progress block, clear to end of screen
    if _prev_progress_lines > 0:
        sys.stdout.write(f"\033[{_prev_progress_lines}A\033[J")

    output = "\n".join(lines)
    sys.stdout.write(f"{output}\n")
    sys.stdout.flush()
    _prev_progress_lines = len(lines)

    return done_count


def show_status(msg):
    """Append a status message to the progress block."""
    global _prev_progress_lines
    sys.stdout.write(f"{msg}\n")
    sys.stdout.flush()
    _prev_progress_lines += 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True, help="Instance name")
    parser.add_argument("--study", required=True, help="Study directory name")
    parser.add_argument("--sims", required=True, help="Sim range (e.g. 0-249)")
    parser.add_argument("--jobs", type=int, default=10, help="Number of parallel jobs (default: 10)")
    parser.add_argument("--sync-every", type=int, default=120, help="Seconds between syncs (default: 120)")
    parser.add_argument("--resume", action="store_true", help="Skip launching make, just monitor")
    parser.add_argument("--no-dphy-files", action="store_true", help="Don't rsync delphy.dphy files back")
    args = parser.parse_args()

    # Load metadata
    meta_path = os.path.join(METADATA_DIR, f"{args.name}.json")
    if not os.path.exists(meta_path):
        print(f"ERROR: {meta_path} not found. Run cloud_setup.py first.", file=sys.stderr)
        sys.exit(1)
    with open(meta_path) as f:
        meta = json.load(f)
    ip = meta["ip"]

    # Set up SSH options with key file from metadata
    global _ssh_opts
    key_file = meta.get("ssh_key_file")
    if key_file:
        _ssh_opts = list(BASE_SSH_OPTS) + ["-i", key_file]

    sim_indices = parse_sims_range(args.sims)
    total_steps = parse_total_steps(args.study)
    total_sims = len(sim_indices)
    study = args.study

    print(f"Instance:    {args.name} ({ip})")
    print(f"Study:       {study}")
    print(f"Sims:        {sim_indices[0]}-{sim_indices[-1]} ({total_sims} total)")
    print(f"Total steps: {format_steps(total_steps)}")
    print()

    # Step 2: Start runs (unless --resume)
    if not args.resume:
        # Check lock file
        result = ssh_run(ip, f"test -f /root/wcss/{study}/make.lock")
        if result.returncode == 0:
            print(f"ERROR: Lock file exists on instance — make is already running for {study}.",
                  file=sys.stderr)
            print("Use --resume to monitor the existing run.", file=sys.stderr)
            sys.exit(1)

        # Write a launcher script on the instance to avoid SSH hanging
        # (inline commands with make -jN keep the SSH channel open because
        # child processes inherit file descriptors from the SSH session)
        targets = " ".join(f"{sim_dir_name(i)}/.done" for i in sim_indices)
        launcher = (
            f"#!/bin/bash\n"
            f"touch /root/wcss/{study}/make.lock\n"
            f"cd /root/wcss/{study}/sims\n"
            f"make -j{args.jobs} {targets}\n"
            f"rm -f /root/wcss/{study}/make.lock\n"
        )
        timestamp = time.strftime("%Y-%m-%d-%H-%M-%S")
        launcher_path = f"/root/wcss/{study}/run_{timestamp}_{os.getpid()}.sh"
        result = ssh_run(ip, f"cat > {launcher_path} << 'LAUNCHER_EOF'\n{launcher}LAUNCHER_EOF\nchmod +x {launcher_path}")
        if result.returncode != 0:
            print(f"ERROR: Failed to create launcher script: {result.stderr}", file=sys.stderr)
            sys.exit(1)

        print("Starting runs...")
        launch_cmd = f"{launcher_path} </dev/null >/root/wcss/{study}/make.log 2>&1 & disown"
        result = ssh_run(ip, launch_cmd)
        print(f"Launched make -j{args.jobs} with {total_sims} targets.")
        print()

    # Build rsync include list
    rsync_includes = [
        "--include=*/",
        "--include=sim_*/delphy.log",
        "--include=sim_*/delphy.trees",
        "--include=sim_*/.done",
    ]
    if not args.no_dphy_files:
        rsync_includes.append("--include=sim_*/delphy.dphy")
    rsync_includes.append("--exclude=*")

    local_sims = os.path.join(study, "sims")

    # Handle Ctrl-C gracefully
    def on_sigint(sig, frame):
        print("\n\nCtrl-C received. Runs continue on the instance.")
        print(f"Re-attach with: ./cloud_run.py --name {args.name} --study {study} --sims {args.sims} --resume")
        sys.exit(0)
    signal.signal(signal.SIGINT, on_sigint)

    # Step 3: Monitor loop
    global _prev_progress_lines
    print("Monitoring progress...")
    while True:
        # 3a: Report progress
        status = gather_progress(ip, study, sim_indices)
        if status is not None:
            done_count = display_progress(status, total_steps, sim_indices)

            # 3c: Check if all done
            if done_count == total_sims:
                print(f"\nAll {total_sims} sims complete!")
                # Final rsync
                print("Final rsync...")
                subprocess.run([
                    "rsync", "-avz",
                    *rsync_includes,
                    "-e", f"ssh {' '.join(ssh_opts())}",
                    f"root@{rsync_host(ip)}:/root/wcss/{study}/sims/",
                    f"{local_sims}/",
                ], check=False)
                print("Done.")
                break
        else:
            print("  (could not reach instance, will retry)")

        # 3b: Rsync results back
        show_status("Syncing results...")
        rsync_proc = subprocess.Popen(
            ["rsync", "-avz",
             *rsync_includes,
             "-e", f"ssh {' '.join(ssh_opts())}",
             f"root@{rsync_host(ip)}:/root/wcss/{study}/sims/",
             f"{local_sims}/"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        # Show a rolling window of the latest rsync output lines
        rsync_lines = []
        max_rsync_lines = 10
        prev_shown = 0
        for line in rsync_proc.stdout:
            line = line.rstrip()
            if not line:
                continue
            rsync_lines.append(line)
            if len(rsync_lines) > max_rsync_lines:
                rsync_lines.pop(0)
            # Overwrite previous rsync output
            if prev_shown > 0:
                sys.stdout.write(f"\033[{prev_shown}A\033[J")
            for rl in rsync_lines:
                sys.stdout.write(f"  {rl}\n")
            sys.stdout.flush()
            prev_shown = len(rsync_lines)
        rsync_proc.wait()
        # Clear "Syncing results..." and rsync output lines
        total_to_clear = prev_shown + 1  # +1 for "Syncing results..."
        sys.stdout.write(f"\033[{total_to_clear}A\033[J")
        sys.stdout.flush()
        _prev_progress_lines -= 1  # undo the show_status increment

        next_update = time.strftime("%H:%M:%S", time.localtime(time.time() + args.sync_every))
        show_status(f"Next update at {next_update}")
        time.sleep(args.sync_every)


if __name__ == "__main__":
    main()
