#!/usr/bin/env python3
"""Commission a Hetzner Cloud instance and upload a WCSS study."""

import argparse
import glob
import json
import os
import subprocess
import sys
import time

METADATA_DIR = ".wcss-cloud"
BASE_SSH_OPTS = ["-o", "StrictHostKeyChecking=accept-new", "-o", "BatchMode=yes"]


def ssh_opts(key_file=None):
    opts = list(BASE_SSH_OPTS)
    if key_file:
        opts.extend(["-i", key_file])
    return opts


def run(cmd, **kwargs):
    print(f"  $ {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, **kwargs)


def get_server_ipv6(name):
    """Get the server's IPv6 address (first address in its /64 network)."""
    result = subprocess.run(
        ["hcloud", "server", "describe", name, "-o", "json"],
        capture_output=True, text=True, check=True,
    )
    info = json.loads(result.stdout)
    network = info["public_net"]["ipv6"]["ip"]  # e.g. "2a01:4f8:1c1e:da48::/64"
    prefix = network.split("/")[0]  # "2a01:4f8:1c1e:da48::"
    if prefix.endswith("::"):
        return prefix + "1"  # "2a01:4f8:1c1e:da48::1"
    return prefix


def rsync_host(ip):
    """Format IP for rsync/scp host:path syntax (IPv6 needs brackets)."""
    if ":" in ip:
        return f"[{ip}]"
    return ip


def wait_for_ssh(ip, key_file=None, timeout=120):
    print(f"Waiting for SSH on {ip}...")
    opts = ssh_opts(key_file)
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            result = subprocess.run(
                ["ssh", *opts, f"root@{ip}", "true"],
                capture_output=True, timeout=10,
            )
            if result.returncode == 0:
                print("  SSH ready.")
                return
        except subprocess.TimeoutExpired:
            pass
        time.sleep(2)
    print(f"ERROR: SSH not available on {ip} after {timeout}s", file=sys.stderr)
    sys.exit(1)


def load_metadata(name):
    path = os.path.join(METADATA_DIR, f"{name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def save_metadata(name, meta):
    os.makedirs(METADATA_DIR, exist_ok=True)
    path = os.path.join(METADATA_DIR, f"{name}.json")
    with open(path, "w") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, help="Study directory name (e.g. 06_missing_data)")
    parser.add_argument("--name", required=True, help="Hetzner instance name")
    parser.add_argument("--type", default="cpx52", help="Server type (default: cpx52)")
    parser.add_argument("--location", default="nbg1", help="Server location (default: nbg1)")
    parser.add_argument("--ssh-key", help="SSH key name registered with Hetzner")
    parser.add_argument("--ssh-key-file", help="Path to local SSH private key file (e.g. ~/.ssh/id_hetzner)")
    parser.add_argument("--delphy", default=os.path.expanduser("~/now/delphy/build/release/delphy"),
                        help="Path to delphy binary (default: ~/now/delphy/build/release/delphy)")
    args = parser.parse_args()

    # Validate local study directory
    study_sims = os.path.join(args.study, "sims")
    if not os.path.isdir(study_sims):
        print(f"ERROR: {study_sims}/ not found. Run from the wcss/ directory.", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(os.path.join(study_sims, "Makefile")):
        print(f"ERROR: {study_sims}/Makefile not found.", file=sys.stderr)
        sys.exit(1)

    # Check if instance already exists
    meta = load_metadata(args.name)
    is_new_instance = meta is None

    # Resolve SSH key file path
    key_file = args.ssh_key_file
    if key_file:
        key_file = os.path.expanduser(key_file)

    if is_new_instance:
        # Validate delphy binary
        if not os.path.isfile(args.delphy):
            print(f"ERROR: Delphy binary not found at {args.delphy}", file=sys.stderr)
            sys.exit(1)

        if not args.ssh_key:
            print("ERROR: --ssh-key is required when creating a new instance.", file=sys.stderr)
            sys.exit(1)

        # Step 1: Create instance
        print(f"Creating instance {args.name}...")
        run([
            "hcloud", "server", "create",
            "--name", args.name,
            "--type", args.type,
            "--location", args.location,
            "--image", "ubuntu-24.04",
            "--ssh-key", args.ssh_key,
            "--without-ipv4",
        ])
        ip = get_server_ipv6(args.name)
        meta = {"name": args.name, "ip": ip, "studies": []}
        if key_file:
            meta["ssh_key_file"] = key_file
    else:
        ip = meta["ip"]
        if not key_file:
            key_file = meta.get("ssh_key_file")
        if args.study in meta.get("studies", []):
            print(f"Study {args.study} already set up on {args.name} ({ip}).")
            print("Re-uploading files...")
        else:
            print(f"Instance {args.name} ({ip}) already exists, adding study {args.study}...")

    opts = ssh_opts(key_file)

    # Remove any stale host key for this IP (common when IPs are reused
    # by new instances), so StrictHostKeyChecking=accept-new works.
    # Remove both bare and bracketed forms for IPv6 addresses.
    if is_new_instance:
        subprocess.run(["ssh-keygen", "-R", ip],
                       capture_output=True, check=False)
        if ":" in ip:
            subprocess.run(["ssh-keygen", "-R", f"[{ip}]"],
                           capture_output=True, check=False)

    # Step 2: Wait for SSH
    wait_for_ssh(ip, key_file)

    # Install make if not present (Ubuntu 24.04 minimal doesn't include it)
    if is_new_instance:
        print("Installing make on instance...")
        run(["ssh", *opts, f"root@{ip}",
             "apt-get update -qq && apt-get install -y -qq make"])

    # Step 3: Create study directory
    print(f"Creating /root/wcss/{args.study}/sims on instance...")
    run(["ssh", *opts, f"root@{ip}", f"mkdir -p /root/wcss/{args.study}/sims"])

    # Step 4: Upload files
    print(f"Uploading {study_sims}/ to instance...")
    run([
        "rsync", "-avz",
        "--include=*/",
        "--include=Makefile",
        "--include=sim_*/sim.maple",
        "--exclude=*",
        "-e", f"ssh {' '.join(opts)}",
        f"{study_sims}/",
        f"root@{rsync_host(ip)}:/root/wcss/{args.study}/sims/",
    ])

    if is_new_instance:
        print("Uploading delphy binary...")
        run(["scp", *opts, args.delphy, f"root@{rsync_host(ip)}:/root/wcss/delphy"])
        run(["ssh", *opts, f"root@{ip}", "chmod +x /root/wcss/delphy"])

    # Step 5: Fix up Makefile DELPHY path
    print("Fixing Makefile DELPHY path...")
    run(["ssh", *opts, f"root@{ip}",
         f"sed -i 's|DELPHY = .*|DELPHY = /root/wcss/delphy|' /root/wcss/{args.study}/sims/Makefile"])

    # Update metadata
    if args.study not in meta.get("studies", []):
        meta.setdefault("studies", []).append(args.study)
    save_metadata(args.name, meta)

    # Count sim directories
    sim_dirs = sorted(glob.glob(os.path.join(study_sims, "sim_[0-9]*")))
    last_sim = len(sim_dirs) - 1 if sim_dirs else 0

    # Step 6: Print summary
    print()
    print(f"Instance:  {args.name}")
    print(f"IP:        {ip}")
    print(f"Study:     {args.study}")
    print(f"Sims:      {len(sim_dirs)}")
    print(f"Metadata:  {METADATA_DIR}/{args.name}.json")
    print()
    print("To start runs:")
    print(f"  ./cloud_run.py --name {args.name} --study {args.study} --sims 0-{last_sim}")


if __name__ == "__main__":
    main()
