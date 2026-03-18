#!/usr/bin/env python3
"""Delete a Hetzner Cloud instance and clean up local metadata."""

import argparse
import json
import os
import subprocess
import sys

METADATA_DIR = ".wcss-cloud"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True, help="Instance name to delete")
    args = parser.parse_args()

    meta_path = os.path.join(METADATA_DIR, f"{args.name}.json")
    if not os.path.exists(meta_path):
        print(f"ERROR: {meta_path} not found.", file=sys.stderr)
        sys.exit(1)

    with open(meta_path) as f:
        meta = json.load(f)

    print(f"Deleting instance {args.name} ({meta['ip']})...")
    result = subprocess.run(["hcloud", "server", "delete", args.name])
    if result.returncode != 0:
        print("ERROR: hcloud server delete failed.", file=sys.stderr)
        sys.exit(1)

    os.remove(meta_path)
    print(f"Removed {meta_path}")

    # Remove .wcss-cloud/ if empty
    try:
        os.rmdir(METADATA_DIR)
        print(f"Removed empty {METADATA_DIR}/")
    except OSError:
        pass  # directory not empty

    print(f"Instance {args.name} deleted.")


if __name__ == "__main__":
    main()
