#!/usr/bin/env python3
"""Run all Delphy jobs via the generated Makefile.

Usage:
    ./02_run.py             # uses half the available CPUs
    ./02_run.py --jobs 20   # specify parallelism explicitly
"""

import argparse
import os
import subprocess
import sys


def main():
    default_jobs = max(1, os.cpu_count() // 2)

    parser = argparse.ArgumentParser(description="Run all Delphy jobs")
    parser.add_argument("--jobs", "-j", type=int, default=default_jobs,
                        help=f"Number of parallel jobs (default: {default_jobs})")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    sims_dir = os.path.join(script_dir, "sims")

    print(f"Running make -j{args.jobs} in {sims_dir} ...")
    result = subprocess.run(["make", f"-j{args.jobs}"], cwd=sims_dir)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
