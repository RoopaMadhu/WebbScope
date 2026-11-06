#!/usr/bin/env bash
# Make a small tarball of the source to copy to the HPC (no venv / builds / caches): scp webbscope-src.tar.gz <hpc>:
set -euo pipefail
cd "$(dirname "$0")"
tar --exclude='__pycache__' --exclude='.DS_Store' -czf webbscope-src.tar.gz webbscope run_app.py build_linux.sh requirements-linux.txt README.md docs
echo "made webbscope-src.tar.gz ($(du -h webbscope-src.tar.gz | cut -f1))"
