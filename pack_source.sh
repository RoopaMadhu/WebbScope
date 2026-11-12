#!/usr/bin/env bash
# Make a small tarball of the source to copy to the HPC (no venv / builds / caches): scp WebbScope-v<version>-src.tar.gz <hpc>:
# Everything unpacks into one top-level WebbScope/ folder:  tar xzf WebbScope-v<version>-src.tar.gz && cd WebbScope
# The version comes from webbscope/__init__.py.
set -euo pipefail
cd "$(dirname "$0")"
NAME=WebbScope
VERSION="$(sed -n 's/^__version__ = "\(.*\)"$/\1/p' webbscope/__init__.py)"
OUT="$NAME-v$VERSION-src.tar.gz"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
mkdir "$STAGE/$NAME"
rsync -a --exclude='__pycache__' --exclude='.DS_Store' \
  webbscope run_app.py build_linux.sh requirements-linux.txt README.md docs "$STAGE/$NAME/"
COPYFILE_DISABLE=1 tar -C "$STAGE" -czf "$OUT" "$NAME"
echo "made $OUT ($(du -h "$OUT" | cut -f1))"
