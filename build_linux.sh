#!/usr/bin/env bash
# Build WebbScope for Linux (x86_64). Run this ON Linux (PyInstaller cannot cross-compile), e.g. on an HPC compute node:
#     bash build_linux.sh            # uses `python3` (3.10+); override with PYTHON=/path/to/python3.12
# Output: dist/WebbScope-linux-x86_64.tar.gz  (unpack anywhere, run ./WebbScope/WebbScope PROJECT_FOLDER)
# Needs internet or a pip mirror for the first run. Everything lands in .venv-linux/, dist/ and build-linux/.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
"$PYTHON" -c 'import sys; assert sys.version_info >= (3, 10), "need Python >= 3.10 (try: module load python/3.12)"'
# PyInstaller needs objdump (package binutils) on Linux. Try the cluster's module system first, then stop with instructions.
if ! command -v objdump >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then module load binutils >/dev/null 2>&1 || true; fi
if ! command -v objdump >/dev/null 2>&1; then
  echo "ERROR: objdump not found (PyInstaller needs it; it comes with 'binutils')." >&2
  echo "  try:  module avail binutils   then   module load binutils" >&2
  echo "  or:   conda install -c conda-forge binutils   (then check:  ls \$CONDA_PREFIX/bin | grep objdump)" >&2
  echo "  or ask the cluster admins for binutils. Then run this script again." >&2
  exit 1
fi
[ -d .venv-linux ] || "$PYTHON" -m venv .venv-linux
.venv-linux/bin/python -m pip install --upgrade pip
.venv-linux/bin/python -m pip install -r requirements-linux.txt
.venv-linux/bin/pyinstaller --noconfirm --clean --name WebbScope \
  --distpath dist --workpath build-linux --specpath build-linux \
  --add-data "$PWD/webbscope/web:webbscope/web" \
  --collect-submodules uvicorn \
  --hidden-import orjson --hidden-import shapely --hidden-import scipy.ndimage --hidden-import PIL --hidden-import readimc --hidden-import imageio.v2 --collect-all pylibCZIrw --collect-submodules webbscope \
  --exclude-module webview \
  run_app.py
echo "--- smoke test"
dist/WebbScope/WebbScope --version
tar -C dist -czf dist/WebbScope-linux-x86_64.tar.gz WebbScope
echo "built dist/WebbScope-linux-x86_64.tar.gz"
