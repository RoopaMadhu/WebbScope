#!/usr/bin/env bash
# Build a standalone app with its own Python inside. Must be run on the target OS (no cross-compiling).
# orjson is a hidden import: FastAPI imports it lazily, so PyInstaller does not see it (the first build lacked it).
# shapely is imported inside a function (region import), so it is listed explicitly too.
set -euo pipefail
cd "$(dirname "$0")"
.venv/bin/python -m PyInstaller --noconfirm --clean --name WebbScope --windowed \
  --distpath dist --workpath build --specpath build \
  --add-data "$PWD/webbscope/web:webbscope/web" --add-data "$PWD/webbscope/engine/data:webbscope/engine/data" \
  --collect-submodules uvicorn --collect-submodules webview \
  --hidden-import orjson --hidden-import shapely --hidden-import scipy.ndimage --hidden-import PIL --hidden-import readimc --hidden-import h5py --hidden-import imageio.v2 --collect-all pylibCZIrw --collect-submodules webbscope \
  run_app.py
# Zip the app for download; the version comes from webbscope/__init__.py.
VERSION="$(sed -n 's/^__version__ = "\(.*\)"$/\1/p' webbscope/__init__.py)"
ditto -c -k --keepParent dist/WebbScope.app "dist/WebbScope-v$VERSION-mac.zip"
echo "built dist/WebbScope-v$VERSION-mac.zip"
