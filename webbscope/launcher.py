"""Starts the local server, then shows the page in a native window (pywebview) or, as a fallback, the default browser."""
import argparse, socket, sys, threading, time, webbrowser
import uvicorn
from . import __version__, server


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class Api:
    """Exposed to the page as window.pywebview.api (native folder picker, opening report pages in the default browser)."""
    def __init__(self, base=""):
        self.base = base

    def open_url(self, path):
        """QC report / Profile / validation report: pages the native window cannot open in a new tab. Only paths of this server are opened."""
        if not isinstance(path, str) or not path.startswith("/") or path.startswith("//") or "\\" in path or not self.base:
            return False
        return webbrowser.open(self.base + path)

    def pick_folder(self):
        import webview
        kind = getattr(getattr(webview, "FileDialog", None), "FOLDER", None) or webview.FOLDER_DIALOG
        r = webview.windows[0].create_file_dialog(kind)
        return r[0] if r else None

    def pick_file(self, types=None):
        """types: optional list of extensions without the dot, e.g. ["czi"]."""
        import webview
        kind = getattr(getattr(webview, "FileDialog", None), "OPEN", None) or webview.OPEN_DIALOG
        ft = (f"Files ({';'.join('*.' + t for t in types)})", "All files (*.*)") if types else ()
        r = webview.windows[0].create_file_dialog(kind, file_types=ft)
        return r[0] if r else None


def main(argv=None):
    ap = argparse.ArgumentParser(prog="webbscope", description="Spatial transcriptomics viewer (proof of concept)")
    ap.add_argument("project", nargs="?", help="project folder to open at start")
    ap.add_argument("--port", type=int, default=0, help="default: a free port")
    ap.add_argument("--browser", action="store_true", help="use the default web browser instead of a native window")
    ap.add_argument("--headless", action="store_true", help="server only, no window and no browser (default on Linux)")
    ap.add_argument("--window", action="store_true", help="Linux only: try the native window (needs a display and pywebview) instead of server-only mode")
    ap.add_argument("--host", default="127.0.0.1", help="address to listen on (default: 127.0.0.1, reachable only from this machine / through an SSH tunnel)")
    ap.add_argument("--force", action="store_true", help="redo the import even if the cache is current")
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args(argv)

    if sys.platform.startswith("linux") and not a.window and not a.browser:
        a.headless = True                                # HPC / servers: no display, so serve and let the user connect with a browser
    port = a.port or _free_port()
    srv = uvicorn.Server(uvicorn.Config(server.app, host=a.host, port=port, log_level="warning"))
    threading.Thread(target=srv.run, daemon=True).start()
    t0 = time.time()
    while not srv.started:
        if time.time() - t0 > 30:
            sys.exit("server did not start")
        time.sleep(0.02)
    url = f"http://127.0.0.1:{port}"
    print(f"WebbScope {__version__} → {url}", flush=True)
    if a.headless and sys.platform.startswith("linux"):
        host = socket.gethostname()
        print(f"\nOn your own computer run (jumps through the login node to this machine, {host}):\n"
              f"    ssh -N -L {port}:localhost:{port} -J <you>@<login-node> <you>@{host}\n"
              f"(on the login node itself: ssh -N -L {port}:localhost:{port} <you>@<login-node>)\n"
              f"then open http://localhost:{port} in a browser. Stop with Ctrl-C.\n", flush=True)
    if a.project:
        server.open_project(a.project, a.force)

    if a.headless or a.browser:
        if a.browser:
            webbrowser.open(url)
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            pass
        return
    try:
        import webview
    except ImportError:
        print("pywebview not available; opening the browser instead")
        webbrowser.open(url)
        while True:
            time.sleep(1)
    webview.settings["ALLOW_DOWNLOADS"] = True      # Save PNG / Export CSV in the native window
    webview.create_window("WebbScope", url, js_api=Api(url), width=1400, height=900)
    webview.start()
