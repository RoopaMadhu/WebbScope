"""Starts the local server, then shows the page in a native window (pywebview) or, as a fallback, the default browser."""
import argparse, getpass, secrets, signal, socket, sys, threading, time, webbrowser
import uvicorn
from . import __version__, access, server


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
        return webbrowser.open(access.link(self.base, path))    # the default browser has no cookie of this server yet

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


def _set_token():
    """--set-token: the user's own token, saved privately. Typed, never on the command line, where other users of a node could see it."""
    if sys.stdin.isatty():
        t = getpass.getpass("New access token: ")
        if t != getpass.getpass("Again: "):
            sys.exit("the two entries differ; nothing was changed")
    else:
        t = sys.stdin.readline()
    try:
        f = access.save_token(t)
    except (ValueError, OSError) as e:
        sys.exit(f"token not saved: {e}")
    print(f"saved in {f}. Links opened with the previous token stop working when the server restarts.")


def _link_lines(port):
    """What to show for opening the server in a browser. A saved or chosen token stays valid across jobs, so it is not written where a log file could keep it."""
    base = f"http://localhost:{port}"
    if access.SOURCE == "random":
        return f"    {access.link(base)}\n(this token is for this start only: after a restart, open the new link)"
    where = "~/.webbscope/access_token" if access.SOURCE == "saved" else "WEBBSCOPE_TOKEN"
    keeps = "it stays the same if the job restarts, so a page left open keeps working"
    if sys.stdout.isatty():
        return f"    {access.link(base)}\n(treat it like a password; {keeps})"
    return f"    {base}/?token=<your token>\n(your token is in {where}; not printed because this output may be a log file; {keeps})"


def main(argv=None):
    ap = argparse.ArgumentParser(prog="webbscope", description="Spatial transcriptomics viewer (proof of concept)")
    ap.add_argument("project", nargs="?", help="project folder to open at start")
    ap.add_argument("--port", type=int, default=0, help="default: a free port")
    ap.add_argument("--browser", action="store_true", help="use the default web browser instead of a native window")
    ap.add_argument("--headless", action="store_true", help="server only, no window and no browser (default on Linux)")
    ap.add_argument("--window", action="store_true", help="Linux only: try the native window (needs a display and pywebview) instead of server-only mode")
    ap.add_argument("--host", default="127.0.0.1", help="address to listen on (default: 127.0.0.1, reachable only from this machine / through an SSH tunnel)")
    ap.add_argument("--force", action="store_true", help="redo the import even if the cache is current")
    tk = ap.add_mutually_exclusive_group()
    tk.add_argument("--set-token", action="store_true", help="choose the access token kept for server (HPC) mode, then exit (it is typed, not put on the command line)")
    tk.add_argument("--new-token", action="store_true", help="replace the kept access token with a new random one, then exit (old links stop working)")
    tk.add_argument("--random-token", action="store_true", help="use a new token for this start only, not kept across restarts")
    ap.add_argument("--version", action="version", version=f"v{__version__}")
    a = ap.parse_args(argv)
    if a.set_token:
        return _set_token()
    if a.new_token:
        try:
            f = access.save_token(secrets.token_urlsafe(32))
        except OSError as e:
            sys.exit(f"token not saved: {e}")
        return print(f"new token saved in {f}. Old links stop working when the server restarts; the next start shows how to open it.")

    if sys.platform.startswith("linux") and not a.window and not a.browser:
        a.headless = True                                # HPC / servers: no display, so serve and let the user connect with a browser
    port = a.port or _free_port()
    try:
        access.configure(fixed=a.headless, one_off=a.random_token)   # before the server accepts anything; server (HPC) mode keeps its token
    except ValueError as e:
        sys.exit(str(e))
    srv = uvicorn.Server(uvicorn.Config(server.app, host=a.host, port=port, log_level="warning"))
    threading.Thread(target=srv.run, daemon=True).start()
    t0 = time.time()
    while not srv.started:
        if time.time() - t0 > 30:
            sys.exit("server did not start")
        time.sleep(0.02)
    url = f"http://127.0.0.1:{port}"
    entry = access.link(url)                             # the address with the token: the app's own window and browser start from it
    access.write_connection_file(port, __version__)
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))   # a job scheduler stops the server this way: exit normally, so the connection file is removed
    print(f"WebbScope v{__version__} → {url}", flush=True)
    for w in access.WARNINGS:
        print(f"WARNING: {w}", flush=True)
    if a.headless and sys.platform.startswith("linux"):
        host = socket.gethostname()
        print(f"\nOn your own computer run (jumps through the login node to this machine, {host}):\n"
              f"    ssh -N -L {port}:localhost:{port} -J <you>@<login-node> <you>@{host}\n"
              f"(on the login node itself: ssh -N -L {port}:localhost:{port} <you>@<login-node>)\n"
              f"then open this link in a browser (if you tunnel to another local port, change the port in it):\n"
              f"{_link_lines(port)}\n"
              f"Stop with Ctrl-C.\n", flush=True)
    elif a.headless:
        print(f"Open this link in a browser:\n{_link_lines(port)}\n", flush=True)
    if a.project:
        server.open_project(a.project, a.force)

    if a.headless or a.browser:
        if a.browser:
            webbrowser.open(entry)
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
        webbrowser.open(entry)
        while True:
            time.sleep(1)
    webview.settings["ALLOW_DOWNLOADS"] = True      # Save PNG / Export CSV in the native window
    webview.create_window("WebbScope", entry, js_api=Api(url), width=1400, height=900)
    webview.start()
