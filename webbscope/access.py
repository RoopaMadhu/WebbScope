"""Access token: only someone holding the server's token can use it (pages, data, API).

127.0.0.1 is "this machine", and on a shared compute node that includes every other user of the node, so the port alone protects nothing there.
Every request must carry the server's token. Which token (first that applies):
  --random-token           a new one for this start only
  WEBBSCOPE_TOKEN          the user's own
  server mode (HPC)        the token saved in ~/.webbscope/access_token (private to the account, created at the first start; --set-token chooses it,
                           --new-token replaces it). It stays the same across job restarts, so a page left open keeps working and keeps its unsaved edits.
  app window / --browser   a new one at each start (the app hands it to its own window)
It is accepted in one of:
  a cookie                 set when the access link (http://host:port/?token=...) is opened; the token is then removed from the address bar
  Authorization: Bearer    for programs
  X-WebbScope-Token        the same, as a plain header
The cookie is HttpOnly and SameSite=Strict, so scripts and other websites cannot read or send it. Programs on the same account find the server's address in
~/.webbscope/connections/<host>-<port>.json (readable only by that account, removed when the server stops), with the agent key, not the token.
Claude's connector (webbscope/mcp.py) gets the agent key instead (agent_key(): derived from the token), which reaches only /api/agent/...: what it may
see there is the sharing level of the project, which only the token (the user's browser) can change (agent.py).
"""
import atexit, hashlib, hmac, json, os, re, secrets, socket, stat, time
from pathlib import Path
from urllib.parse import urlencode
from fastapi.responses import HTMLResponse, JSONResponse

MIN_LENGTH = 16
TOKEN = None                                 # None = not started through the launcher (the app used as a library, tests): no check
COOKIE = "webbscope_access"
SOURCE = None                                # "saved", "env" or "random"
WARNINGS = []                                # for the launcher to print
_conn_file = None


def token_file():
    return Path.home() / ".webbscope" / "access_token"


def validate(t):
    """The token as given, stripped; ValueError if it is too short or has characters that cannot travel in a header."""
    t = (t or "").strip()
    if len(t) < MIN_LENGTH or not re.fullmatch(r"[\x21-\x7e]+", t):
        raise ValueError(f"an access token must be at least {MIN_LENGTH} characters, plain letters, digits or symbols, without spaces")
    return t


def _write_private(path, token):
    """A new private file holding the token, next to `path`; returns its name (the caller moves it into place)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.{secrets.token_hex(4)}")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "w") as fh:
        fh.write(token + "\n")
    return tmp


def save_token(token):
    """Replace the saved token (--set-token, --new-token). Returns the file."""
    f = token_file()
    os.replace(_write_private(f, validate(token)), f)
    return f


def _read_saved():
    """The saved token, or None if there is none. Refuses a link or a file of another user; a file others could read is made private, with a warning."""
    f = token_file()
    try:
        fd = os.open(f, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except FileNotFoundError:
        return None
    with os.fdopen(fd) as fh:
        st = os.fstat(fh.fileno())
        if st.st_uid != os.getuid():
            raise OSError(f"{f} belongs to another user")
        if stat.S_IMODE(st.st_mode) & 0o077:
            os.fchmod(fh.fileno(), 0o600)
            WARNINGS.append(f"{f} could be read by other users; it is private now. If they may have seen it, run with --new-token.")
        t = fh.read()
    try:
        return validate(t)
    except ValueError:
        raise ValueError(f"{f} does not hold a valid token; choose one with --set-token or replace it with --new-token")


def _saved_or_new():
    """The saved token, created at the first start. Two jobs starting at once end up with the same one (the file appears atomically)."""
    t = _read_saved()
    if t is not None:
        return t, "saved"
    f, new = token_file(), secrets.token_urlsafe(32)
    tmp = _write_private(f, new)
    try:
        os.link(tmp, f)                      # fails if another job has just created it
        return new, "saved"
    except FileExistsError:
        return _read_saved(), "saved"
    finally:
        tmp.unlink(missing_ok=True)


def configure(fixed=False, one_off=False):
    """Set the token before the server starts (the order is at the top of this file). fixed = server mode. Returns it."""
    global TOKEN, COOKIE, SOURCE
    WARNINGS.clear()
    env = os.environ.get("WEBBSCOPE_TOKEN")
    if one_off:
        TOKEN, SOURCE = secrets.token_urlsafe(32), "random"
    elif env is not None:
        try:
            TOKEN, SOURCE = validate(env), "env"
        except ValueError as e:
            raise ValueError(f"WEBBSCOPE_TOKEN: {e}")
    elif fixed:
        try:
            TOKEN, SOURCE = _saved_or_new()
        except OSError as e:
            TOKEN, SOURCE = secrets.token_urlsafe(32), "random"
            WARNINGS.append(f"the saved token could not be used ({e}); this start has its own token, which changes if the job restarts.")
    else:
        TOKEN, SOURCE = secrets.token_urlsafe(32), "random"
    COOKIE = "webbscope_access_" + hashlib.sha256(b"webbscope-cookie:" + TOKEN.encode()).hexdigest()[:12]   # one per token, so servers with different
    return TOKEN                             # tokens on one machine do not overwrite each other, and the same token still works after moving to another port


def agent_key():
    """The key of Claude's connector: derived from the token (it changes with it), and of no use beyond /api/agent/. None without a token."""
    if TOKEN is None:
        return None
    return "wsa_" + hmac.new(TOKEN.encode(), b"webbscope-agent-key", hashlib.sha256).hexdigest()[:40]


def link(base, path="/"):
    """The address that lets a browser in: `path` of the server at `base`, with the token."""
    p, _, frag = path.partition("#")
    return f"{base}{p}{'&' if '?' in p else '?'}{urlencode({'token': TOKEN})}" + (f"#{frag}" if frag else "")


def write_connection_file(port, version=""):
    """Where the server is and its token, for programs of this account (the MCP connector). In the home folder, not the cache, which may be a shared
    disk. Readable only by this account. Returns the path, or None if it cannot be written (the server works without it)."""
    global _conn_file
    host = re.sub(r"[^A-Za-z0-9._-]", "_", socket.gethostname())
    d = Path.home() / ".webbscope" / "connections"
    f, tmp = d / f"{host}-{port}.json", d / f".{host}-{port}.{os.getpid()}.tmp"
    try:
        d.mkdir(parents=True, exist_ok=True)
        os.chmod(d, 0o700)
        fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w") as fh:
            json.dump({"url": f"http://127.0.0.1:{port}", "port": port, "host": socket.gethostname(), "agent_key": agent_key(), "pid": os.getpid(),
                       "version": version, "started": time.strftime("%Y-%m-%dT%H:%M:%S%z")}, fh, indent=1)
        os.replace(tmp, f)
    except OSError:
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass
        return None
    _conn_file = f
    atexit.register(remove_connection_file)
    return f


def remove_connection_file():
    """Remove this server's connection file (only if it is still this server's)."""
    global _conn_file
    f, _conn_file = _conn_file, None
    if f is None:
        return
    try:
        if json.loads(f.read_text()).get("pid") == os.getpid():
            f.unlink()
    except (OSError, ValueError):
        pass


def _ok(value):
    return bool(value) and hmac.compare_digest(value.encode("utf-8", "ignore"), TOKEN.encode())


_PAGE = """<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>WebbScope: access</title>
<style>body{margin:0;font:14px/1.5 -apple-system,system-ui,sans-serif;background:#16181d;color:#dfe3ea}#b{max-width:560px;margin:16vh auto;padding:0 20px}
h1{font-size:20px;margin:0 0 8px}p{color:#9fb0cc}code{color:#c9d4e6}.e{color:#ff8a80}form{display:flex;gap:6px;margin-top:14px}
input{flex:1;min-width:0;font:inherit;background:#262a33;color:#dfe3ea;border:1px solid #3a4050;border-radius:4px;padding:6px 10px}
button{font:inherit;background:#2b5fb8;color:#fff;border:1px solid #3a74d6;border-radius:4px;padding:6px 12px;cursor:pointer}</style></head><body><div id="b">
<h1>This WebbScope server needs its access link</h1>%s
<p>Open the link the server printed when it started (<code>http://localhost:&lt;port&gt;/?token=…</code>), or paste its token here. If the server was restarted, it has a new token; a page you already have open keeps working once you open the new link in another tab.</p>
<form method="get"><input name="token" type="password" placeholder="access token" autocomplete="off" autofocus aria-label="access token"><button>Open</button></form>
</div></body></html>"""


def _enter(request):
    """A valid token in the address: set the cookie and go to the same page without it. A page script rather than an HTTP redirect, because a browser
    does not send a SameSite=Strict cookie on a redirect that started from another site (a link clicked in a chat or an email)."""
    rest = [(k, v) for k, v in request.query_params.multi_items() if k != "token"]
    target = "/" + request.url.path.lstrip("/\\") + (("?" + urlencode(rest)) if rest else "")
    js = json.dumps(target).replace("<", "\\u003c")
    r = HTMLResponse(f'<!doctype html><meta charset="utf-8"><title>WebbScope</title><script>location.replace({js})</script>')
    r.set_cookie(COOKIE, TOKEN, httponly=True, samesite="strict", path="/")
    r.headers["Cache-Control"] = "no-store"
    r.headers["Referrer-Policy"] = "no-referrer"
    return r


def gate(request):
    """None if the request may go on, else the response to send instead."""
    if TOKEN is None:
        return None
    given = request.query_params.get("token") if request.method in ("GET", "HEAD") else None
    if given is not None and _ok(given):
        return _enter(request)
    h = request.headers
    auth = h.get("authorization", "")
    bearer = auth[7:].strip() if auth[:7].lower() == "bearer " else ""
    if _ok(bearer) or _ok(h.get("x-webbscope-token", "")) or _ok(request.cookies.get(COOKIE, "")):
        return None
    ak = agent_key()
    if ak and any(hmac.compare_digest(x.encode("utf-8", "ignore"), ak.encode()) for x in (bearer, h.get("x-webbscope-token", "")) if x):
        if request.url.path.startswith("/api/agent/"):
            request.scope["webbscope_agent"] = True
            return None
        return JSONResponse({"detail": "this key is Claude's connector key: it reaches only /api/agent/"}, status_code=403)
    if request.url.path.startswith("/api/") or "text/html" not in h.get("accept", ""):
        return JSONResponse({"detail": "this WebbScope server needs its access token"}, status_code=401, headers={"WWW-Authenticate": "Bearer"})
    msg = '<p class="e">That token is not this server\'s.</p>' if given is not None else ""
    return HTMLResponse(_PAGE % msg, status_code=401, headers={"Cache-Control": "no-store", "Referrer-Policy": "no-referrer"})
