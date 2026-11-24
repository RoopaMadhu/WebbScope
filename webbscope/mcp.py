"""WebbScope's MCP connector: lets Claude (Claude Code, Claude Desktop or any MCP client) work with the project open in WebbScope, within the
sharing level the user set for that project in the viewer (About > Claude access; agent.py).

Started by the MCP client and spoken to over stdio (JSON-RPC 2.0, one message per line); it needs only Python's standard library:
    python3 -m webbscope.mcp [--url http://127.0.0.1:PORT --key KEY]          (from the WebbScope folder)
    /Applications/WebbScope.app/Contents/MacOS/WebbScope --mcp [...]            (the Mac app)
Without --url / --key (or WEBBSCOPE_URL / WEBBSCOPE_AGENT_KEY) it uses the newest WebbScope server of this account on this machine
(~/.webbscope/connections/*.json), looked up again at every call, so a restarted server is found. For a server on the HPC reached through an
SSH tunnel, give the tunnel's local address and the key shown in the viewer (About > Claude access > Connect Claude).
The key reaches only WebbScope's /api/agent/ endpoints; the user's access token stays with the user.
"""
import argparse, json, os, sys, time, urllib.error, urllib.request
from pathlib import Path

PROTOCOL = "2025-06-18"
NAME = "webbscope"

INSTRUCTIONS = """WebbScope is a viewer for spatial transcriptomics (cells with x/y positions in tissue sections, gene expression, cell types, tissues,
donors and groups). These tools reach the project the user has open in WebbScope, as far as the user's sharing level for that project allows
(off, structure, summaries, full, run code). Only the user can change the level, in WebbScope (About > Claude access); if a tool is refused,
say which level it needs and let the user decide. Start with webbscope_describe. Donors are the statistical unit: compare groups with one value
per donor (e.g. a donor's mean over its cells), not cells, and say how many donors each group has. Every call is shown to the user in an activity log."""

TOOLS = [
    {"name": "webbscope_describe",
     "description": "What the open WebbScope project holds and what the user lets you see: the sharing level, cell and gene counts, the genes, the fields "
                    "(sample = tissue, donor, group, cell_type, cell_state, meta:<column> ...) with their values (and cell counts from the 'summaries' level), "
                    "the groups the user chose, embeddings, images. Call this first.",
     "inputSchema": {"type": "object", "properties": {}, "additionalProperties": False}},
    {"name": "webbscope_query",
     "description": "Numbers from the open project. kind='counts': cells per combination of the `by` fields (level 'summaries'). kind='expression': mean "
                    "expression and % of cells expressing each gene per combination (level 'summaries'; combinations of fewer than 10 cells show their count only; "
                    "at most 200 genes). kind='cells': single cells with id, tissue, x, y, the chosen `columns` and `genes` (level 'full'; `limit` <= 50000, page "
                    "with `offset`). `filters` keeps only cells whose field has one of the listed values, e.g. {\"group\": [\"B\"], \"cell_type\": [\"Cluster 3\"]}. "
                    "Group by donor (and group) to get one value per donor for statistics.",
     "inputSchema": {"type": "object", "properties": {
         "kind": {"type": "string", "enum": ["counts", "expression", "cells"]},
         "by": {"type": "array", "items": {"type": "string"}, "description": "fields to group by (counts / expression), at most 4"},
         "genes": {"type": "array", "items": {"type": "string"}},
         "columns": {"type": "array", "items": {"type": "string"}, "description": "cells: fields or numeric columns to include"},
         "layer": {"type": "string", "enum": ["logcpx", "counts"], "description": "logcpx (default) = the viewer's log-normalised expression; counts = raw counts"},
         "filters": {"type": "object", "additionalProperties": {"type": "array", "items": {"type": "string"}}},
         "limit": {"type": "integer"}, "offset": {"type": "integer"}},
         "required": ["kind"], "additionalProperties": False}},
    {"name": "webbscope_run_code",
     "description": "Run Python or R code you write on the project's data, with the user's own interpreter, on the machine running WebbScope (level 'run code'). "
                    "The data folder is in the environment variable WEBBSCOPE_DATA: cells.parquet (one row per cell: the cell table plus webbscope_sample, "
                    "webbscope_donor, webbscope_group), genes.txt, and the expression matrix as csr_data.npy / csr_indices.npy / csr_indptr.npy (genes x cells, CSR); "
                    "its README.txt says how to load them and how the viewer normalises. Only packages the user already has are available; nothing is installed. "
                    "Print results; PNG files the code writes in its working folder are returned as pictures (up to 4).",
     "inputSchema": {"type": "object", "properties": {
         "language": {"type": "string", "enum": ["python", "r"]}, "code": {"type": "string"},
         "timeout": {"type": "integer", "description": "seconds, default 300, at most 1800"}},
         "required": ["code"], "additionalProperties": False}},
    {"name": "webbscope_viewer",
     "description": "The WebbScope viewer the user is looking at. action='get': what it shows (tissue, colouring, filters, open windows). action='show' (level "
                    "'structure'): show a tissue (`sample`), colour by a field or a gene (`color_by`) and/or set the viewer's filters (`filters`: viewer fields only, "
                    "e.g. cell_type). action='snapshot' (level 'full'): a picture of the main view with its legend. The user sees what you show.",
     "inputSchema": {"type": "object", "properties": {
         "action": {"type": "string", "enum": ["get", "show", "snapshot"]}, "sample": {"type": "string"}, "color_by": {"type": "string"},
         "filters": {"type": "object", "additionalProperties": {"type": "array", "items": {"type": "string"}}}},
         "required": ["action"], "additionalProperties": False}},
]


class Server:
    def __init__(self, url=None, key=None):
        self.url, self.key = url, key

    def target(self):
        """(url, key): as given, else the newest running WebbScope server of this account on this machine."""
        if self.url and self.key:
            return self.url.rstrip("/"), self.key
        d = Path.home() / ".webbscope" / "connections"
        best = None
        for f in d.glob("*.json") if d.is_dir() else []:
            try:
                c = json.loads(f.read_text())
                os.kill(int(c["pid"]), 0)                       # still running
            except (OSError, ValueError, KeyError):
                continue
            if c.get("agent_key") and (best is None or f.stat().st_mtime > best[0]):
                best = (f.stat().st_mtime, c)
        if best is None:
            raise RuntimeError("no running WebbScope server found on this machine (start WebbScope and open a project; for a server on the HPC "
                               "give the connector --url and --key, see About > Claude access > Connect Claude in the viewer)")
        return (self.url or best[1]["url"]).rstrip("/"), self.key or best[1]["agent_key"]

    def call(self, method, path, body=None, timeout=120):
        url, key = self.target()
        req = urllib.request.Request(url + path, json.dumps(body).encode() if body is not None else None,
                                     {"Content-Type": "application/json", "Authorization": f"Bearer {key}"}, method=method)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read() or b"null")
        except urllib.error.HTTPError as e:
            try:
                detail = json.loads(e.read()).get("detail")
            except Exception:
                detail = None
            raise RuntimeError(f"WebbScope answered {e.code}: {detail or e.reason}")
        except urllib.error.URLError as e:
            raise RuntimeError(f"cannot reach WebbScope at {url} ({e.reason}); is it running (and the SSH tunnel open, for the HPC)?")


def _text(obj):
    return {"type": "text", "text": obj if isinstance(obj, str) else json.dumps(obj, separators=(",", ":"), ensure_ascii=False, default=str)}


def run_tool(srv, name, a):
    """The MCP content list for one tool call."""
    if name == "webbscope_describe":
        st = srv.call("GET", "/api/agent/status")
        if st.get("level", "off") == "off":
            return [_text(st)]
        return [_text(srv.call("GET", "/api/agent/describe"))]
    if name == "webbscope_query":
        kind = a.get("kind")
        if kind in ("counts", "expression"):
            body = {k: a[k] for k in ("by", "genes", "layer", "filters") if k in a}
            return [_text(srv.call("POST", "/api/agent/summary", {"kind": kind, **body}, timeout=600))]
        if kind == "cells":
            body = {k: a[k] for k in ("columns", "genes", "layer", "filters", "limit", "offset") if k in a}
            return [_text(srv.call("POST", "/api/agent/cells", body, timeout=600))]
        raise RuntimeError("kind must be counts, expression or cells")
    if name == "webbscope_run_code":
        t = int(a.get("timeout") or 300)
        r = srv.call("POST", "/api/agent/run", {"language": a.get("language", "python"), "code": a["code"], "timeout": t}, timeout=t + 120)
        pics = r.pop("pictures", [])
        out = [_text(r)]
        out += [{"type": "image", "data": p["png_base64"], "mimeType": "image/png"} for p in pics]
        return out
    if name == "webbscope_viewer":
        act = a.get("action")
        if act == "get":
            return [_text(srv.call("GET", "/api/agent/view"))]
        if act == "show":
            body = {k: a[k] for k in ("sample", "color_by", "filters") if a.get(k)}
            return [_text(srv.call("POST", "/api/agent/view", body))]
        if act == "snapshot":
            r = srv.call("POST", "/api/agent/snapshot", timeout=60)
            return [{"type": "image", "data": r["png_base64"], "mimeType": "image/png"}]
        raise RuntimeError("action must be get, show or snapshot")
    raise RuntimeError(f"unknown tool {name}")


def handle(srv, msg):
    """The response to one JSON-RPC message (None for a notification)."""
    mid, method, params = msg.get("id"), msg.get("method"), msg.get("params") or {}
    if mid is None:
        return None                                             # notifications (initialized, cancelled ...) need no answer
    def ok(result):
        return {"jsonrpc": "2.0", "id": mid, "result": result}
    if method == "initialize":
        from . import __version__
        return ok({"protocolVersion": params.get("protocolVersion") or PROTOCOL, "capabilities": {"tools": {"listChanged": False}},
                   "serverInfo": {"name": NAME, "title": "WebbScope", "version": __version__}, "instructions": INSTRUCTIONS})
    if method == "ping":
        return ok({})
    if method == "tools/list":
        return ok({"tools": TOOLS})
    if method == "tools/call":
        name, args = params.get("name"), params.get("arguments") or {}
        try:
            return ok({"content": run_tool(srv, name, args), "isError": False})
        except Exception as e:                                  # reported to Claude as a tool error, not a protocol error
            return ok({"content": [_text(str(e))], "isError": True})
    return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32601, "message": f"method not found: {method}"}}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="webbscope-mcp", description="WebbScope's MCP connector (stdio)")
    ap.add_argument("--url", default=os.environ.get("WEBBSCOPE_URL"), help="the WebbScope server, e.g. http://localhost:8765 (default: the newest one on this machine)")
    ap.add_argument("--key", default=os.environ.get("WEBBSCOPE_AGENT_KEY"), help="its agent key (viewer: About > Claude access > Connect Claude)")
    a = ap.parse_args(argv)
    srv = Server(a.url, a.key)
    inp = sys.stdin if sys.stdin is not None else os.fdopen(0, "r", encoding="utf-8")       # a windowed (app) build may start without them
    out = sys.stdout if sys.stdout is not None else os.fdopen(1, "w", encoding="utf-8")
    for line in inp:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except ValueError:
            out.write(json.dumps({"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "parse error"}}) + "\n"); out.flush()
            continue
        for m in msg if isinstance(msg, list) else [msg]:
            r = handle(srv, m)
            if r is not None:
                out.write(json.dumps(r) + "\n"); out.flush()


if __name__ == "__main__":
    main()
