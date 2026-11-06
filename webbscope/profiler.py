"""Lightweight built-in profiling: what happens while a project opens (every log line with its time and memory, plus the engine set-up steps)
and how long each API endpoint takes afterwards. Always on (a dict update per request); read it at /profile or /api/profile."""
import os, resource, sys, threading, time
from contextlib import contextmanager


def peak_gb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3 if sys.platform == "darwin" else 1024 ** 2)           # bytes on macOS, KiB on Linux


def rss_gb():
    """Current resident memory; None where it cannot be read cheaply (macOS without /proc)."""
    try:
        with open("/proc/self/statm") as f:
            return int(f.read().split()[1]) * os.sysconf("SC_PAGE_SIZE") / 1024 ** 3
    except Exception:
        return None


class Profiler:
    def __init__(self):
        self.lock = threading.Lock()
        self.reset()

    def reset(self):
        with self.lock:
            self.t0 = None
            self.project = None
            self.events = []          # [{t, msg, dt, rss, peak}] while a project opens (t = seconds since the open started, dt = since the previous event)
            self.steps = []           # [{name, seconds, peak}] named steps of the open
            self.requests = {}        # endpoint -> {n, total, max, last, bytes}

    def start_open(self, project):
        with self.lock:
            self.t0, self.project, self.events, self.steps = time.time(), str(project), [], []
        self.mark("open started")

    def mark(self, msg):
        with self.lock:
            if self.t0 is None:
                return
            now = time.time() - self.t0
            prev = self.events[-1]["t"] if self.events else 0.0
            self.events.append({"t": round(now, 3), "dt": round(now - prev, 3), "msg": str(msg), "rss": rss_gb(), "peak": round(peak_gb(), 3)})

    @contextmanager
    def step(self, name):
        t = time.time()
        try:
            yield
        finally:
            with self.lock:
                self.steps.append({"name": name, "seconds": round(time.time() - t, 3), "peak": round(peak_gb(), 3)})

    def request(self, endpoint, seconds, nbytes):
        with self.lock:
            r = self.requests.setdefault(endpoint, {"n": 0, "total": 0.0, "max": 0.0, "last": 0.0, "bytes": 0})
            r["n"] += 1; r["total"] += seconds; r["max"] = max(r["max"], seconds); r["last"] = seconds; r["bytes"] += nbytes

    def report(self):
        with self.lock:
            reqs = [{"endpoint": k, "n": v["n"], "total_ms": round(v["total"] * 1000, 1), "mean_ms": round(v["total"] * 1000 / v["n"], 1),
                     "max_ms": round(v["max"] * 1000, 1), "last_ms": round(v["last"] * 1000, 1), "mb": round(v["bytes"] / 1e6, 2)}
                    for k, v in self.requests.items()]
            return {"project": self.project, "open_seconds": self.events[-1]["t"] if self.events else None, "peak_gb": round(peak_gb(), 3), "rss_gb": rss_gb(),
                    "events": list(self.events), "steps": list(self.steps), "requests": sorted(reqs, key=lambda r: -r["total_ms"])}


PROF = Profiler()
