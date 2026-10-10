"""Measure an actual gtdb_hgt run including Python startup; no shell execution."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

try:
    import resource
except ImportError:
    resource = None


def option(args, key, default=None):
    for i, arg in enumerate(args):
        if arg == key:
            if i + 1 >= len(args):
                raise ValueError(f"Missing value for {key}")
            return args[i + 1]
        if arg.startswith(key + "="):
            return arg[len(key) + 1:]
    return default


def command_info(command):
    p = subprocess.run(command, capture_output=True, text=True, timeout=30)
    return {"command": command, "returncode": p.returncode,
            "stdout": p.stdout.strip(), "stderr": p.stderr.strip()}


def machine_info():
    model = platform.processor()
    if sys.platform.startswith("linux"):
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    elif sys.platform == "darwin":
        model = command_info(["sysctl", "-n", "machdep.cpu.brand_string"])["stdout"]
    result = {"host": platform.node(), "system": platform.platform(),
            "architecture": platform.machine(), "cpu": model,
            "logical_cpus": os.cpu_count(),
            "load_average_before_run": list(os.getloadavg()) if hasattr(os, "getloadavg") else None}
    if sys.platform.startswith("linux"):
        memory = {}
        for line in Path("/proc/meminfo").read_text().splitlines():
            key, value = line.split(":", 1)
            if key in ("MemTotal", "MemAvailable", "SwapTotal", "SwapFree"):
                memory[key] = int(value.split()[0]) * 1024
        result["memory_before_run_bytes"] = memory
    return result


def source_fingerprints():
    root = Path(__file__).resolve().parent.parent
    paths = sorted(p for p in (root / "gtdb_hgt").rglob("*")
                   if p.is_file() and not p.name.startswith(".")
                   and p.suffix in (".py", ".js"))
    paths.append(Path(__file__).resolve())
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--expected-genomes", type=int,
                   help="Require this genome count in the DB's database.json before measuring")
    p.add_argument("run_args", nargs=argparse.REMAINDER,
                   help="After --, supply the arguments for HGT_Finder run")
    a = p.parse_args()
    args = a.run_args[1:] if a.run_args[:1] == ["--"] else a.run_args
    out_arg, db_arg = option(args, "--out"), option(args, "--db")
    if not out_arg or not db_arg:
        p.error("--out and --db must be present after --")
    out, db = Path(out_arg), Path(db_arg)
    if out.exists() and any(out.iterdir()):
        p.error("Use a fresh output directory for every measured run")
    if not db.exists() and Path(str(db) + ".dmnd").exists():
        db = Path(str(db) + ".dmnd")
    if not db.is_file():
        p.error(f"Database file not found: {db}")
    db_manifest = db.parent / "database.json"
    metadata = json.loads(db_manifest.read_text()) if db_manifest.exists() else None
    if a.expected_genomes is not None and (not metadata or metadata.get("genomes") != a.expected_genomes):
        p.error("Database manifest missing or genome count differs from --expected-genomes")
    diamond = option(args, "--diamond", "diamond")
    info = {"preflight_started_utc": datetime.now(timezone.utc).isoformat(),
            "machine": machine_info(), "python": a.python,
            "source_sha256": source_fingerprints(),
            "diamond_version": command_info([diamond, "version"]),
            "database": {"path": str(db.resolve()), "size_bytes": db.stat().st_size,
                         "metadata": metadata, "dbinfo": command_info([diamond, "dbinfo", "--db", str(db)])},
            "cache_state": "Not controlled; no OS cache flushing",
            "database_construction_included": False}
    if info["database"]["dbinfo"]["returncode"]:
        p.error("DIAMOND could not read the database; see diamond dbinfo")
    command = [a.python, "-m", "gtdb_hgt", "run", *args]
    before = resource.getrusage(resource.RUSAGE_CHILDREN) if resource else None
    info["started_utc"] = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    if hasattr(os, "wait4"):
        child = subprocess.Popen(command)
        _, status, usage = os.wait4(child.pid, 0)
        code = os.waitstatus_to_exitcode(status)
        child.returncode = code
    else:
        code = subprocess.run(command).returncode
        usage = resource.getrusage(resource.RUSAGE_CHILDREN) if resource else None
    info.update(command=command, exit_code=code,
                whole_command_wall_seconds=time.perf_counter() - started,
                finished_utc=datetime.now(timezone.utc).isoformat())
    if usage:
        multiplier = 1 if sys.platform == "darwin" else 1024
        subtract = before if not hasattr(os, "wait4") else None
        info["resources"] = {
            "user_cpu_seconds": usage.ru_utime - (subtract.ru_utime if subtract else 0),
            "system_cpu_seconds": usage.ru_stime - (subtract.ru_stime if subtract else 0),
            "max_process_rss_bytes": usage.ru_maxrss * multiplier,
            "rss_definition": "Largest individual process RSS in the waited command tree; not simultaneous summed RSS"}
    stages = out / "stage_timings.json"
    if stages.exists():
        info["pipeline_timings"] = json.loads(stages.read_text())
    manifest = out / "run.json"
    if code == 0 and manifest.exists():
        info["result"] = json.loads(manifest.read_text())["counts"]
    out.mkdir(parents=True, exist_ok=True)
    with (out / "benchmark.json").open("x") as f:
        json.dump(info, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")
    print(f"Measured command wall time: {info['whole_command_wall_seconds']:.3f} seconds")
    print(f"Benchmark: {(out / 'benchmark.json').resolve()}")
    raise SystemExit(code)


if __name__ == "__main__":
    main()
