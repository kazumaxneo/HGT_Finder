"""Server driver: construct verified full R232 DB, then measure one LKSZ1 run."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proteins", type=Path, required=True)
    p.add_argument("--taxonomy", type=Path, required=True)
    p.add_argument("--genbank", type=Path, required=True)
    p.add_argument("--workspace", type=Path, required=True)
    p.add_argument("--state", type=Path, required=True)
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--expected-genomes", type=int, default=199923)
    p.add_argument("--diamond", default="/usr/bin/diamond")
    a=p.parse_args()
    state={"created_utc":datetime.now(timezone.utc).isoformat(),
           "workspace":str(a.workspace),"expected_genomes":a.expected_genomes}
    def save(status, **extra):
        state.update(status=status, updated_utc=datetime.now(timezone.utc).isoformat(), **extra)
        a.state.parent.mkdir(parents=True,exist_ok=True)
        temp=a.state.with_suffix(".tmp")
        temp.write_text(json.dumps(state,ensure_ascii=False,indent=2)+"\n")
        temp.replace(a.state)
        print(json.dumps(state,ensure_ascii=False),flush=True)
    try:
        if a.workspace.exists() and any(a.workspace.iterdir()):
            raise ValueError("Workspace must be new and empty; existing data are never replaced")
        a.workspace.mkdir(parents=True,exist_ok=True)
        db=a.workspace/"database"
        save("building_database")
        began=time.perf_counter()
        subprocess.run([sys.executable,"-u","-m","gtdb_hgt","build-db",
                        "--proteins",str(a.proteins),"--taxonomy",str(a.taxonomy),
                        "--out",str(db),"--threads",str(a.threads),"--diamond",a.diamond,
                        "--database-label","GTDB R232 all representatives"],check=True)
        setup=time.perf_counter()-began
        manifest=json.loads((db/"database.json").read_text())
        if manifest["genomes"]!=a.expected_genomes or manifest["skipped_files"]:
            raise ValueError(f"Full source verification failed: {manifest['genomes']} genomes; "
                             f"{len(manifest['skipped_files'])} skipped files")
        save("measuring_full_run",database_construction_seconds=setup,database=manifest)
        run=a.workspace/"LKSZ1_full_R232_sankey"
        subprocess.run([sys.executable,"-u","scripts/benchmark_run.py",
                        "--expected-genomes",str(a.expected_genomes),"--",
                        "--genbank",str(a.genbank),"--db",str(db/"proteins.dmnd"),
                        "--taxonomy",str(db/"taxonomy.tsv"),"--query-genome","GCF_040436315.1",
                        "--out",str(run),"--sankey","--diamond",a.diamond,"--threads",str(a.threads),
                        "--database-label","GTDB R232 all representatives"],check=True)
        save("complete",benchmark=str(run/"benchmark.json"),report=str(run/"report.html"),
             sankey=str(run/"sankey.html"))
    except Exception as e:
        save("failed",error=str(e))
        raise


if __name__=="__main__":
    main()
