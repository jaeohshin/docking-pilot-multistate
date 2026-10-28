"""Fetch wwPDB validation reports (XML) for the holo entries of the 79 fixed candidates, for ligand RSCC/RSR.
Public RCSB files only; existing files are kept. Run on the login node (46 small files)."""
import json, time, urllib.request, concurrent.futures
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
C = json.loads((HERE/'candidates.json').read_text())
DST = HERE/'validation_xml'; DST.mkdir(exist_ok=True)


def one(pid):
    dst = DST/f'{pid}_validation.xml.gz'
    if dst.exists():
        return pid, 'cached'
    url = f'https://files.rcsb.org/pub/pdb/validation_reports/{pid[1:3]}/{pid}/{pid}_validation.xml.gz'
    for n in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as f:
                data = f.read()
            assert data[:2] == b'\x1f\x8b', 'not gzip'
            dst.write_bytes(data)
            return pid, len(data)
        except Exception as e:
            if n == 2:
                return pid, f'failed {e!r}'
            time.sleep(2)


with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    for r in ex.map(one, sorted({c['holo'] for c in C})):
        print(*r, flush=True)
