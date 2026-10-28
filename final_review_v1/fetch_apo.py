from pathlib import Path
import urllib.request,time,concurrent.futures
P=Path(__file__).resolve().parent
def one(pid):
    dst=P/'cifs'/(pid+'.cif.gz')
    if dst.exists():return pid,'cached'
    for n in range(3):
        try:
            with urllib.request.urlopen('https://files.rcsb.org/download/'+pid+'.cif.gz',timeout=60) as f:data=f.read()
            assert data[:2]==b'\x1f\x8b'
            dst.write_bytes(data);return pid,len(data)
        except Exception:
            if n==2:raise
            time.sleep(2)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    for r in ex.map(one,(P/'pdb_ids.txt').read_text().split()): print(*r,flush=True)
