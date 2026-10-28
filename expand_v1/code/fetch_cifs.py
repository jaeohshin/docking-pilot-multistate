"""확대 cohort의 deposited mmCIF를 한 폴더에 모은다. 기존 캐시에서 복사하고, 없는 것만 RCSB에서 받는다."""
import concurrent.futures, json, shutil, time, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
DST = HERE/'cifs'; DST.mkdir(exist_ok=True)
CACHES = [PILOT.parents[1]/'work/binding_site/cif', PILOT/'run_v1/download_cache', PILOT/'final_review_v1/cifs']

need = sorted({x for c in json.loads((HERE/'candidates.json').read_text()) for x in (c['apo'], c['holo'])})
have = {}
for d in CACHES:
    if d.is_dir():
        for f in d.glob('*.cif.gz'):
            have.setdefault(f.name.split('.')[0].lower(), f)

copied = [p for p in need if p in have and not (DST/f'{p}.cif.gz').exists()
          and (shutil.copy2(have[p], DST/f'{p}.cif.gz') or True)]
missing = [p for p in need if not (DST/f'{p}.cif.gz').exists()]
print(f'need {len(need)}  copied {len(copied)}  to download {len(missing)}', flush=True)


def fetch(pid):
    dst = DST/f'{pid}.cif.gz'
    for n in range(3):
        try:
            with urllib.request.urlopen(f'https://files.rcsb.org/download/{pid}.cif.gz', timeout=60) as f:
                data = f.read()
            assert data[:2] == b'\x1f\x8b', 'not gzip'
            dst.write_bytes(data)
            return pid, len(data)
        except Exception as e:
            if n == 2:
                return pid, f'FAILED {e!r}'
            time.sleep(2)


fails = []
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    for pid, r in ex.map(fetch, missing):
        if isinstance(r, str):
            fails.append(pid); print(pid, r, flush=True)
print('downloaded', len(missing) - len(fails), 'failed', len(fails))
assert not fails, fails
assert all((DST/f'{p}.cif.gz').exists() for p in need), 'missing structures after fetch'
print('cifs ready:', len(need))
