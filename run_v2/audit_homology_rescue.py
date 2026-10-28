import csv,gzip,json,re,hashlib
from pathlib import Path
from collections import defaultdict,Counter
r=Path('/gpfs/jaeohshin/work/sidechain-multistate')
split=r/'runs/2026-10-19_scratch-v2_01/split'
manifest=list(csv.DictReader((split/'manifest.tsv').open(),delimiter='\t'))
certificate=json.loads((split/'verified.json').read_text())
assert hashlib.sha256((split/'manifest.tsv').read_bytes()).hexdigest()==certificate['manifest_sha256']
group_splits=defaultdict(set); entry_splits=defaultdict(set)
for m in manifest:
    group_splits[m['group_id']].add(m['split'])
    pid=m['pdb_id'] or m['stem'].split('_')[0]
    entry_splits[pid.lower()].add(m['split'])
src=(r/'scripts/178_pocket_variation.py').read_text()
cof=set(re.search(r'COFACTORS = set\("""(.*?)"""\.split',src,re.S).group(1).split())
cases={}; pairs=set(); groups=set(); rows=0
with gzip.open(r/'results/apo_holo/residues.csv.gz','rt') as h:
    for d in csv.DictReader(h):
        rows+=1; pairs.add((d['apo'],d['holo'])); groups.add(d['group'])
        if d['pocket']!='1': continue
        for ligand in d['ligands'].split(';'):
            if not ligand: continue
            k=(d['group'],d['apo'],d['holo'],d['chain'],ligand)
            if k not in cases:
                cases[k]={'group':k[0],'apo':k[1],'holo':k[2],'holo_chain':k[3],'ligand_id':k[4],
                          'cofactor':ligand.split(':')[0] in cof,'res_apo':d['res_apo'],'res_holo':d['res_holo'],
                          'residues':set(),'apo_mc_residues':set(),'selected_residues':set(),'changed_residues':set()}
            c=cases[k]
            assert (c['res_apo'],c['res_holo'])==(d['res_apo'],d['res_holo'])
            assert d['i'] not in c['residues'], 'duplicate holo residue/ligand key'
            c['residues'].add(d['i'])
            if d['apo_mc']=='1': c['apo_mc_residues'].add(d['i'])
            if d['cls']=='selected': c['selected_residues'].add(d['i'])
            if d['cls']!='same': c['changed_residues'].add(d['i'])
def quality(c):
    try: return 0<float(c['res_apo'])<=2.5 and 0<float(c['res_holo'])<=2.5
    except ValueError: return False
def n(cs):
    return {'ligand_chain_cases':len(cs),'pdb_pairs':len({(c['apo'],c['holo']) for c in cs}),'identical_sequence_groups':len({c['group'] for c in cs})}
stages={}; cur=list(cases.values()); stages['pocket_ligand_rows']=n(cur)
cur=[c for c in cur if not c['cofactor']]; stages['noncofactor']=n(cur)
cur=[c for c in cur if quality(c)]; stages['resolution_both_le_2_5']=n(cur)
overlap=Counter(tuple(sorted(group_splits[c['group']])) for c in cur)
cur=[c for c in cur if 'train' not in group_splits[c['group']]]; stages['exclude_train_identical_group']=n(cur)
cur=[c for c in cur if 'validation' not in group_splits[c['group']]]; stages['exclude_validation_identical_group']=n(cur)
cur=[c for c in cur if not ({'train','validation'} & (entry_splits[c['apo']]|entry_splits[c['holo']]))]; stages['exclude_trainval_entry_overlap']=n(cur)
strata={}
for name,fn in [('apo_mc',lambda c:bool(c['apo_mc_residues'])),('selected',lambda c:bool(c['selected_residues'])),('changed',lambda c:bool(c['changed_residues'])),('known_exploratory_group',lambda c:'exploratory' in group_splits[c['group']]),('unknown_group',lambda c:not group_splits[c['group']])]:
    strata[name]=n([c for c in cur if fn(c)])

import sys,time,subprocess,os
from Bio import Align
outdir=Path(__file__).resolve().parent
assert not (outdir/'homology.json').exists(), 'preserve previous results'
def dump(name,obj):
    (outdir/name).write_text(json.dumps(obj,indent=2)+'\n')
def fasta_read(path):
    d={}; key=None
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for line in f:
            if line.startswith('>'):
                key=line[1:].split()[0]
                if key in d: raise ValueError('duplicate FASTA '+key)
                d[key]=''
            elif key: d[key]+=line.strip()
    return d
def write_fa(name,d):
    p=outdir/name
    p.write_text(''.join('>'+k+'\n'+s+'\n' for k,s in sorted(d.items())))
    return p
print('metadata',stages,flush=True)
wanted={c[x] for c in cur for x in ('apo','holo')}
entries={}
for p in sorted((r/'reference/frozen/data/qfit_at_scale/replicate_cache').glob('*.jsonl')):
    with p.open() as f:
        for line in f:
            e=json.loads(line)
            if e['pdb_id'] in wanted:
                if e['pdb_id'] in entries and entries[e['pdb_id']]!=e: raise ValueError('conflicting cache')
                entries[e['pdb_id']]=e
keys={pid+'_'+ch for pid,e in entries.items() for ch in e['chains']}
seqs={}; key=None
with gzip.open(r/'data/pdb_seqres.txt.gz','rt') as f:
    for line in f:
        if line.startswith('>'):
            token=line[1:].split()[0]
            key=token[:4].lower()+token[4:]
            if key not in keys or 'mol:protein' not in line: key=None
            if key:
                if key in seqs: raise ValueError('duplicate SEQRES '+key)
                seqs[key]=''
        elif key: seqs[key]+=line.strip()

# Rescue missing full SEQRES only from deposited mmCIF polymer entities; never use observed sequence.
import gemmi,urllib.request
missing_cases=[c for c in cur if c['holo']+'_'+c['holo_chain'] not in seqs]
rescue_pids={c[k] for c in missing_cases for k in ('apo','holo')}
rescue_cache=outdir/'seqres_cache'; rescue_cache.mkdir(exist_ok=True)
rescue=[]
for pid in sorted(rescue_pids):
    try:
        p=ROOT_UNUSED if False else r/'work/binding_site/cif'/(pid+'.cif.gz')
        if not p.exists():
            p=rescue_cache/(pid+'.cif.gz')
            if not p.exists():
                with urllib.request.urlopen('https://files.rcsb.org/download/'+pid+'.cif.gz',timeout=30) as f: data=f.read()
                p.write_bytes(data)
        st=gemmi.read_structure(str(p)); st.setup_entities()
        added=[]
        for ch in st[0]:
            poly=ch.get_polymer()
            if not len(poly): continue
            entity=st.get_entity_of(poly)
            if entity is None or not entity.full_sequence: continue
            s=gemmi.one_letter_code(entity.full_sequence)
            key=pid+'_'+ch.name
            if key in keys and key not in seqs:
                if not s or 'X' in s: continue
                seqs[key]=s;added.append(key)
        rescue.append({'pdb':pid,'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'added':added})
    except Exception as exc:
        rescue.append({'pdb':pid,'error':repr(exc)})
dump('seqres_rescue.json',rescue)
print('full_seqres_rescue',sum(len(x.get('added',[])) for x in rescue),'chains',flush=True)

qseq={}; records=[]
for c in sorted(cur,key=lambda c:(c['group'],c['apo'],c['holo'],c['holo_chain'],c['ligand_id'])):
    c={k:sorted(v) if isinstance(v,set) else v for k,v in c.items()}
    c['case_id']=hashlib.sha256('|'.join(c[k] for k in ('group','apo','holo','holo_chain','ligand_id')).encode()).hexdigest()[:16]
    hk=c['holo']+'_'+c['holo_chain']
    seq=seqs.get(hk)
    ae=entries.get(c['apo'])
    reason=None
    if not seq: reason='missing_holo_full_seqres'
    elif not ae: reason='missing_apo_cache'
    else:
        matches=[ch for ch in ae['chains'] if seqs.get(c['apo']+'_'+ch)==seq]
        if not matches: reason='no_apo_chain_with_identical_full_seqres'
        else:
            # Sequence-equivalent apo copies: prefer longest observed chain, then matching chain name, then lexical.
            matches.sort(key=lambda ch:(-len(ae['chains'][ch]['seq']),ch!=c['holo_chain'],ch))
            c['apo_chain']=matches[0]; c['equivalent_apo_chains']=matches
            qid='Q'+hashlib.sha256(seq.encode()).hexdigest()[:16]
            if qid in qseq: assert qseq[qid]==seq
            qseq[qid]=seq; c['query_id']=qid
    c['sequence_status']=reason or 'full_seqres_pair_verified'
    records.append(c)
dump('metadata_candidates.json',records)
dump('prescreen.json',{'stages':stages,'sequence_status':dict(Counter(c['sequence_status'] for c in records)),
                       'query_sequences':len(qseq),'source_rows':rows})
print('sequence_status',Counter(c['sequence_status'] for c in records),'queries',len(qseq),flush=True)
refs={}
for subset in ('train','validation'):
    for key,s in fasta_read(split/(subset+'.fasta')).items():
        rid=subset+'__'+key
        assert rid not in refs; refs[rid]=s
query=write_fa('candidate_full_seqres.fasta',qseq)
reference=write_fa('train_validation.fasta',refs)
mm='/gpfs/deepfold/others/mmseqs/bin/mmseqs'
commands=[]
def search(name,a,b,ntarget):
    output=outdir/(name+'.tsv')
    cmd=[mm,'easy-search',str(a),str(b),str(output),str(outdir/(name+'_tmp')),
         '--threads','4','-s','7.5','--max-seqs',str(max(2000,ntarget)),'-e','1e-3',
         '--min-seq-id','0.3','--alignment-mode','3','--seq-id-mode','0',
         '--split-memory-limit','4G','--format-output','query,target,fident,alnlen,evalue,bits,qcov,tcov,qlen,tlen',
         '--remove-tmp-files','1']
    t=time.monotonic()
    with (outdir/(name+'.log')).open('w') as f: subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    commands.append({'command':cmd,'seconds':time.monotonic()-t,'output_sha256':hashlib.sha256(output.read_bytes()).hexdigest()})
    hits=[]
    with output.open() as f:
        for line in f:
            a=line.rstrip().split('\t')
            assert len(a)==10
            if float(a[2])>=.30 and float(a[4])<=.001 and max(float(a[6]),float(a[7]))>=.5:
                hits.append(a)
    print(name,'qualifying',len(hits),'seconds',commands[-1]['seconds'],flush=True)
    return hits
hits1=search('candidate_vs_trainval',query,reference,len(refs))
hits2=search('trainval_vs_candidate',reference,query,len(qseq))
excluded={h[0] for h in hits1}|{h[1] for h in hits2}
refseq=set(refs.values())
excluded|={q for q,s in qseq.items() if s in refseq}
clean={q:s for q,s in qseq.items() if q not in excluded}
edges=[]
if clean:
    cleanfa=write_fa('homology_clean.fasta',clean)
    edges=search('clean_all_vs_all',cleanfa,cleanfa,len(clean))
for c in records:
    c['homology_status']='not_screened' if 'query_id' not in c else 'trainval_homology' if c['query_id'] in excluded else 'no_qualifying_trainval_hit'
dump('homology_candidates.json',records)
cleanrecords=[c for c in records if c['homology_status']=='no_qualifying_trainval_hit']
dump('homology.json',{'status':'complete','policy':{'identity':.3,'coverage_either':.5,'evalue':.001,'reciprocal':True},
    'candidate_sequences':len(qseq),'excluded_sequences':len(excluded),'remaining_sequences':len(clean),
    'remaining':n(cleanrecords),'edges':edges,'commands':commands,
    'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [query,reference,split/'manifest.tsv',r/'results/apo_holo/residues.csv.gz']}})
# Only small entry records needed for subsequent residue mapping. No structures copied yet.
needed={c[k] for c in cleanrecords for k in ('apo','holo')}
dump('clean_cache.json',{k:entries[k] for k in needed})
print('HOMOLOGY_DONE',n(cleanrecords),'unique_sequences',len(clean),flush=True)


