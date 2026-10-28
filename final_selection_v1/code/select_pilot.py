"""Final accept/exclude/defer calls for the 79 candidates, one pilot representative per homology component, and the
apo-minor == holo-major counts. Inputs: candidates79.json, assembly_contact_audit.json, ligand_roles.py,
minor_state_residues.tsv. The minor-state signal is joined only after the calls and representatives are fixed.

Case rules, in order (first failing rule decides):
  1 ligand role exclude -> exclude; role defer -> defer                                   (ligand_roles.py)
  2 holo crystal_dependence 'major' in the ligand's assembly -> exclude (crystal-context stratum)
  3 cofactor/co-ligand dependence (>= 10% ligand atoms within 4 A of a cofactor or >= 6-atom HET that is absent at
    the apo site) -> defer (needs an explicit receptor recipe; holo HETs are not transplanted by default)
  4 an in-assembly partner chain contacts the ligand but the apo has no same-protein copy there -> defer
  5 validation RSCC < 0.80 -> defer (reference pose unclear)
  otherwise accept (eligible).
Representative per component among eligible cases, ordered by
  (pose not clear [RSCC >= 0.90 and mean occupancy >= 0.90; missing RSCC = not clear],
   environment tier max(holo, apo crystal dependence) none < minor < major,
   worse of apo/holo resolution, case_id).
The apo side enters the tier only (an apo crystal mate is not in the docking receptor but may shape apo side chains);
this tier definition was fixed after reading the audit summary and before any representative was listed. A holo-only
tier is reported as a sensitivity check.
"""
import csv, hashlib, json, sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_roles import ROLES  # noqa: E402

C = {c['case_id']: c for c in json.loads((PILOT/'final_review_v1/candidates79.json').read_text())}
AUD = json.loads((HERE/'assembly_contact_audit.json').read_text())
A = {o['case_id']: o for o in AUD['cases']}
assert not AUD['errors'] and set(A) == set(C)
TIER = {'none': 0, 'minor': 1, 'major': 2}


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


rows = []
for cid, c in C.items():
    o = A[cid]
    role, verdict, level, evidence = ROLES[(c['holo'], c['ligand_id'])]
    h = o['holo_contacts'][o['holo_assembly_used']]
    a = o['apo_contacts'][o['apo_assembly_used']]
    dep = [d for d in o['holo_neighbour_dependence'] if not d['present_at_apo_site']]
    missing_partner = [p for p in o['holo_partner_chains'] if not p['apo_in_assembly_same_sequence_copies']]
    rscc = float(o['validation']['rscc']) if o['validation'].get('rscc') else None
    occ = o['ligand_mean_occupancy']
    if verdict != 'accept':
        status, reason = verdict, f'ligand role: {role}'
    elif h['crystal_dependence'] == 'major':
        status, reason = 'exclude', (f"crystal-context dependent: {h['outside_share']:.2f} of protein contact atoms and up to "
                                     f"{h['outside_max_ligand_atom_fraction']:.2f} of ligand atoms from non-assembly copies")
    elif dep:
        status, reason = 'defer', 'needs absent ' + ', '.join(f"{'/'.join(d['components'])} ({d['ligand_atom_fraction_4A']:.2f} of ligand atoms <= 4 A)" for d in dep)
    elif missing_partner:
        status, reason = 'defer', 'partner chain absent at the apo site'
    elif rscc is not None and rscc < 0.80:
        status, reason = 'defer', f'reference pose unclear (RSCC {rscc:.2f})'
    else:
        status, reason = 'accept', 'eligible'
    receptor = sorted({f"{r['chain']}{'' if r['copy'] == 'asu' else '(sym op' + str(r['op']) + str(tuple(r['shift'])) + ')'}"
                       for r in o['apo_records'] if r['kind'] == 'Polymer' and r['contact_atoms'] and r['in_assembly'].get(o['apo_assembly_used'])})
    removed = sorted({x for n in o['holo_nonpolymer_neighbours'] for x in n['components']})
    rows.append(dict(
        case_id=cid, homology_component=o['homology_component'], source_group=c['group'], apo=c['apo'], holo=c['holo'],
        apo_chain=c['apo_deposited_chain'], holo_chain=c['holo_deposited_chain'], ligand=c['ligand_id'], ligand_name=c['ccd']['name'],
        heavy_atoms=c['ligand_heavy_atoms'], res_apo=c['res_apo'], res_holo=c['res_holo'], ligand_role=role, role_verdict=verdict,
        evidence_level=level, role_evidence=evidence, rscc=rscc, ligand_mean_occupancy=round(occ, 3),
        holo_assembly=f"{o['holo_assembly_used']}:{o['holo_assemblies'][o['holo_assembly_used']]['details']}",
        holo_crystal_dependence=h['crystal_dependence'], holo_outside_share=round(h['outside_share'], 3),
        holo_outside_ligand_fraction=round(h['outside_max_ligand_atom_fraction'], 3), holo_outside_copies=';'.join(h['outside_copies']),
        apo_crystal_dependence=a['crystal_dependence'], apo_outside_share=round(a['outside_share'], 3),
        apo_outside_ligand_fraction=round(a['outside_max_ligand_atom_fraction'], 3),
        partner_chains=';'.join(f"{p['holo_copy']}->{','.join(p['apo_in_assembly_same_sequence_copies']) or 'MISSING'}" for p in o['holo_partner_chains']),
        cofactor_or_coligand_dependence=';'.join(f"{'/'.join(d['components'])}:{d['ligand_atom_fraction_4A']:.2f}:{'apo' if d['present_at_apo_site'] else 'absent'}" for d in o['holo_neighbour_dependence']),
        apo_receptor_chains=';'.join(receptor), holo_het_neighbours_not_in_receptor=';'.join(removed),
        pocket_backbone_rmsd_deposited=round(o['deposited_pocket_backbone_rmsd'], 3),
        transfer_qfit_vs_deposited_rmsd=round(o['qfit_vs_deposited_transfer_rmsd'], 3),
        status=status, reason=reason))

# representatives, fixed before the minor-state join
def rank(r, holo_only=False):
    clear = r['rscc'] is not None and r['rscc'] >= 0.90 and r['ligand_mean_occupancy'] >= 0.90
    tier = TIER[r['holo_crystal_dependence']] if holo_only else max(TIER[r['holo_crystal_dependence']], TIER[r['apo_crystal_dependence']])
    return (not clear, tier, max(float(r['res_apo']), float(r['res_holo'])), r['case_id'])


eligible = [r for r in rows if r['status'] == 'accept']
by_comp = defaultdict(list)
for r in eligible:
    by_comp[r['homology_component']].append(r)
reps = {k: min(v, key=rank)['case_id'] for k, v in by_comp.items()}
reps_holo_only = {k: min(v, key=lambda r: rank(r, True))['case_id'] for k, v in by_comp.items()}
for r in rows:
    r['pilot_representative'] = int(r['case_id'] in reps.values())
    r['eligible_reserve'] = int(r['status'] == 'accept' and not r['pilot_representative'])

# minor-state join (after selection)
M = list(csv.DictReader(open(HERE/'minor_state_residues.tsv'), delimiter='\t'))
check = json.loads((HERE/'minor_state_check.json').read_text())
assert check['all_match'], 'recomputed minor-state residues differ from the candidate records'
pos = defaultdict(list)
for m in M:
    pos[m['case_id']].append(f"{m['resname']} holo {m['holo_residue']} / apo {m['apo_residue']} (source i={m['source_index_0based']})")
for r in rows:
    r['minor_state_residues'] = '; '.join(pos.get(r['case_id'], []))

fields = list(rows[0])
order = sorted(rows, key=lambda r: (r['homology_component'], {'accept': 0, 'defer': 1, 'exclude': 2}[r['status']], -r['pilot_representative'], r['case_id']))
for name, sel in [('case_review_all79.tsv', order), ('eligible_cases.tsv', [r for r in order if r['status'] == 'accept']),
                  ('pilot_manifest.tsv', [r for r in order if r['pilot_representative']]),
                  ('mechanism_subset.tsv', [r for r in order if r['minor_state_residues']])]:
    with open(HERE/name, 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t'); w.writeheader(); w.writerows(sel)


def count(sel):
    s = [r for r in sel if r['minor_state_residues']]
    return dict(cases=len(sel), components=len({r['homology_component'] for r in sel}),
                minor_state_cases=len(s), minor_state_components=len({r['homology_component'] for r in s}),
                minor_state_residues=sum(len(pos[r['case_id']]) for r in s), minor_state_case_ids=[r['case_id'] for r in s])


summary = dict(
    status=Counter(r['status'] for r in rows),
    reasons=Counter(r['reason'] if r['status'] != 'accept' else 'eligible' for r in rows),
    all_79=count(rows), eligible=count(eligible), pilot_representatives=count([r for r in rows if r['pilot_representative']]),
    representatives=reps, representative_sensitivity_holo_only_tier=dict(
        same=reps == reps_holo_only, changed={k: [reps[k], reps_holo_only[k]] for k in reps if reps[k] != reps_holo_only[k]}),
    components_all=len({r['homology_component'] for r in rows}),
    components_without_eligible=sorted({r['homology_component'] for r in rows} - set(by_comp)),
    inputs_sha256={p: sha(HERE/p if not p.startswith('..') else HERE/p) for p in
                   ['assembly_contact_audit.json', 'minor_state_residues.tsv', 'minor_state_check.json',
                    '../final_review_v1/candidates79.json', 'code/ligand_roles.py', 'code/select_pilot.py']})
(HERE/'selection_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
print(json.dumps({k: v for k, v in summary.items() if k != 'inputs_sha256'}, indent=1))
