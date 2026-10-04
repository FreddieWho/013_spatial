#!/usr/bin/env python3
"""Bounded acceptance checks for the completed R-17 run, not scientific calibration."""
from pathlib import Path
import csv,hashlib,json,shutil
import numpy as np

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_peripheral_20261002'
def table(name):
 with (OUT/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
contract=json.loads((OUT/'contract.json').read_text())
assert sha(OUT/'contract.json')=='0ef70069927b150afae67d0c2e2d16c2dc8839864c67c6d6d22f69f70695ccf6'
assert sha(OUT/'panel.json')==contract['panel_sha256']
input_hashes=json.loads((OUT/'source_hashes.json').read_text())
for p,record in input_hashes.items():assert sha(ROOT/p)==record['sha256'],p
elig=table('eligibility.tsv');valid={(r['section_id'],r['mask_um']) for r in elig if r['status']=='ELIGIBLE'}
assert len(elig)==52 and len(valid)==48
panel=json.loads((OUT/'panel.json').read_text());effects=table('effects.tsv')
expected={(sid,rad,name) for sid,rad in valid for name in panel}
assert len(effects)==2448
assert {(r['section_id'],r['mask_um'],r['set_id']) for r in effects}==expected
assert all(np.isfinite(float(r['effect_sd'])) and r['formal_p']=='' and r['claim']=='NOT_CALIBRATED_PATIENT_ID_UNKNOWN' for r in effects)
assert len(table('curves.tsv'))==24480
intervention=table('hidden_intervention.tsv');assert len(intervention)==24 and all(float(r['max_abs_change'])==0 for r in intervention)
cal=table('calibration.tsv');assert len(cal)==1440
assert len({(r['section_id'],r['mask_um'],r['family'],r['kind'],r['effect_sd'],r['width_um']) for r in cal})==1440
assert all(int(r['n'])==(400 if r['kind']=='null' else 200) for r in cal)
gate=json.loads((OUT/'calibration_gate.json').read_text());assert gate['status']=='FAILED'
assert all(not r['patient_id'] for r in table('identity.tsv'))
assert shutil.disk_usage(ROOT).free>=1.2e12
receipt=dict(status='EXECUTION_COMPLETE_DESCRIPTIVE_ONLY',input_assets_sha_verified=len(input_hashes),
 source_sections=26,eligible_sections=24,eligible_section_masks=48,readouts=51,effect_rows=2448,curve_rows=24480,
 hidden_intervention_sections=24,hidden_max_abs_change=0,tests=(OUT/'tests.txt').read_text().strip(),
 calibration=gate,calibration_cells=1440,formal_p_blank=True,patient_inference='NOT_RUN_IDENTITY_HARD_BLOCK',
 localization='NOT_RUN_GATE_NOT_MET',external_downloads=False,gpu_used=False,storage_free_bytes=shutil.disk_usage(ROOT).free,
 contract_sha256=sha(OUT/'contract.json'),source_code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'scripts/tls_peripheral.py',ROOT/'scripts/tls_peripheral_figures.py',ROOT/'scripts/validate_tls_peripheral.py']})
(OUT/'validation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(receipt['status'],len(input_hashes),'source assets verified; 2448 effects; 1440 complete calibration cells')
