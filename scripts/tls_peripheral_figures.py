#!/usr/bin/env python3
"""R-17 figures. Full data tables accompany every plot; no patient CI implied."""
from pathlib import Path
import csv,json,importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_peripheral_20261002';FIG=OUT/'figures'
FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
def table(path):
 with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
def save(fig,name):
 fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');fig.savefig(FIG/(name+'.svg'),bbox_inches='tight');fig.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)
def short(s):return s.replace('gse175540-','').replace('_frozen','').replace('_ffpe','')

def masks():
 spec=importlib.util.spec_from_file_location('p',ROOT/'scripts/tls_peripheral.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 contract=json.loads((OUT/'contract.json').read_text());elig=table(OUT/'eligibility.tsv')
 for cohort in ['USZ','GSE175540']:
  items=[r for r in contract['sections'] if r['cohort']==cohort];nc=3 if cohort=='GSE175540' else 2;nr=int(np.ceil(len(items)/nc))
  fig,axs=plt.subplots(nr,nc,figsize=(7.2,nr*2),squeeze=False)
  for ax,info in zip(axs.ravel(),items):
   sid=info['section_id'];g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];lab=g['labels'];d,v,n,f=m.geometry(xy,lab,130)
   ax.scatter(*xy.T,s=1,c='#dadada',rasterized=True)
   for mask,color in [(f,'#88b8d0'),(n,'#47778e'),((d>130)&(d<=230),'#e7bd80'),((d>0)&(d<=130),'#c8b09a'),(lab=='TLS','#a75552')]:
    ax.scatter(*xy[mask].T,s=2,c=color,rasterized=True)
   q=next(r for r in elig if r['section_id']==sid and r['mask_um']=='130')
   ax.set_title(short(sid)+'\n'+('eligible' if q['status']=='ELIGIBLE' else q['reason']),fontsize=6.3)
   ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
   x,y=xy.min(axis=0);ax.plot([x,x+1000],[y-150,y-150],c='black',lw=1);ax.text(x,y-300,'1 mm',fontsize=5,va='top');ax.set_ylim(y-700,xy[:,1].max()+200)
  for ax in axs.ravel()[len(items):]:ax.axis('off')
  fig.suptitle(cohort+' | known TLS mask and fixed distance bands\nred: TLS; tan: excluded first ring; amber: extra ring; dark blue: near; pale blue: far',fontsize=8)
  fig.tight_layout(rect=(0,0,1,.96));save(fig,'mask_'+cohort)

def results():
 effects=table(OUT/'effects.tsv');curves=table(OUT/'curves.tsv');panel=json.loads((OUT/'panel.json').read_text())
 # All 51 rows, both masks; cells are observed section effects, not tests.
 sets=list(panel);sections=list(dict.fromkeys(r['section_id'] for r in effects))
 fig,axes=plt.subplots(1,2,figsize=(11,12),sharey=True)
 for ax,rad in zip(axes,[130,230]):
  mat=np.full((len(sets),len(sections)),np.nan)
  for r in effects:
   if int(r['mask_um'])==rad and r['status']=='DESCRIPTIVE_ONLY':mat[sets.index(r['set_id']),sections.index(r['section_id'])]=float(r['effect_sd'])
  im=ax.imshow(mat,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto',interpolation='nearest')
  ax.set_title(f'Mask {rad} um | section-level, descriptive')
  ax.set_xticks(range(len(sections)));ax.set_xticklabels([short(s) for s in sections],rotation=90,fontsize=5.5)
  ax.set_yticks(range(len(sets)));ax.set_yticklabels([s.replace('HALLMARK_','') for s in sets],fontsize=6)
 cax=fig.add_axes([.94,.35,.015,.30]);fig.colorbar(im,cax=cax,label='Near - far / visible SD (clipped at +/-1)')
 fig.subplots_adjust(left=.30,bottom=.19,right=.91,wspace=.10,top=.96);save(fig,'all_effects')
 # Complete curves, one page per readout; no data-driven panel selection.
 with PdfPages(FIG/'all_readout_curves.pdf') as pdf:
  for name in sets:
   fig,axs=plt.subplots(1,2,figsize=(7.2,4.3),sharey=True)
   for ax,cohort in zip(axs,['USZ','GSE175540']):
    sub=[r for r in curves if r['set_id']==name and r['mask_um']=='130' and r['cohort']==cohort and r['readout_status']=='DESCRIPTIVE_ONLY']
    for si,sid in enumerate(dict.fromkeys(r['section_id'] for r in sub)):
     vals=[r for r in sub if r['section_id']==sid and int(r['n_spots'])>=30 and r['mean_score']!='']
     ax.plot([float(r['lo_um'])+50 for r in vals],[float(r['mean_score']) for r in vals],alpha=.8,lw=.8,color=plt.get_cmap('tab20')(si),label=sid.replace('usz-','').replace('gse175540-','').split('_')[0])
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,fontsize=5,frameon=False);ax.set_title(cohort+' | individual sections');ax.set_xlabel('Distance to labelled TLS (um)');ax.set_xlim(100,1000)
   axs[0].set_ylabel('Visible-spot AUCell mean');fig.suptitle(name.replace('HALLMARK_',''),fontsize=8);fig.subplots_adjust(left=.10,right=.98,top=.88,bottom=.39,wspace=.18);pdf.savefig(fig)
   if name=='B_AXIS':fig.savefig(FIG/'B_axis_curves.png',dpi=180,bbox_inches='tight');fig.savefig(FIG/'B_axis_curves.svg',bbox_inches='tight')
   plt.close(fig)

def calibration():
 rows=table(OUT/'calibration.tsv');names=['white','short150','long800','mixture','non_gaussian','random_trend']
 fig,axes=plt.subplots(1,2,figsize=(7.2,3.3))
 for k,name in enumerate(names):
  sub=[r for r in rows if r['kind']=='null' and r['family']==name]
  axes[0].scatter(np.full(len(sub),k)+np.linspace(-.18,.18,len(sub)),[float(r['rate']) for r in sub],s=8,alpha=.5,color='#47778e')
  sub=[r for r in rows if r['kind']=='power' and r['effect_sd']=='0.5' and r['family']==name]
  axes[1].scatter(np.full(len(sub),k)+np.linspace(-.18,.18,len(sub)),[float(r['rate']) for r in sub],s=8,alpha=.4,color='#9b665e')
 axes[0].axhline(.05,color='black',ls='--',lw=.8);axes[0].axhline(.1,color='gray',ls=':',lw=.8);axes[0].set_ylabel('Null rejection fraction');axes[0].set_title('400 null draws / section, mask, family')
 axes[1].axhline(.8,color='black',ls='--',lw=.8);axes[1].set_ylabel('Detection fraction');axes[1].set_title('0.5 SD injected contrast | 200 draws')
 for ax in axes:ax.set_xticks(range(6));ax.set_xticklabels(names,rotation=45,ha='right',fontsize=6);ax.set_ylim(-.01,1.01)
 axes[0].set_ylim(-.005,.22)
 fig.suptitle('Known-model calibration on complete eligible tissue supports',fontsize=8);fig.tight_layout();save(fig,'calibration')

if __name__=='__main__':
 import sys
 {'masks':masks,'results':results,'calibration':calibration}[sys.argv[1]]()
