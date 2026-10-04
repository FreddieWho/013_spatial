#!/usr/bin/env python3
"""Full-scope descriptive multiscale and reverse figures, Python backend only."""
from pathlib import Path
import csv,gzip,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002';FIG=OUT/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7,'legend.frameon':False})
COLORS=['#9b9b9b','#8a718b','#477e98'];METHODS=['geometry_depth','B_AXIS','train_selected_GO'];LABELS=['Geometry / depth','B axis','GO combination']


def save(fig,name):
    for ext in ['svg','pdf']:fig.savefig(FIG/f'{name}.{ext}',facecolor='white')
    fig.savefig(FIG/f'{name}.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
    fig.savefig(FIG/f'{name}.png',dpi=160,facecolor='white');plt.close(fig)


def support():
    rows=list(csv.DictReader((OUT/'summary/length_width_support.tsv').open(),delimiter='\t'))
    mat=np.array([int(r['evaluable'])/int(r['attempted']) if int(r['attempted']) else np.nan for r in rows]).reshape(6,4)
    fig,ax=plt.subplots(figsize=(7.2,4));im=ax.imshow(mat,cmap='Blues',vmin=0,vmax=1,aspect='auto',interpolation='nearest')
    for i in range(6):
        for j in range(4):
            r=rows[i*4+j];ax.text(j,i,f"{int(r['evaluable']):,} / {int(r['attempted']):,}",ha='center',va='center',fontsize=7,color='white' if mat[i,j]>.65 else '#222222')
    ax.set_xticks(range(4),['100','250','500','1000']);ax.set_yticks(range(6),['<500','500–1000','1000–2000','2000–4000','4000–8000','≥8000'])
    ax.set_xlabel('Corridor half-width (μm)');ax.set_ylabel('Endpoint path length (μm)')
    ax.set_title('Multi-scale support | evaluable / attempted geometries',fontsize=8,pad=10)
    fig.colorbar(im,ax=ax,fraction=.035,pad=.04,label='Evaluable fraction')
    fig.text(.15,.048,'Both paths and masks; counts are geometric tests, not bridges or independent samples.',fontsize=6.5)
    fig.text(.15,.018,'The 480 disconnected combinations have no path-length bin and are listed separately.',fontsize=6.5)
    fig.subplots_adjust(left=.15,right=.94,bottom=.17,top=.88);save(fig,'multiscale_support')


def reverse():
    values={}
    with gzip.open(OUT/'reverse/summary/paired_section_metrics.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['program'] not in METHODS or not r['auc']:continue
            if r['program']=='B_AXIS' and r['arm']!='full':continue
            key=(r['split'],int(r['mask_um']),r['cohort'],r['section_id'])
            values.setdefault(key,{})[r['program']]=float(r['auc'])
    for split in ['leave_section','leave_cohort']:
        fig,axes=plt.subplots(2,3,figsize=(7.2,5),sharey=True)
        for row,cohort in enumerate(['USZ','GSE175540']):
            for col,radius in enumerate([230,500,1000]):
                ax=axes[row,col];sub=[v for k,v in sorted(values.items()) if k[:3]==(split,radius,cohort)]
                for i,v in enumerate(sub):
                    y=np.array([v.get(m,np.nan) for m in METHODS]);jitter=(i-(len(sub)-1)/2)*.017
                    ax.plot(np.arange(3)+jitter,y,color='#cccccc',lw=.5,zorder=1)
                    for j,color in enumerate(COLORS):ax.scatter(j+jitter,y[j],s=11,c=color,zorder=2,edgecolor='white',linewidth=.2)
                if sub:
                    for j,color in enumerate(COLORS):
                        mean=np.nanmean([v.get(METHODS[j],np.nan) for v in sub]);ax.scatter(j,mean,s=32,c=color,marker='D',edgecolor='black',linewidth=.4,zorder=3)
                else:ax.text(.5,.65,'Not testable\n(no eligible two-class evaluation)',transform=ax.transAxes,ha='center',va='center',fontsize=6.5,color='#666666')
                ax.axhline(.5,ls='--',color='#777777',lw=.6,zorder=0);ax.set_ylim(-.02,1.02);ax.set_xlim(-.35,2.35)
                ax.set_xticks(range(3),['Geometry','B axis','GO comb.'],rotation=20,ha='right',fontsize=6.3)
                ax.set_title(f'{cohort} | mask {radius} μm | n = {len(sub)}',fontsize=7)
                if col==0:ax.set_ylabel('Held-out section AUC')
        handles=[Line2D([],[],color=c,marker='o',lw=0,markersize=4,label=l) for c,l in zip(COLORS,LABELS)]
        fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.055),ncol=3,fontsize=7)
        title='Leave-one-section-out' if split=='leave_section' else 'Leave-one-source-out'
        fig.suptitle(title+' | internal validation, all eligible sections',fontsize=9,y=.98)
        fig.text(.5,.018,'Points: sections; diamonds: section means. No patient-independence or calibrated localization claim.',ha='center',fontsize=6.2)
        fig.subplots_adjust(left=.08,right=.99,top=.9,bottom=.19,wspace=.15,hspace=.43);save(fig,'reverse_'+split)


def all_programs():
    names=np.load(OUT/'scores/usz-KC1.npz')['set_ids'].tolist();lookup={name:i for i,name in enumerate(names)}
    conditions=[(str(radius),split,cohort) for split in ['leave_section','leave_cohort'] for radius in [230,500,1000] for cohort in ['USZ','GSE175540']]
    cond_lookup={key:i for i,key in enumerate(conditions)};arms=['radial','angular','full'];matrix=np.full((3,len(names),12),np.nan)
    for r in csv.DictReader((OUT/'reverse/summary/all_programs.tsv').open(),delimiter='\t'):
        if r['arm'] not in arms or r['program'] not in lookup or not r['mean_auc']:continue
        matrix[arms.index(r['arm']),lookup[r['program']],cond_lookup[r['mask_um'],r['split'],r['cohort']]]=float(r['mean_auc'])
    np.savez_compressed(FIG/'inverse_heatmap_source.npz',auc=matrix,program=np.array(names),conditions=np.array(conditions),arms=np.array(arms))
    with (FIG/'program_row_index.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['panel_row','program','panel']);w.writerows((i+1,name,'GO BP' if i<1691 else 'fixed control') for i,name in enumerate(names))
    cmap=plt.get_cmap('RdBu').copy();cmap.set_bad('#dedede')
    for panel,indices in [('GO',np.arange(1691)),('controls',np.arange(1691,1742))]:
        fig,axes=plt.subplots(1,3,figsize=(7.2,9),sharey=True)
        for a,ax in enumerate(axes):
            im=ax.imshow(matrix[a,indices],cmap=cmap,vmin=0,vmax=1,aspect='auto',interpolation='nearest')
            ax.set_title(arms[a].capitalize()+' information',fontsize=8)
            ax.set_xticks(range(12),[f'{"S" if s=="leave_section" else "C"} {r} {"U" if c=="USZ" else "G"}' for r,s,c in conditions],rotation=90,fontsize=5.5)
            ax.axvline(5.5,color='white',lw=1)
            if panel=='GO':
                ax.set_yticks([0,499,999,1499,1690],['1','500','1000','1500','1691']);axes[0].set_ylabel('GO BP row in frozen panel order')
            else:
                ax.set_yticks(range(51),[names[i].replace('HALLMARK_','') for i in indices],fontsize=4.6)
        cax=fig.add_axes([.30,.063,.55,.012]);fig.colorbar(im,cax=cax,orientation='horizontal',label='Section-mean AUC (gray: not testable)')
        fig.suptitle(('All 1691 GO BP sets' if panel=='GO' else 'All 51 fixed controls')+' | no outcome-based row selection',fontsize=9,y=.975)
        fig.text(.5,.014,'S: section holdout; C: source holdout; U: USZ; G: GSE175540. Column number: mask radius in μm.',ha='center',fontsize=6)
        fig.subplots_adjust(left=.11 if panel=='GO' else .29,right=.99,top=.94,bottom=.15,wspace=.10);save(fig,'inverse_all_'+panel)


if __name__=='__main__':
    assert (FIG/'figure_contract.json').exists()
    support();reverse();all_programs()
    (FIG/'reverse_figures_complete.json').write_text(json.dumps({'status':'RENDERED_PENDING_VISUAL_QA','figures':5,'GO_rows':1691,'control_rows':51,'backend':'Python','formal_p':False},indent=2)+'\n')
