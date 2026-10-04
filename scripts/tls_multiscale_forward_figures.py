#!/usr/bin/env python3
"""Full program atlas and explicitly exploratory maps, Python backend only."""
from pathlib import Path
import argparse,csv,importlib.util,json,textwrap
import numpy as np
from scipy import sparse,spatial
from scipy.sparse.csgraph import dijkstra
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002';FIG=OUT/'figures';ATLAS=FIG/'forward_atlas'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'legend.frameon':False})
def tab(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def short(sid):return sid.replace('usz-','').replace('gse175540-','').split('_')[0]
def save(fig,path,dpi=300):
    fig.savefig(path.with_suffix('.pdf'),dpi=dpi);fig.savefig(path.with_suffix('.svg'),dpi=dpi)
    fig.savefig(path.with_suffix('.png'),dpi=160);fig.savefig(path.with_suffix('.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)


def atlas_data():
    z=np.load(OUT/'integrated/halo_edge_bins.npz');rows=tab(OUT/'integrated/all_bins.tsv');sources=json.loads((OUT/'contract.json').read_text())['sources'];index={}
    for i,r in enumerate(rows):index.setdefault((r['family'],r['cohort'],r.get('compartment',''),r['section_id'],int(r['mask_um'])),[]).append(i)
    colors={}
    for cohort in ['USZ','GSE175540']:
        sids=[r['section_id'] for r in sources if r['cohort']==cohort]
        colors.update({sid:plt.get_cmap('tab20')(i) for i,sid in enumerate(sids)})
    return z,rows,sources,index,colors


def atlas_page(j,data):
    z,rows,sources,index,colors=data;values=z['means'];name=str(z['set_ids'][j]);fig,axes=plt.subplots(2,2,figsize=(7.2,8.4))
    groups=[('halo','USZ',''),('halo','GSE175540',''),('edge','USZ','TUM'),('edge','USZ','EXPLICIT_NON_TUMOR')]
    titles=['TLS distance | USZ','TLS distance | GSE175540','Labelled interface | TUM','Labelled interface | explicit non-TUM']
    for ax,(family,cohort,comp),title in zip(axes.ravel(),groups,titles):
        total=0
        for source in sources:
            sid=source['section_id']
            if source['cohort']!=cohort:continue
            label_used=False
            for radius,style in [(130,'-'),(230,'--')]:
                ids=index.get((family,cohort,comp,sid,radius),[])
                if not ids:continue
                ids=sorted(ids,key=lambda i:float(rows[i]['lo_um']));x=np.array([(float(rows[i]['lo_um'])+float(rows[i]['hi_um']))/2 for i in ids]);y=values[ids,j]
                if not np.isfinite(y).any():continue
                total+=np.isfinite(y).sum();ax.plot(x,y,style,color=colors[sid],lw=.7,marker='o',markersize=1.6,alpha=.85,label=short(sid) if not label_used else None);label_used=True
        xmax=max(11000,max((float(r['lo_um'])+float(r['hi_um']))/2 for r in rows)*1.05)
        ax.set_xscale('log',base=2);ax.set_xlim(80,xmax);ax.set_xticks([100,400,1600,6400],['100','400','1600','6400'],fontsize=6)
        ax.set_title(title,fontsize=8);ax.set_xlabel('Distance (μm; log axis)',fontsize=6.5);ax.set_ylabel('Mean visible-spot AUCell',fontsize=6.5)
        ax.tick_params(axis='y',labelsize=6);ax.yaxis.offsetText.set_fontsize(6)
        if total:ax.legend(loc='upper center',bbox_to_anchor=(.5,-.28),ncol=4 if cohort=='GSE175540' else 4,fontsize=4.8,handlelength=1.5,columnspacing=.65,labelspacing=.35)
        else:ax.text(.5,.5,'No evaluable distance bins',transform=ax.transAxes,ha='center',color='#777777')
    # Use common y limits within each pair; levels across programs are not compared.
    for axpair in axes:
        bounds=[ax.get_ylim() for ax in axpair if ax.has_data()]
        if not bounds:continue
        lo=min(b[0] for b in bounds);hi=max(b[1] for b in bounds)
        for ax in axpair:ax.set_ylim(lo,hi)
    fig.suptitle('\n'.join(textwrap.wrap(name,88)),fontsize=8,y=.982)
    fig.text(.5,.035,'Solid: 130 μm core margin; dashed: 230 μm. Individual sections; unsupported bins break lines.',ha='center',fontsize=6)
    fig.text(.5,.014,'Available labelled interface fragments only. No fitted field, spatial significance or patient CI.',ha='center',fontsize=6)
    fig.subplots_adjust(left=.09,right=.985,bottom=.17,top=.89,wspace=.29,hspace=1.05)
    return fig


def atlas_book(task):
    start,end,name=task;ATLAS.mkdir(parents=True,exist_ok=True);receipt=ATLAS/f'{name}.json'
    if receipt.exists():return json.loads(receipt.read_text())
    data=atlas_data();index=[]
    with PdfPages(ATLAS/f'{name}.pdf') as pdf:
        for page,j in enumerate(range(start,end),1):
            fig=atlas_page(j,data);pdf.savefig(fig);plt.close(fig)
            index.append(dict(program_index=j,program=str(data[0]['set_ids'][j]),book=f'{name}.pdf',page=page,panel='GO_BP' if j<1691 else 'FIXED_CONTROL'))
    receipt.write_text(json.dumps(dict(status='COMPLETE',start=start,end=end,pages=end-start,index=index),indent=2)+'\n')
    print('atlas book',name,end-start,'pages',flush=True);return json.loads(receipt.read_text())


def atlas(workers):
    tasks=[(a,min(a+128,1691),f'GO_{a//128+1:02}') for a in range(0,1691,128)]+[(1691,1742,'fixed_controls')]
    if workers==1:reports=[atlas_book(t) for t in tasks]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as pool:reports=list(pool.map(atlas_book,tasks))
    index=[row for report in reports for row in report['index']];assert sorted(r['program_index'] for r in index)==list(range(1742))
    with (ATLAS/'index.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(index[0]),delimiter='\t');w.writeheader();w.writerows(index)
    data=atlas_data();save(atlas_page(1741,data),ATLAS/'B_AXIS_preview')
    (ATLAS/'complete.json').write_text(json.dumps(dict(status='ALL1742_CURVE_PAGES_RENDERED',books=len(reports),pages=len(index),GO_pages=1691,control_pages=51),indent=2)+'\n')


def overview():
    d=OUT/'integrated';z=np.load(d/'bridge_all_programs.npz');names=z['set_ids'];groups=tab(d/'bridge_groups.tsv');reports=[]
    units={'raw_median':'Visible-field SD','adjusted_median':'Visible-field SD','native_z_median':'Native-reference SD (not a p-value)',
           'positive_continuity':'Median within-TLS fraction','negative_continuity':'Median within-TLS fraction','side_asymmetry':'Visible-field SD'}
    for metric,title,limit,cmap in [('raw_median','Raw centre − sides',1,'RdBu_r'),('adjusted_median','After endpoint halo baseline',1,'RdBu_r'),('native_z_median','Relative to qualified native references',3,'RdBu_r'),('positive_continuity','Positive ≥4/5 segments',1,'Blues'),('negative_continuity','Negative ≥4/5 segments',1,'Blues'),('side_asymmetry','Absolute left − right contrast',1,'Blues')]:
        matrix=z[metric].T;cm=plt.get_cmap(cmap).copy();cm.set_bad('#dedede');diverging=metric in ['raw_median','adjusted_median','native_z_median']
        fig,ax=plt.subplots(figsize=(9,9));im=ax.imshow(matrix,aspect='auto',interpolation='nearest',cmap=cm,vmin=-limit if diverging else 0,vmax=limit)
        ax.axhline(1690.5,color='black',lw=.6);ax.set_yticks([0,499,999,1499,1690,1741],['1','500','1000','1500','1691','51 controls'])
        ax.set_xticks([0,249,499,749,999,1249,len(groups)-1],[str(i) for i in [1,250,500,750,1000,1250,len(groups)]])
        ax.set_xlabel('Complete geometry-group index (source / family / mask / width / length / path)');ax.set_ylabel('Frozen GO panel order; fixed controls below black line')
        ax.set_title(title+' | all programs and all groups',fontsize=10)
        cax=fig.add_axes([.25,.09,.55,.014]);fig.colorbar(im,cax=cax,orientation='horizontal',extend='both' if diverging else 'max',label=f'{units[metric]}; display clipped at {"±" if diverging else ""}{limit}; gray = not testable')
        fig.text(.5,.026,'Not calibrated spatial evidence. Broad patches and all tested scales remain in the indexed matrix.',ha='center',fontsize=7)
        fig.subplots_adjust(left=.105,right=.985,bottom=.18,top=.95);save(fig,FIG/f'bridge_all_{metric}',dpi=300)
        reports.append(dict(metric=metric,programs=len(names),groups=len(groups),status='DESCRIPTIVE_ONLY',units=units[metric],display_limit=limit))
    (FIG/'bridge_overviews_complete.json').write_text(json.dumps(reports,indent=2)+'\n')


def map_geometry(sid):
    spec=importlib.util.spec_from_file_location('f',ROOT/'scripts/tls_multiscale_forward.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];d=OUT/'forward'/sid;receipt=json.loads((d/'operator_complete.json').read_text());paths=[]
    if not receipt['evaluable_geometries']:return g,paths
    rows=tab(d/'effect_rows.tsv');pairs={(r['pair_id'],r['mode']):r for r in tab(d/'pairs.tsv')};endpoints={r['endpoint_id']:r for r in tab(OUT/'geometry'/f'{sid}_endpoints.tsv')}
    graph=sparse.load_npz(OUT/'geometry'/f'{sid}_graph.npz')
    for family in ['TLS_TLS','TLS_REGION','TLS_EDGE']:
        eligible=sorted([r for r in rows if r['family']==family and r['mask_um']=='130' and r['width_um']=='250' and r['mode']=='tissue_graph_shortest'],key=lambda r:(float(r['length_um']),int(r['pair_id'])))
        if not eligible:continue
        row=eligible[len(eligible)//2];pair=pairs[row['pair_id'],row['mode']];a=int(endpoints[pair['source_endpoint']]['spot_index']);b=int(endpoints[pair['target_endpoint']]['spot_index'])
        _,pred=dijkstra(graph,directed=False,indices=a,return_predecessors=True);path=m.trace_path(pred,a,b)
        assert path is not None
        paths.append(dict(family=family,xy=xy[path],effect_row=int(row['effect_row']),pair_id=row['pair_id'],length_um=float(row['length_um'])))
    return g,paths


def maps():
    representatives=json.loads((OUT/'integrated/representatives.json').read_text())['programs'];sources=json.loads((OUT/'contract.json').read_text())['sources'];cache={};manifest=[]
    qualification={r['section_id']:r['geometry_status'] for r in tab(OUT/'sections.tsv')}
    for source in sources:
        sid=source['section_id'];g,paths=map_geometry(sid);z=np.load(OUT/'scores'/f'{sid}.npz');names=z['set_ids'].tolist();program_indices=[names.index(name) for name in representatives];S=z['scores'][:,program_indices]
        xy=g['xy'];labels=g['labels'];visible=spatial.cKDTree(xy[labels=='TLS']).query(xy)[0]>130
        mean=np.nanmean(S[visible],axis=0);sd=np.nanstd(S[visible],axis=0,ddof=1);values=np.divide(S-mean,sd,out=np.full_like(S,np.nan),where=sd>1e-12)
        values[~visible]=np.nan;values[:,~z['valid'][program_indices]]=np.nan
        if qualification[sid]!='SOURCE_GEOMETRY_AVAILABLE':values[:]=np.nan
        cache[sid]=(g,paths,values)
        for path in paths:manifest.append(dict(section_id=sid,**{k:v for k,v in path.items() if k!='xy'},mask_um=130,halfwidth_um=250,mode='tissue_graph_shortest',selection='median length within family; no expression selection'))
    with (FIG/'representative_candidate_geometry.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
    styles={'TLS_TLS':':','TLS_REGION':'--','TLS_EDGE':'-.'};pages=[]
    for j,name in enumerate(representatives):
        for cohort in ['USZ','GSE175540']:
            sids=[r['section_id'] for r in sources if r['cohort']==cohort];nc=4 if cohort=='USZ' else 3;nr=int(np.ceil(len(sids)/nc))
            fig,axes=plt.subplots(nr,nc,figsize=(7.2,nr*1.9+1.25),squeeze=False)
            for ax,sid in zip(axes.ravel(),sids):
                g,paths,values=cache[sid];xy=g['xy'];labels=g['labels'];ax.scatter(*xy.T,s=.6,c='#e7e7e7',rasterized=True)
                ok=np.isfinite(values[:,j]);im=ax.scatter(*xy[ok].T,c=values[ok,j],s=1.4,cmap='RdBu_r',vmin=-2.5,vmax=2.5,rasterized=True)
                ax.scatter(*xy[labels=='TLS'].T,s=3.5,facecolors='none',edgecolors='black',linewidths=.3,rasterized=True)
                for path in paths:ax.plot(*path['xy'].T,ls=styles[path['family']],color='#222222',lw=.6)
                origin=xy.min(axis=0);ax.plot([origin[0],origin[0]+1000],[origin[1]-150]*2,c='black',lw=.8);ax.text(origin[0],origin[1]-320,'1 mm',fontsize=4.5,va='top')
                ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_ylim(origin[1]-650,xy[:,1].max()+150);ax.set_title(short(sid),fontsize=6.5)
                if not ok.any():ax.text(.5,.5,'Mask / program\nnot evaluable',ha='center',va='center',transform=ax.transAxes,fontsize=5.5)
                for spine in ax.spines.values():spine.set_visible(False)
            for ax in axes.ravel()[len(sids):]:ax.axis('off')
            fig.suptitle('\n'.join(textwrap.wrap(str(name),86))+'\n'+cohort+' | visible molecular scores; candidate geometry only',fontsize=8,y=.985)
            height=fig.get_figheight()
            cax=fig.add_axes([.24,.60/height,.53,.06/height]);fig.colorbar(im,cax=cax,orientation='horizontal',extend='both',label='Visible-spot score z (display clipped at ±2.5)')
            fig.text(.5,.12/height,'Known TLS cores +130 μm hidden; circles: TLS labels. Lines: metadata-selected candidates, NOT confirmed bridges.',ha='center',fontsize=5.8)
            handles=[Line2D([],[],c='#222222',ls=style,lw=.7,label=family.replace('_','–')) for family,style in styles.items()]
            fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.84/height),ncol=3,fontsize=6)
            fig.subplots_adjust(left=.02,right=.99,top=.86 if nr==2 else .93,bottom=1.10/height,wspace=.05,hspace=.13)
            filename=f'representative_{j+1:02}_{cohort}';save(fig,FIG/filename);pages.append(dict(program=name,cohort=cohort,file=filename+'.pdf',sections=sids,section_status={sid:qualification[sid] if np.isfinite(cache[sid][2][:,j]).any() else 'NOT_EVALUABLE' for sid in sids},exploratory=True))
    (FIG/'representative_maps_complete.json').write_text(json.dumps(pages,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['atlas','overview','maps']);ap.add_argument('--workers',type=int,default=1);a=ap.parse_args()
    assert (FIG/'forward_figure_contract.json').exists()
    {'atlas':lambda:atlas(a.workers),'overview':overview,'maps':maps}[a.stage]()
