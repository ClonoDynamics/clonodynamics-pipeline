#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate one semi-synthetic TCR repertoire series on fixed empirical support.

Purpose and workflow position
-----------------------------
Apply one oracle-calibrated combination of sigma_fast and q_bio to the frozen
empirical design. run_dose_series.sh can call this program for all calibrated
doses. This script generates counts and separate truth files; it does not run
production state inference, transition assembly or temporal-scaling analysis.

Required inputs
---------------
--calibration-dir contains:
  00_manifest.json
  02_analysis_sampling_design.csv
  03_step2_pair_parameters.csv
  04_observation_structure_targets.json
  11_subject_union_references.csv
  12_source_repertoire_audit.csv
  the subject reference files indexed by 11_subject_union_references.csv.
The corrected original repertoire files must remain accessible at their recorded
paths, or via --source-dir. Their aaSeqCDR3 IDs define positive support; these
source tables must therefore contain observed positive clonotypes, not padded
zero-count rows. Original count vectors are not copied into simulated outputs.

Generation
----------
For each subject, normalize f_reference_mean_repertoire and generate all nominal
visits. log weights equal log(reference frequency) + sigma_fast*z_fast(t) +
sqrt(q_bio)*W(t); W starts at zero and has independent increments with variance
delta_t. Log-softmax converts weights into normalized frequencies. Biological
fast shocks and unit random-walk paths are fixed by seed, subject and time across
doses. Equal count-RNG seeds do not guarantee identical count allocations or
identical subsequent RNG states across different doses.
For each eligible subject/visit/replicate, retain the empirical positive set A
and total read depth N. Draw independent Gamma(kappa, scale=1/kappa) factors for
active clones, normalize weights f_true*Gamma, allocate N-|A| remaining reads
multinomially, and add one mandatory read to each active clone. Inactive clones
receive zero counts. kappa is the pair-specific empirical Step-2 k parameter.
The stored positivity mask, richness and depth are exact. This is conditional
Gamma-multinomial allocation with mandatory positive reads, not unconstrained
independent negative-binomial sampling.

Scenarios
---------
--scenario null: sigma_fast=0, q_bio=0 (no simulated frequency variation).
--scenario plateau: sigma_fast=0.20, q_bio=0 (non-accumulating variation).
--scenario accumulation: sigma_fast=0.20, q_bio=0.01 (uncalibrated preset).
--scenario custom: explicitly pass the oracle-selected parameters.
Default seed: 20260918. R0p00 in the publication dose series is the PLATEAU,
not the script's strict null preset.

Outputs under --out-dir
-----------------------
repertoires/<subject>_<time>-<replica>.parquet:
  aaSeqCDR3, readCount and readFraction; production-analysis inputs.
ground_truth/subject_<subject>.parquet:
  reference_frequency and logf_bio_t<time>; NEVER pass these to production Step 2.
00_manifest.json: parameters, seeds, hashes, support checks and scope.
01_generated_sampling_summary.csv: one row per generated repertoire.
02_oracle_pair_summary.csv: full-reference interval means/variances (ddof=0),
  including a flag identifying analysis-eligible visit pairs.
03_support_audit_by_subject.csv: observed-support state/transition counts.

Checks and limitations
----------------------
The program verifies exact positivity support, total depth and richness, and
aggregate union/BOTH/common4 support against the empirical targets. Counts are
integer units with no artificial even-count scaling. Only eligible visits are
written as repertoires; all nominal visits have simulated truth. This does not
model the emergence of occupancy/detectability or establish multi-seed power.
An existing output directory is refused. Dependencies: Python >=3.9, NumPy,
pandas and a pandas-compatible Parquet engine (e.g. pyarrow).

Example
-------
python 03_generate_support_conditioned_repertoires.py \
    --calibration-dir CALIBRATION_DIR --source-dir SOURCE_REPERTOIRES \
    --out-dir NEW_SCENARIO_DIR --scenario custom \
    --sigma-fast 0.20 --q-bio 0.01142857142857143 --seed 20260918

Provenance
----------
Documentation-only edition of 7b-generate_empirical_support_conditioned_v5_1.py.
Numerical code, defaults, output names and internal versions are unchanged.
New executions record the new script hash; historical manifests stay unchanged.
"""
from __future__ import annotations
import argparse, hashlib, json, math, platform
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd

VERSION="5.1.0-empirical-support-conditioned-common-random-numbers-2026-09-19"
SCHEMA="clonodynamics_empirical_support_conditioned_v5_1"
SUFFIXES=(".tsv",".tsv.gz",".csv",".csv.gz",".txt",".txt.gz",".parquet",".pq",".feather",".arrow")

def require(c,m):
    if not c: raise ValueError(m)

def sha(path,bs=1024*1024):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(bs),b""): h.update(b)
    return h.hexdigest()

def clean(x):
    if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [clean(v) for v in x]
    if isinstance(x,Path): return str(x)
    if isinstance(x,(np.integer,)): return int(x)
    if isinstance(x,(np.bool_,)): return bool(x)
    if isinstance(x,(float,np.floating)): return float(x) if math.isfinite(float(x)) else None
    return x

def write_json(path,obj):
    Path(path).write_text(json.dumps(clean(obj),indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")

def stable_log_softmax(x):
    x=np.asarray(x,np.float64); m=float(np.max(x))
    return x-(m+math.log(float(np.exp(x-m).sum(dtype=np.float64))))

def parse_bool(s,label):
    if s.dtype==bool: return s.copy()
    m={"true":True,"false":False,"1":True,"0":False,"1.0":True,"0.0":False}
    x=s.astype(str).str.strip().str.lower().map(m)
    require(x.notna().all(),f"Unrecognized boolean in {label}")
    return x.astype(bool)

def strip_suffix(name):
    x=Path(name).name
    changed=True
    while changed:
        changed=False; low=x.lower()
        for suf in (".gz",".tsv",".csv",".txt",".parquet",".pq",".feather",".arrow"):
            if low.endswith(suf):
                x=x[:-len(suf)]; changed=True; break
    return x

class Resolver:
    def __init__(self,root=None):
        self.root=Path(root).resolve() if root else None
        self.by_name={}; self.by_stem={}
        if self.root is None: return
        require(self.root.is_dir(),f"source-dir not found: {self.root}")
        for p in self.root.rglob("*"):
            if not p.is_file() or not any(p.name.lower().endswith(s) for s in SUFFIXES): continue
            self.by_name.setdefault(p.name,[]).append(p.resolve())
            self.by_stem.setdefault(strip_suffix(p.name),[]).append(p.resolve())
    def resolve(self,recorded,requested):
        p=Path(str(recorded))
        if p.is_file(): return p.resolve()
        require(self.root is not None,f"Recorded source missing: {recorded}; supply --source-dir")
        name=Path(str(requested)).name
        hits=self.by_name.get(name,[])
        if len(hits)==1: return hits[0]
        hits=self.by_stem.get(strip_suffix(name),[])
        if len(hits)==1: return hits[0]
        if not hits: raise FileNotFoundError(requested)
        raise RuntimeError(f"Ambiguous source {requested}: {hits}")

def read_ids(path):
    path=Path(path); name=path.name.lower()
    if name.endswith((".parquet",".pq")):
        d=pd.read_parquet(path,columns=["aaSeqCDR3"])
    elif name.endswith((".feather",".arrow")):
        d=pd.read_feather(path,columns=["aaSeqCDR3"])
    else:
        sep="\t" if name.endswith((".tsv",".tsv.gz",".txt",".txt.gz")) else ","
        d=pd.read_csv(path,sep=sep,usecols=["aaSeqCDR3"],keep_default_na=False,dtype={"aaSeqCDR3":str})
    require(not d.empty,f"Empty repertoire: {path}")
    ids=d["aaSeqCDR3"].astype(str)
    require(((ids.str.len()>0)&(ids==ids.str.strip())).all(),f"Invalid aaSeqCDR3: {path}")
    return pd.Index(pd.unique(ids),dtype=object)

def scenario_params(args):
    p={
        "null":{"sigma_fast":0.0,"q_bio":0.0},
        "plateau":{"sigma_fast":0.20,"q_bio":0.0},
        "accumulation":{"sigma_fast":0.20,"q_bio":0.01},
        "custom":{"sigma_fast":0.0,"q_bio":0.0},
    }[args.scenario].copy()
    if args.sigma_fast is not None: p["sigma_fast"]=float(args.sigma_fast)
    if args.q_bio is not None: p["q_bio"]=float(args.q_bio)
    for k,v in p.items(): require(np.isfinite(v) and v>=0,f"{k} invalid")
    return p

def fast_visit_seed(master,subject,time):
    return np.random.SeedSequence([int(master),int(subject),117031,501,int(time)])

def rw_increment_seed(master,subject,time):
    return np.random.SeedSequence([int(master),int(subject),117031,502,int(time)])

def count_seed(master,subject):
    return np.random.SeedSequence([int(master),int(subject),910247,500])

def generate_bio(logf,times,sigma_fast,q_bio,master_seed,subject):
    """Generate truth with common random numbers across all q_bio values.

    z_fast(t) and the unit-Brownian increments are determined only by
    master_seed, subject and time. Therefore changing q_bio rescales the same
    RW path by sqrt(q_bio), and changing sigma_fast rescales the same fast
    shocks. In particular q_bio=0 and q_bio>0 share identical fast shocks.
    """
    out=np.empty((len(times),len(logf)),np.float64)
    rw_unit=np.zeros(len(logf),np.float64)
    prev=times[0]
    sqrt_q=math.sqrt(float(q_bio))
    for j,t in enumerate(times):
        if j>0:
            dt=float(t-prev); require(dt>0,"times not increasing")
            zr=np.random.default_rng(rw_increment_seed(master_seed,subject,t)).normal(0.0,1.0,size=len(logf))
            rw_unit += math.sqrt(dt)*zr
            prev=t
        if sigma_fast>0:
            zf=np.random.default_rng(fast_visit_seed(master_seed,subject,t)).normal(0.0,1.0,size=len(logf))
            fast=float(sigma_fast)*zf
        else:
            fast=0.0
        out[j]=stable_log_softmax(logf + sqrt_q*rw_unit + fast)
    return out

def support_counts(rng,p,active,depth,kappa):
    idx=np.flatnonzero(active); n=len(idx)
    require(n>0 and n<=depth,f"Invalid support n={n}, depth={depth}")
    require(np.isfinite(kappa) and kappa>0,f"invalid kappa={kappa}")
    w=np.asarray(p[idx],np.float64)*rng.gamma(shape=float(kappa),scale=1.0/float(kappa),size=n)
    sw=float(w.sum(dtype=np.float64))
    if not np.isfinite(sw) or sw<=0:
        w=np.asarray(p[idx],np.float64); sw=float(w.sum(dtype=np.float64))
    w/=sw
    extra=depth-n
    extras=rng.multinomial(extra,w).astype(np.int64,copy=False) if extra>0 else np.zeros(n,np.int64)
    c=np.zeros(len(p),np.int64); c[idx]=1+extras
    require(int(c.sum(dtype=np.int64))==depth,"depth audit failed")
    require(np.array_equal(c>0,active),"support audit failed")
    return c

def comb2(x):
    x=np.asarray(x,np.int64); return x*(x-1)//2

def oracle_rows(subject,times,bio,observed):
    rows=[]
    for a in range(len(times)-1):
        for b in range(a+1,len(times)):
            dx=bio[b]-bio[a]
            rows.append(dict(
                subject=subject,time0=int(times[a]),time1=int(times[b]),
                lag=int(times[b]-times[a]),
                both_visits_analysis_eligible=bool(times[a] in observed and times[b] in observed),
                n_reference_clonotypes=int(len(dx)),
                mean_dx_bio=float(np.mean(dx)),
                var_dx_bio=float(np.var(dx,ddof=0)),
            ))
    return rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calibration-dir",required=True,type=Path)
    ap.add_argument("--out-dir",required=True,type=Path)
    ap.add_argument("--scenario",required=True,choices=["null","plateau","accumulation","custom"])
    ap.add_argument("--seed",type=int,default=20260918)
    ap.add_argument("--sigma-fast",type=float,default=None)
    ap.add_argument("--q-bio",type=float,default=None)
    ap.add_argument("--source-dir",type=Path,default=None)
    args=ap.parse_args()

    cal=args.calibration_dir.resolve(); out=args.out_dir.resolve()
    if out.exists(): raise FileExistsError(f"Refusing overwrite: {out}")
    out.mkdir(parents=True); (out/"repertoires").mkdir(); (out/"ground_truth").mkdir()
    pars=scenario_params(args)

    design=pd.read_csv(cal/"02_analysis_sampling_design.csv")
    pair_params=pd.read_csv(cal/"03_step2_pair_parameters.csv")
    refs=pd.read_csv(cal/"11_subject_union_references.csv")
    src=pd.read_csv(cal/"12_source_repertoire_audit.csv")
    targets=json.loads((cal/"04_observation_structure_targets.json").read_text())

    design["analysis_eligible_pair"]=parse_bool(design["analysis_eligible_pair"],"design")
    pair_params["analysis_eligible_pair"]=parse_bool(pair_params["analysis_eligible_pair"],"params")
    ad=design.loc[design["analysis_eligible_pair"]].copy()
    pp=pair_params.loc[pair_params["analysis_eligible_pair"],["subject","time","kappa"]].copy()
    ad=ad.merge(pp,on=["subject","time"],how="left",validate="one_to_one")
    require(ad["kappa"].notna().all(),"missing kappa")

    resolver=Resolver(args.source_dir)
    subjects=sorted(int(x) for x in ad.subject.unique())
    all_times=sorted(int(x) for x in design.time.unique()); tpos={t:i for i,t in enumerate(all_times)}
    gen_rows=[]; ora_rows=[]; supp_rows=[]; truth_files=[]
    totals=dict(positive_repertoire_states=0,union_timepoint_states=0,both_timepoint_states=0,
                step5_transition_rows=0,common4_rows=0)

    for subject in subjects:
        rr=refs.loc[refs.subject==subject]; require(len(rr)==1,f"reference {subject}")
        ref=pd.read_csv(cal/str(rr.iloc[0]["reference_file"]),
                        usecols=["aaSeqCDR3","f_reference_mean_repertoire"])
        ids=ref.aaSeqCDR3.astype(str).to_numpy()
        freq=ref.f_reference_mean_repertoire.to_numpy(np.float64)
        freq/=freq.sum(dtype=np.float64)
        idxmap=pd.Series(np.arange(len(ids),dtype=np.int64),index=pd.Index(ids,dtype=object))
        bio=generate_bio(np.log(freq),all_times,pars["sigma_fast"],pars["q_bio"],
                         args.seed,subject)
        sd=ad.loc[ad.subject==subject].sort_values("time")
        observed=set(int(x) for x in sd.time)
        ora_rows.extend(oracle_rows(subject,all_times,bio,observed))

        support={}
        for row in sd.itertuples(index=False):
            t=int(row.time)
            for rep in (1,2):
                ar=src.loc[(src.subject==subject)&(src.time==t)&(src.replica==rep)]
                require(len(ar)==1,f"source row {(subject,t,rep)}")
                r=ar.iloc[0]
                path=resolver.resolve(str(r.resolved_file),str(r.requested_file))
                present=read_ids(path)
                pos=idxmap.reindex(present)
                require(pos.notna().all(),f"support clone absent reference {(subject,t,rep)}")
                active=np.zeros(len(ids),bool); active[pos.to_numpy(np.int64)]=True
                target=int(getattr(row,f"n_clonotypes_rep{rep}"))
                require(int(active.sum())==target,f"richness mismatch {(subject,t,rep)}")
                support[(t,rep)]=active

        U=[]; B=[]
        for row in sd.itertuples(index=False):
            t=int(row.time); A=support[(t,1)]; C=support[(t,2)]
            U.append(A|C); B.append(A&C)
        Um=np.column_stack(U); Bm=np.column_stack(B)
        k=Um.sum(axis=1,dtype=np.int16); b=Bm.sum(axis=1,dtype=np.int16)
        repstates=sum(int(support[(int(t),r)].sum()) for t in sd.time for r in (1,2))
        ust=int(k.sum(dtype=np.int64)); bst=int(b.sum(dtype=np.int64))
        st5=int(comb2(k).sum(dtype=np.int64)); c4=int(comb2(b).sum(dtype=np.int64))
        totals["positive_repertoire_states"]+=repstates
        totals["union_timepoint_states"]+=ust
        totals["both_timepoint_states"]+=bst
        totals["step5_transition_rows"]+=st5
        totals["common4_rows"]+=c4
        supp_rows.append(dict(subject=subject,n_analysis_timepoints=len(sd),
                              positive_repertoire_states=repstates,union_timepoint_states=ust,
                              both_timepoint_states=bst,step5_rows=st5,common4_rows=c4,
                              common4_fraction=c4/st5 if st5 else np.nan))

        gt={"aaSeqCDR3":ids,"reference_frequency":freq}
        for j,t in enumerate(all_times): gt[f"logf_bio_t{t}"]=bio[j]
        gtpath=out/"ground_truth"/f"subject_{subject}.parquet"
        pd.DataFrame(gt).to_parquet(gtpath,index=False,compression="zstd")
        truth_files.append({"path":str(gtpath.relative_to(out)),"sha256":sha(gtpath),"subject":subject})

        rngc=np.random.default_rng(count_seed(args.seed,subject))
        for row in sd.itertuples(index=False):
            t=int(row.time); ptrue=np.exp(bio[tpos[t]]); kappa=float(row.kappa)
            for rep in (1,2):
                active=support[(t,rep)]; depth=int(getattr(row,f"depth_rep{rep}"))
                c=support_counts(rngc,ptrue,active,depth,kappa)
                pos=c>0; cc=c[pos]
                fn=f"{subject}_{t}-{rep}.parquet"
                pd.DataFrame({"aaSeqCDR3":ids[pos],"readCount":cc,
                              "readFraction":cc.astype(np.float64)/float(depth)}
                            ).to_parquet(out/"repertoires"/fn,index=False,compression="zstd")
                gen_rows.append(dict(subject=subject,time=t,replica=rep,file=f"repertoires/{fn}",
                                     nominal_depth=depth,realized_depth=int(cc.sum(dtype=np.int64)),
                                     target_empirical_richness=int(active.sum()),
                                     synthetic_richness=int(len(cc)),kappa=kappa,
                                     min_positive_count=int(cc.min()),max_count=int(cc.max()),
                                     n_count_1=int(np.count_nonzero(cc==1)),
                                     n_odd_positive_counts=int(np.count_nonzero(cc%2))))
        print(f"[SUBJECT {subject}] Step5={st5:,} common4={c4:,} reps={2*len(sd)}",flush=True)

    expected={
        "positive_repertoire_states":int(targets["positive_repertoire_states"]),
        "union_timepoint_states":int(targets["union_timepoint_states"]),
        "both_timepoint_states":int(targets["both_timepoint_states"]),
        "step5_transition_rows":int(targets["step5_transition_rows"]),
        "common4_rows":int(targets["common4_rows"]),
    }
    matches={k:bool(totals[k]==expected[k]) for k in expected}
    require(all(matches.values()),f"support mismatch: {totals} vs {expected}")

    gen=pd.DataFrame(gen_rows); ora=pd.DataFrame(ora_rows); supp=pd.DataFrame(supp_rows)
    gen.to_csv(out/"01_generated_sampling_summary.csv",index=False)
    ora.to_csv(out/"02_oracle_pair_summary.csv",index=False)
    supp.to_csv(out/"03_support_audit_by_subject.csv",index=False)
    require(len(gen)==int(targets["analysis_eligible_repertoires"]),"repertoire count mismatch")
    require(np.array_equal(gen.realized_depth.to_numpy(np.int64),gen.nominal_depth.to_numpy(np.int64)),"depth mismatch")
    require(np.array_equal(gen.synthetic_richness.to_numpy(np.int64),gen.target_empirical_richness.to_numpy(np.int64)),"richness mismatch")

    manifest=dict(
        schema_version=SCHEMA,script_version=VERSION,script_sha256=sha(Path(__file__).resolve()),
        python_version=platform.python_version(),numpy_version=np.__version__,pandas_version=pd.__version__,
        calibration_dir=str(cal),calibration_manifest_sha256=sha(cal/"00_manifest.json"),
        scenario=args.scenario,seed=args.seed,scenario_parameters=pars,
        validation_scope="semi-synthetic conditional on exact corrected empirical clone-specific support",
        support=dict(source="exact corrected empirical A/B positivity masks",
                     clone_specific_support_preserved=True,support_randomized_across_clones=False,
                     temporal_support_order_preserved=True,analysis_mask_preserved=True,
                     frozen_support_audit=matches,frozen_support_totals=totals),
        counts=dict(count_scale=1,exact_depth=True,exact_support=True,
                    model="one mandatory read per empirical-positive clone + Gamma-multinomial extras weighted by true biological frequency",
                    A_B_gamma_draws_independent=True),
        biology=dict(reference="f_reference_mean_repertoire",sigma_fast=pars["sigma_fast"],q_bio=pars["q_bio"],
                     ground_truth_written=True,
                     common_random_numbers=True,
                     fast_shocks_fixed_across_scenarios=True,
                     rw_unit_path_fixed_across_scenarios=True),
        anti_leakage=dict(real_cross_cov_used_to_generate=False,real_temporal_slope_used_to_generate=False,
                          real_counts_reused=False,real_positive_support_reused=True),
        truth_files=truth_files,
    )
    write_json(out/"00_manifest.json",manifest)

    print("\nEmpirical-support-conditioned synthetic experiment generated.")
    print(f"  scenario: {args.scenario}; sigma_fast={pars['sigma_fast']}; q_bio={pars['q_bio']}")
    print(f"  generated repertoires: {len(gen)}")
    print(f"  support audit: {matches}")
    print(f"  Step5/common4: {totals['step5_transition_rows']:,} / {totals['common4_rows']:,} "
          f"({totals['common4_rows']/totals['step5_transition_rows']:.6%})")
    print(f"  odd positive count entries: {int(gen.n_odd_positive_counts.sum()):,}")
    print(f"  wrote: {out}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
