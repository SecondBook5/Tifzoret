#!/usr/bin/env python3
"""Approximate mediation-component power requirements from observed scores."""
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path
from statistics import NormalDist

def main():
    p=argparse.ArgumentParser(); p.add_argument("--inputs",required=True); p.add_argument("--outdir",required=True); a=p.parse_args(); rows=list(csv.DictReader(Path(a.inputs).open(),delimiter="\t")); out=Path(a.outdir); (out/"tables").mkdir(parents=True,exist_ok=True)
    result=[]
    for outcome in ("mediator","outcome"):
        groups={0:[],1:[]}
        for row in rows:
            try: groups[int(float(row["treatment"]))].append(float(row[outcome]))
            except (ValueError,KeyError): pass
        if all(groups.values()):
            means=[sum(values)/len(values) for values in groups.values()]; pooled=math.sqrt(sum(sum((x-sum(v)/len(v))**2 for x in v) for v in groups.values())/max(sum(map(len,groups.values()))-2,1)); effect=abs(means[1]-means[0])/pooled if pooled else 0; needed=math.ceil(2*((NormalDist().inv_cdf(.975)+NormalDist().inv_cdf(.8))/effect)**2) if effect else None
            result.append({"component":outcome,"standardized_effect":effect,"estimated_n_per_group_80_percent_power":needed})
    fields=["component","standardized_effect","estimated_n_per_group_80_percent_power"]
    with (out/"tables/mediation_power.tsv").open("w",newline="") as h: w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(result)
    (out/"mediation_power_summary.json").write_text(json.dumps({"schema_version":1,"method":"normal-approximation planning estimate based on observed standardized effects","results":result,"warnings":["Observed-effect power estimates are exploratory and can be optimistic in small samples."]},indent=2)+"\n")
if __name__=="__main__":main()
