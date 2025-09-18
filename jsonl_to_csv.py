#!/usr/bin/env python3
import json, csv, argparse
from typing import Any, Dict, List

EVIDENCE_KEYS = [
    "GO_MF_json","GO_BP_json","GO_CC_json",
    "Catalytic_Activity_json","Binding_site_json","Cofactor_json","Active_site_json",
    "DNA_binding_json","Pathway_json","Subcellular_location_json","DomainFT_json",
    "Motif_json","Topological_domain_json","EC number","UniPathway_raw","Reactome_raw"
]

def load_jsonl(path: str):
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)

def count_items(v: Any) -> int:
    if v is None: return 0
    if isinstance(v, list): return len(v)
    if isinstance(v, dict): return len(v)
    if isinstance(v, str): return 1 if v.strip() else 0
    return 0

def first_n_terms(go_arr, n=3) -> List[str]:
    out = []
    if isinstance(go_arr, list):
        for g in go_arr:
            if isinstance(g, dict):
                term = g.get("term") or g.get("go_id")
                if isinstance(term, str) and term not in out:
                    out.append(term)
                if len(out) >= n:
                    break
    return out

def first_pathway(ph):
    if isinstance(ph, list) and ph:
        seg = ph[0]
        if isinstance(seg, dict) and isinstance(seg.get("levels"), list) and seg["levels"]:
            return " → ".join(seg["levels"])
    return None

def first_subcell(sub):
    out = []
    if isinstance(sub, list):
        for it in sub:
            if isinstance(it, dict) and isinstance(it.get("location"), str):
                out.append(it["location"])
            if len(out) >= 2:
                break
    return ", ".join(out) if out else None

def first_catalysis(cat):
    if isinstance(cat, list):
        for it in cat:
            if isinstance(it, dict):
                rxn  = it.get("reaction")
                rhea = it.get("rhea_id")
                if isinstance(rxn, str) and rxn.strip():
                    return rxn.strip(), (rhea if isinstance(rhea, str) else None)
    return None, None

def render_question(obj: Dict[str, Any]) -> str:
    prov = obj.get("provenance", {})
    organism = prov.get("organism") or "unknown organism"
    length   = prov.get("length")
    pref = f"Protein from {organism}"
    if isinstance(length, int):
        pref += f" ({length} aa)"
    pref += ". Evidence: "

    ev = obj["question"]["inputs"].get("evidence", {})
    parts = []

    ecs = ev.get("EC number")
    if isinstance(ecs, list) and ecs:
        parts.append(f"EC={'; '.join(ecs)}")

    rxn, rhea = first_catalysis(ev.get("Catalytic_Activity_json"))
    if rxn:
        s = f"Catalysis: '{rxn}'"
        if rhea: s += f" ({rhea})"
        parts.append(s)

    mf = first_n_terms(ev.get("GO_MF_json"))
    if mf:
        parts.append("MF: " + ", ".join(mf))

    pw = first_pathway(ev.get("Pathway_json"))
    if pw:
        parts.append("Pathway: " + pw)

    sub = first_subcell(ev.get("Subcellular_location_json"))
    if sub:
        parts.append("Location: " + sub)

    if not parts:
        parts = ["(no additional evidence provided)"]

    return pref + "; ".join(parts) + " What is its function?"

def build_coverage(obj: Dict[str, Any]) -> Dict[str, Any]:
    ev = obj["question"]["inputs"].get("evidence", {})
    evidence_counts = {k: count_items(ev.get(k)) for k in EVIDENCE_KEYS if k in ev}
    evidence_counts = {k:v for k,v in evidence_counts.items() if v > 0}  # compact

    gs = obj["answer"].get("gold_supports", {})
    supports_counts = {
        "EC_number":        len(gs.get("EC_number", []) or []),
        "RHEA_ids":         len(gs.get("RHEA_ids", []) or []),
        "GO_MF_ids":        len(gs.get("GO_MF_ids", []) or []),
        "GO_BP_ids":        len(gs.get("GO_BP_ids", []) or []),
        "GO_CC_ids":        len(gs.get("GO_CC_ids", []) or []),
        "Pathway_levels":   len(gs.get("Pathway_levels", []) or []),
        "Cofactor_chebi":   len(gs.get("Cofactor_chebi", []) or []),
        "Catalytic_chebi":  len(gs.get("Catalytic_chebi", []) or []),
        "Evidence_codes":   len(gs.get("Evidence_codes", []) or [])
    }
    supports_counts = {k:v for k,v in supports_counts.items() if v > 0}

    return {
        "evidence_counts": evidence_counts,
        "supports_counts": supports_counts,
        "totals": {
            "evidence_items": sum(evidence_counts.values()),
            "support_items":  sum(supports_counts.values())
        }
    }

def convert(jsonl_path: str, csv_path: str):
    with open(csv_path, "w", encoding="utf-8", newline="") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow([
            "id","entry","split","sequence","question",
            "answer","gold_supports_json","coverage_json"
        ])
        n = 0
        for obj in load_jsonl(jsonl_path):
            prov = obj.get("provenance", {})
            entry = prov.get("entry")
            split = obj.get("split", {}).get("name") or "train"
            seq   = prov.get("sequence") or ""
            question = render_question(obj)
            answer   = obj.get("answer", {}).get("function_text") or ""
            gold_supports_json = json.dumps(obj.get("answer", {}).get("gold_supports", {}), ensure_ascii=False)
            coverage_json      = json.dumps(build_coverage(obj), ensure_ascii=False)
            writer.writerow([obj.get("id"), entry, split, seq, question, answer, gold_supports_json, coverage_json])
            n += 1
    return n

def main():
    ap = argparse.ArgumentParser(description="Convert EviFuncQA JSONL to tight CSV")
    ap.add_argument("--in",  dest="inp",  required=True, help="Path to evifuncqa.jsonl")
    ap.add_argument("--out", dest="outp", required=True, help="Output CSV path")
    args = ap.parse_args()
    n = convert(args.inp, args.outp)
    print(f"Wrote {n} rows to {args.outp}")

if __name__ == "__main__":
    main()
