"""Command line: fetch the public data, then score molecules against it.

    toxpred fetch
    toxpred score --smiles "CC(=O)Oc1ccccc1C(=O)O"
    toxpred score --input molecules.smi --cutoff 0.6 --min-members 3
    toxpred score --smiles "..." --reach
    toxpred endpoints
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import data as D
from .consortium import Endpoint, TrafficLight

DEFAULT_ENDPOINTS = [
    "CYP450_CYP3A4", "CYP450_CYP2D6", "CYP450_CYP2C9", "CYP450_CYP2C19",
    "CYP450_CYP1A2", "Mutagenicity_Ames Mutagenicity",
    "Hepatotoxicity_Hepatotoxicity", "Cardiotoxicity_Cardiotoxicity-10",
    "Clinical Toxicity_Clinical toxicity",
]


def _read_smiles(args) -> list[str]:
    if args.smiles:
        return [args.smiles]
    if args.input:
        out = []
        for line in Path(args.input).read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line.split()[0])
        return out
    sys.exit("give --smiles or --input")


def _toxric_files(home: Path):
    root = D.fetch_toxric(home)
    return {p.stem: p for p in root.glob("**/*.csv")
            if not p.name.startswith("._")}


def cmd_fetch(args):
    home = D.cache_dir(args.home)
    print("cache: %s" % home)
    D.fetch_catmos(home)
    D.fetch_toxric(home)
    if args.with_reach:
        D.fetch_reverse_screen(home)
    print("done")


def cmd_endpoints(args):
    home = D.cache_dir(args.home)
    for name in sorted(_toxric_files(home)):
        print("  %s" % name)


def cmd_score(args):
    home = D.cache_dir(args.home)
    smiles = _read_smiles(args)
    print("cache: %s" % home)

    light = TrafficLight.from_catmos(D.fetch_catmos(home))
    print("acute toxicity reference set: %d compounds\n" % len(light))
    calls = light.predict(smiles, cutoff=args.cutoff,
                          min_members=args.min_members,
                          keep_neighbors=args.show_neighbors)
    print("ACUTE ORAL TOXICITY, red at or under 50 mg/kg, green over 2000")
    for p in calls:
        print("  %s" % p)
        for smi, sim, val in p.neighbors:
            print("      %-52s %.3f  %s" % (smi[:52], sim, val))

    files = _toxric_files(home)
    wanted = args.endpoints or DEFAULT_ENDPOINTS
    rows = []
    for name in wanted:
        if name not in files:
            continue
        ep = Endpoint.from_toxric(files[name])
        for q, p in zip(smiles, ep.predict(smiles, cutoff=args.cutoff,
                                           min_members=args.min_members)):
            rows.append((name.split("_", 1)[-1], q, p))
    if not rows:
        return
    print("\nFURTHER ENDPOINTS")
    print("  %-26s %-40s %-8s %9s %8s"
          % ("endpoint", "query", "call", "agreement", "voters"))
    for name, q, p in rows:
        if not p.answered:
            print("  %-26s %-40s %-8s %9s %8s"
                  % (name, q[:40], "no call", "", p.members))
        else:
            print("  %-26s %-40s %-8s %8.0f%% %8d"
                  % (name, q[:40], p.call, 100 * p.confidence, p.members))
    if args.reach:
        _report_reach(home, smiles, args)

    print("\nA call is issued only where the reference set holds a close enough\n"
          "neighbor. 'no call' means the method declines, which is the honest\n"
          "answer for a compound unlike anything measured.")


def _report_reach(home, smiles, args):
    from .reach import PharmCastMissing, fetch_model, read_family_map, reach
    try:
        files = D.fetch_reverse_screen(home)
        model = fetch_model(home)
        fam = read_family_map(args.family_map) if args.family_map else None
        rows = reach(smiles, files["index"], files["sites"], model,
                     family_map=fam, cutoff=args.reach_cutoff)
    except PharmCastMissing as e:
        print("\nCROSS-FAMILY REACH: skipped\n%s" % e)
        return
    print("\nCROSS-FAMILY REACH, index version %s" % files["version"])
    head = "  %-44s %9s %9s %9s" % ("query", "matches", "targets", "best sim")
    if fam:
        head += " %9s" % "families"
    print(head)
    for r in rows:
        line = "  %-44s %9d %9d %9.3f" % (r["smiles"][:44], r["matches"],
                                          r["targets_reached"],
                                          r["best_similarity"])
        if fam:
            line += " %9d" % r["families_reached"]
        print(line)
    print("  Compounds reaching widely are more often toxic. Index and "
          "fingerprint:\n  https://reversescreen.ai and https://pharmcast.ai")


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="toxpred",
        description="Toxicity prediction from structural neighbors")
    ap.add_argument("--home", default=None,
                    help="cache directory, default ~/.toxpred")
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="download the public reference data")
    f.add_argument("--with-reach", action="store_true",
                   help="also fetch the Reverse Screen index")
    f.set_defaults(func=cmd_fetch)

    e = sub.add_parser("endpoints", help="list the endpoints available")
    e.set_defaults(func=cmd_endpoints)

    s = sub.add_parser("score", help="score one or more molecules")
    s.add_argument("--smiles")
    s.add_argument("--input")
    s.add_argument("--cutoff", type=float, default=0.5,
                   help="similarity a reference compound must reach to vote")
    s.add_argument("--min-members", type=int, default=1,
                   help="voters required before a call is issued")
    s.add_argument("--endpoints", nargs="*", default=None)
    s.add_argument("--show-neighbors", type=int, default=0,
                   help="print this many of the voting compounds")
    s.add_argument("--reach", action="store_true",
                   help="also report cross-family reach; needs PharmCast")
    s.add_argument("--reach-cutoff", type=float, default=0.5)
    s.add_argument("--family-map", default=None,
                   help="optional accession,family file to group targets")
    s.set_defaults(func=cmd_score)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    main()
