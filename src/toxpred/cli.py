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
from .contribute import SUBMIT_URL as CONTRIB_URL

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
        D.fetch_reverse_screen(home, version=args.index_version)
    print("done")


# The paper's figures, in the order the paper prints them. The scripts live in
# figures/ beside the measured result files they read, so a figure can be
# redrawn without rerunning the analysis.
FIGURE_SCRIPTS = [
    ("Figures 1 and 2", "make_figures_paper.py"),
    ("Figure 3", "make_figure4_breadth.py"),
    ("Figure 4", "make_figure_structures.py"),
    ("Figure 5", "make_figure7_scaffold_blind.py"),
    ("Figure 6", "make_figure6_useful.py"),
]


def _figures_dir() -> Path:
    """figures/ as shipped in the repository, from wherever this is installed."""
    here = Path(__file__).resolve()
    for base in (here.parents[2], here.parents[3] if len(here.parents) > 3
                 else here.parents[2]):
        cand = base / "figures"
        if (cand / "make_figures_paper.py").exists():
            return cand
    raise SystemExit(
        "the figures directory is not beside this installation. It ships with "
        "the repository, so clone it and run from the clone:\n"
        "  git clone https://github.com/smuskal/toxpred.git")


def cmd_figures(args):
    """Redraw every figure in the paper from the measured result files."""
    import os
    import subprocess

    figs = _figures_dir()
    out = Path(args.out).resolve() if args.out else figs / "out"
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TOXPRED_FIGURE_OUT=str(out))
    if args.data:
        env["TOXPRED_FIGURE_DATA"] = str(Path(args.data).resolve())
    print("figures -> %s" % out)
    failed = []
    for label, script in FIGURE_SCRIPTS:
        print("\n%s: %s" % (label, script))
        r = subprocess.run([sys.executable, str(figs / script)], env=env)
        if r.returncode != 0:
            failed.append("%s (%s)" % (label, script))
    if failed:
        raise SystemExit("\nthese did not draw: %s" % ", ".join(failed))
    print("\nevery figure in the paper is in %s" % out)


def cmd_endpoints(args):
    home = D.cache_dir(args.home)
    for name in sorted(_toxric_files(home)):
        print("  %s" % name)


def cmd_provenance(args):
    import json
    home = D.cache_dir(args.home)
    man = D.read_manifest(home)
    if not man:
        print("nothing fetched yet; run `toxpred fetch`")
        return
    print("cache: %s\n" % home)
    for key in sorted(man):
        r = man[key]
        print("%s" % key)
        print("   from    %s" % r["url"])
        print("   file    %s, %s bytes" % (r["file"], format(r["bytes"], ",")))
        print("   sha256  %s" % r["sha256"])
        print("   fetched %s" % r["fetched"])
        if r.get("note"):
            print("   note    %s" % r["note"])
        print()
    print(json.dumps({k: v["sha256"][:12] for k, v in man.items()},
                     indent=1, sort_keys=True))


def cmd_score(args):
    home = D.cache_dir(args.home)
    smiles = _read_smiles(args)
    print("cache: %s" % home)

    light = TrafficLight.from_catmos(D.fetch_catmos(home))
    print("acute toxicity reference set: %d compounds, fingerprint %s\n"
          % (len(light), light.fingerprint))
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
        files = D.fetch_reverse_screen(home, version=args.index_version)
        model = fetch_model(home)
        # No family map given means the shipped grouping, so the count is a
        # count of families rather than of proteins without anything supplied.
        fam = read_family_map(args.family_map) if args.family_map else None
        rows = reach(smiles, files["index"], files["sites"], model,
                     family_map=fam, cutoff=args.reach_cutoff,
                     top_k=args.reach_top_k)
    except PharmCastMissing as e:
        print("\nCROSS-FAMILY REACH: skipped\n%s" % e)
        return
    print("\nCROSS-FAMILY REACH, index version %s, fingerprint %s"
          % (files["version"], model.name))
    head = ("  %-44s %9s %9s %9s %9s %9s"
            % ("query", "matches", "targets", "families", "antitgts", "best sim"))
    print(head)
    for r in rows:
        line = ("  %-44s %9d %9d %9d %9d %9.3f"
                % (r["smiles"][:44], r["matches"], r["targets_reached"],
                   r.get("families_reached", 0), r["antitargets_reached"],
                   r["best_similarity"]))
        print(line)
    print("  Compounds reaching widely are more often toxic. Index and "
          "fingerprint:\n  https://reversescreen.ai and https://pharmcast.ai")


def cmd_contribute(args):
    from .contribute import build, submit
    meta = build(args.input, args.out, fmt=args.format, home=args.home,
                 cutoff=args.reach_cutoff, family_map=args.family_map)
    print("wrote %s" % meta["written_to"])
    print("  format      %s" % meta["format"])
    print("  compounds   %d" % meta["n_compounds"])
    print("  endpoints   %s" % (", ".join(meta["endpoints"]) or "none found"))
    if meta.get("index_version"):
        print("  index       %s, fingerprint %s"
              % (meta["index_version"], meta["fingerprint"]))
    print("\n%s" % meta["note"])
    if not args.submit:
        print("\nNothing was uploaded. Read the file, then send it with\n"
              "  toxpred contribute ... --submit --email you@example.com")
        return
    if not args.email:
        raise SystemExit("--submit needs --email, which is where the receipt "
                         "and any question about the record come back to")
    reply = submit(args.out, args.email, url=args.submit_url)
    print("\nsent to %s" % args.submit_url)
    print("  receipt     %s" % reply.get("receipt", "(none returned)"))
    print("  queued      %d record%s"
          % (reply.get("n_records", meta["n_compounds"]),
             "" if reply.get("n_records", meta["n_compounds"]) == 1 else "s"))
    print("\nA submitted record changes no prediction until a reviewer moves it\n"
          "into the reference set. Scoring answers from the published set until\n"
          "then, which is the set the paper reports.")


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
    f.add_argument("--index-version", default=None,
                   help="pin the Reverse Screen index release, for example "
                        "2026-09-20, the release the paper reports. Without "
                        "this the current published release is taken.")
    f.set_defaults(func=cmd_fetch)

    g = sub.add_parser("figures",
                       help="redraw every figure in the paper from the "
                            "measured result files")
    g.add_argument("--out", default=None,
                   help="where to write them, default figures/out")
    g.add_argument("--data", default=None,
                   help="a results directory from your own run, default the "
                        "measured files shipped in figures/data")
    g.set_defaults(func=cmd_figures)

    e = sub.add_parser("endpoints", help="list the endpoints available")
    e.set_defaults(func=cmd_endpoints)

    pv = sub.add_parser("provenance",
                        help="what was downloaded, from where, and its checksum")
    pv.set_defaults(func=cmd_provenance)

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
    s.add_argument("--reach-top-k", type=int, default=500,
                   help="most similar indexed ligands kept before the cutoff "
                        "is applied; 500 is the published operating point")
    s.add_argument("--index-version", default=None,
                   help="pin the Reverse Screen index release, for example "
                        "2026-09-20, the release the paper reports")
    s.add_argument("--family-map", default=None,
                   help="your own accession,family file. Without one the "
                        "grouping the paper counted, which ships here, is used")
    s.set_defaults(func=cmd_score)

    c = sub.add_parser("contribute",
                       help="turn your molecules and measurements into a "
                            "poolable record")
    c.add_argument("--input", required=True,
                   help="CSV or TSV with a smiles column; every other column "
                        "is treated as an endpoint measurement")
    c.add_argument("--out", required=True, help="file to write")
    c.add_argument("--format", default="counts", choices=["counts",
                                                          "fingerprint"],
                   help="counts releases four integers and no structure; "
                        "fingerprint carries both signals but published work "
                        "recovers a fraction of structures from it")
    c.add_argument("--reach-cutoff", type=float, default=0.5)
    c.add_argument("--family-map", default=None)
    c.add_argument("--submit", action="store_true",
                   help="send the record to the pool at toxpred.ai once it is "
                        "written, and print the receipt")
    c.add_argument("--email", default=None,
                   help="where the receipt and any question come back to; "
                        "required with --submit")
    c.add_argument("--submit-url", default=CONTRIB_URL,
                   help="where to send it, for a pool of your own")
    c.set_defaults(func=cmd_contribute)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    main()
