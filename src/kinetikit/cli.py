from __future__ import annotations

import argparse
from pathlib import Path

from .fitting import analyze_enzyme_csv, analyze_initial_rates_csv, analyze_rate_csv
from .plotting import save_enzyme_plot, save_json, save_rate_plot


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="kinetikit", description="Fit enzyme kinetics and general chemistry rate-law CSV data.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    enzyme = subparsers.add_parser("enzyme", help="Fit initial-rate enzyme kinetics data")
    enzyme.add_argument("csv", help="CSV with substrate and rate columns")
    enzyme.add_argument("--model", default="michaelis-menten", choices=["michaelis-menten", "substrate-inhibition"], help="enzyme model to fit")
    enzyme.add_argument("--enzyme-conc", type=float, default=None, help="enzyme concentration in the same concentration units as rate/[E] calculations")
    enzyme.add_argument("--plot", default=None, help="optional PNG output path")
    enzyme.add_argument("--json", default=None, help="optional JSON output path")

    rate = subparsers.add_parser("rate", help="Fit zero/first/second-order chemistry rate-law data")
    rate.add_argument("csv", help="CSV with time and concentration columns")
    rate.add_argument("--order", default="auto", choices=["auto", "0", "1", "2"], help="rate-law order or auto")
    rate.add_argument("--plot", default=None, help="optional PNG output path")
    rate.add_argument("--json", default=None, help="optional JSON output path")

    initial = subparsers.add_parser("initial-rates", help="Fit method-of-initial-rates power-law orders")
    initial.add_argument("csv", help="CSV with reactant_* columns and an initial_rate column")
    initial.add_argument("--reactants", nargs="*", default=None, help="optional explicit reactant columns; defaults to reactant_* columns")
    initial.add_argument("--json", default=None, help="optional JSON output path")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "enzyme":
        result = analyze_enzyme_csv(args.csv, model=args.model, enzyme_conc=args.enzyme_conc)
        if args.plot:
            save_enzyme_plot(result, args.plot)
        if args.json:
            save_json(result, args.json)
    elif args.command == "rate":
        order = args.order if args.order == "auto" else int(args.order)
        result = analyze_rate_csv(args.csv, order=order)
        if args.plot:
            save_rate_plot(result, args.plot)
        if args.json:
            save_json(result, args.json)
    elif args.command == "initial-rates":
        result = analyze_initial_rates_csv(args.csv, reactant_columns=args.reactants)
        if args.json:
            save_json(result, args.json)
    else:
        raise AssertionError("unreachable command")
    print(result.summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
