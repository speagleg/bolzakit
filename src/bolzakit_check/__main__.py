"""python -m bolzakit_check FILE... [--deep]: prints 'FILE: VALID...' / 'FILE: INVALID: <reason>'.

Exit code: 0 if every file is unconditionally VALID; 2 if none is INVALID but at least one is only conditionally
VALID (a mode-B count, conditional on lambda_1 >= a; an lg_lower bound, conditional on its gap input); 1 if any file
is INVALID."""
import argparse
import sys

from .check import check_file


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m bolzakit_check", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--deep", action="store_true", help="also re-integrate sampled elements (consistency only)")
    args = ap.parse_args(argv)
    invalid = conditional = False
    for f in args.files:
        ok, msg = check_file(f, deep=args.deep)
        print(f"{f}: {msg}")
        invalid = invalid or not ok
        conditional = conditional or (ok and "(conditional on" in msg)
    return 1 if invalid else 2 if conditional else 0


if __name__ == "__main__":
    sys.exit(main())
