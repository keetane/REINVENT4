from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--delay", type=float, default=0.05)
    parser.add_argument("--rows", type=int, default=2)
    parser.add_argument("--exit-code", type=int, default=0)
    parser.add_argument(
        "--output-mode", choices=["valid", "missing", "empty", "missing-smiles"], default="valid"
    )
    args = parser.parse_args()

    print("fake worker started", flush=True)
    time.sleep(max(0.0, min(args.delay, 60.0)))
    if args.exit_code:
        print("fake worker failed", flush=True)
        raise SystemExit(args.exit_code)
    if args.output_mode == "missing":
        print("fake worker omitted output", flush=True)
        return
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.output_mode == "empty":
        output.write_text("", encoding="utf-8")
        print("fake worker wrote empty output", flush=True)
        return
    with output.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["NLL"] if args.output_mode == "missing-smiles" else ["SMILES", "NLL"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(max(1, min(args.rows, 1_000))):
            row = {"NLL": float(index)}
            if "SMILES" in fieldnames:
                row["SMILES"] = "C" * (index + 1)
            writer.writerow(row)
    print("fake worker finished", flush=True)


if __name__ == "__main__":
    main()
