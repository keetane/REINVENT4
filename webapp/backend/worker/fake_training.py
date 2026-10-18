from __future__ import annotations

import argparse
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--epochs", type=int, required=True)
    parser.add_argument("--delay", type=float, default=0.02)
    args = parser.parse_args()
    for epoch in range(1, args.epochs + 1):
        print(f"Epoch {epoch}/{args.epochs} loss={1 / epoch:.4f}", flush=True)
        time.sleep(args.delay)
    Path(args.output).write_bytes(b"fake trained model")


if __name__ == "__main__":
    main()
