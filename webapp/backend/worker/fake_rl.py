from __future__ import annotations

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    arguments = parser.parse_args()
    arguments.checkpoint.write_bytes(b"fake RL checkpoint")
    arguments.summary.write_text(
        "Agent,Prior,Target,Score,SMILES,SMILES_state,QSAR,QSAR (raw),step\n"
        "1.0,1.1,0.9,0.8,CCO,valid,0.8,6.2,1\n",
        encoding="utf-8",
    )
    print("fake RL completed", flush=True)


if __name__ == "__main__":
    main()
