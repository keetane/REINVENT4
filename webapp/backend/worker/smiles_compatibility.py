"""Check standardized SMILES against a REINVENT checkpoint vocabulary."""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

import torch


def inspect(model_path: Path, smiles: list[str]) -> list[dict[str, object]]:
    checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
    vocabulary = checkpoint.get("vocabulary")
    if vocabulary is None:
        raise ValueError("checkpoint vocabulary is missing")
    module = importlib.import_module(type(vocabulary).__module__)
    tokenizer = module.SMILESTokenizer()
    maximum = int(checkpoint.get("max_sequence_length") or 0)
    rows: list[dict[str, object]] = []
    for value in smiles:
        try:
            tokens = tokenizer.tokenize(value)
            vocabulary.encode(tokens)
            if maximum and len(tokens) > maximum:
                raise ValueError(f"sequence_length_{len(tokens)}_exceeds_{maximum}")
            rows.append({"smiles": value, "compatible": True, "reason": None})
        except (KeyError, ValueError, RuntimeError) as error:
            rows.append({"smiles": value, "compatible": False, "reason": str(error)})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, type=Path)
    arguments = parser.parse_args()
    payload = json.load(sys.stdin)
    print(json.dumps({"items": inspect(arguments.model, payload["smiles"])}))


if __name__ == "__main__":
    main()
