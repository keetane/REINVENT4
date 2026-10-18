from __future__ import annotations

import csv
import io


SUPPORTED_DELIMITERS = ",;\t"


def csv_dict_reader(content: str) -> csv.DictReader:
    """Return a CSV reader with common ChEMBL/export delimiters detected."""
    text = content.lstrip("\ufeff")
    sample = text[:65_536]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=SUPPORTED_DELIMITERS)
    except csv.Error:
        dialect = csv.excel
    return csv.DictReader(io.StringIO(text, newline=""), dialect=dialect)

