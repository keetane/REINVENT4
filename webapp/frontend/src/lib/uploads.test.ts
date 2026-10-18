import { describe, expect, it } from "vitest";

import {
  isDataUploadTooLarge,
  MAX_DATA_UPLOAD_BYTES,
  suggestQsarTargetColumn,
  suggestSmilesColumn
} from "./uploads";

describe("data uploads", () => {
  it("enforces a 50 MiB upload limit", () => {
    expect(MAX_DATA_UPLOAD_BYTES).toBe(52_428_800);
    expect(isDataUploadTooLarge({ size: MAX_DATA_UPLOAD_BYTES } as File)).toBe(false);
    expect(isDataUploadTooLarge({ size: MAX_DATA_UPLOAD_BYTES + 1 } as File)).toBe(true);
  });

  it("prefers the standard ChEMBL columns", () => {
    const columns = ["Molecule ChEMBL ID", "Smiles", "pChEMBL Value", "Standard Value"];
    const smiles = suggestSmilesColumn(columns);

    expect(smiles).toBe("Smiles");
    expect(suggestQsarTargetColumn(columns, smiles)).toBe("pChEMBL Value");
  });
});
