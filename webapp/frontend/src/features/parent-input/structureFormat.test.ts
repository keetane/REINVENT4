import { describe, expect, it } from "vitest";

import { fromEditorStructure, toEditorStructure } from "./structureFormat";

describe("Ketcher structure formatting", () => {
  it("shows LinkInvent warheads as disconnected editor fragments", () => {
    expect(toEditorStructure("*c1ccccc1|*c1ccoc1", true)).toBe(
      "*c1ccccc1.*c1ccoc1"
    );
  });

  it("restores the LinkInvent separator and removes CXSMILES annotations", () => {
    expect(
      fromEditorStructure(
        "OC1:C:C(*):C:N:C:1.*C1:C:O:C:C:1 |$;;;;;;;;;;;;;$|",
        true
      )
    ).toBe("OC1:C:C(*):C:N:C:1|*C1:C:O:C:C:1");
  });

  it("removes CXSMILES annotations for a single molecule", () => {
    expect(fromEditorStructure("CCO |$;;;$|", false)).toBe("CCO");
  });
});
