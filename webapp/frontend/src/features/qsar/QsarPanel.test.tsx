import { cleanup, render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { afterEach, describe, expect, it, vi } from "vitest";

import i18n from "../../i18n";
import { QsarPanel } from "./QsarPanel";

function response(body: unknown) {
  return { ok: true, status: 200, json: async () => body };
}

describe("QsarPanel model evaluation", () => {
  afterEach(async () => {
    cleanup();
    await i18n.changeLanguage("ja");
    vi.restoreAllMocks();
  });

  it("predicts a molecule drawn or entered for a registered model", async () => {
    const model = {
      id: "e3e31ec2-a1f1-40f2-89d6-8eb4af04f9d6",
      name: "activity_model",
      algorithm: "lightgbm",
      task_type: "regression",
      target_name: "pChEMBL Value",
      target_unit: "as_provided",
      created_at: "2026-10-18T00:00:00Z",
      training_data_count: 80,
      test_data_count: 20,
      metrics: { test: { r2: 0.5, rmse: 0.4, mae: 0.3 } },
      metadata: {},
      ad_threshold: 0.4,
      in_use: false
    };
    const fetchMock = vi.fn(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url === "/api/qsar/capabilities") return response({ algorithms: [], devices: [] });
      if (url === "/api/qsar-models") return response({ items: [model] });
      if (url === "/api/parents/validate") return response({ input_smiles: "CCO", canonical_smiles: "CCO", molblock: "", atom_count: 3, molecular_formula: "C2H6O" });
      if (url === "/api/qsar/predict") return response({ model_id: model.id, predictions: [{ smiles: "CCO", standardized_smiles: "CCO", molecular_weight: 46.07, formal_charge: 0, valid: true, prediction: 6.25, scoring_prediction: 6.25, nearest_tanimoto: 0.72, in_ad: true, error: null }] });
      throw new Error(`unexpected request: ${url}`);
    });
    vi.stubGlobal("fetch", fetchMock);
    const user = userEvent.setup();

    render(<QsarPanel onCreated={() => undefined} />);
    await user.click(await screen.findByRole("button", { name: "新規構造を予測" }));
    await user.type(screen.getByLabelText("評価対象 SMILES"), "CCO");
    await user.click(screen.getByRole("button", { name: "予測値を計算" }));

    expect(await screen.findByText("6.25000")).toBeInTheDocument();
    expect(screen.getByText("In")).toBeInTheDocument();
    expect(screen.getByText(/最近傍Tanimoto 0.720/)).toBeInTheDocument();
  });
});
