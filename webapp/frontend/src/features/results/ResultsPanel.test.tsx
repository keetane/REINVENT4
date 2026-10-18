import { cleanup, render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { afterEach, describe, expect, it, vi } from "vitest";

import type { Job } from "../../api/client";
import i18n from "../../i18n";
import { ResultsPanel } from "./ResultsPanel";

function response(body: unknown) {
  return { ok: true, status: 200, json: async () => body };
}

const comparisonJob: Job = {
  id: "comparison-job-1",
  job_type: "qsar_comparison",
  status: "completed",
  comment: "P2X3 model comparison",
  created_at: "2026-10-19T00:00:00Z",
  queued_at: "2026-10-19T00:00:00Z",
  started_at: "2026-10-19T00:00:01Z",
  finished_at: "2026-10-19T00:01:00Z",
  updated_at: "2026-10-19T00:01:00Z",
  parent_smiles: "",
  parent_smiles_canonical: null,
  model_id: null,
  progress: { current: 3, total: 3, percent: 100, message: "completed" },
  metrics: {},
  exit_code: 0,
  error_summary: null
};

describe("QSAR comparison results", () => {
  afterEach(async () => {
    cleanup();
    await i18n.changeLanguage("ja");
    vi.restoreAllMocks();
  });

  it("reopens a completed comparison in result analysis", async () => {
    const fetchMock = vi.fn(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url === "/api/qsar-comparisons/comparison-job-1") return response({
        primary_metric: "rmse",
        recommended_algorithm: "extra_trees",
        selection_basis: "lowest training-only cross-validation rmse",
        successful_count: 2,
        failed_count: 1,
        results: [
          {
            algorithm: "extra_trees",
            status: "completed",
            duration_seconds: 4.25,
            metrics: {
              test: { rmse: 0.4321, r2: 0.7219 },
              cross_validation: { rmse: 0.4123, rmse_std: 0.0123, rmse_values: [0.4, 0.41, 0.427] },
              performance_gap: { rmse: 0.0198 },
              applicability_domain_fraction: 0.82,
              inference_molecules_per_second: 1234.5
            },
            plots: { actual_vs_predicted: [[5.1, 5.0], [6.2, 6.3]] }
          },
          {
            algorithm: "random_forest",
            status: "completed",
            duration_seconds: 3.1,
            metrics: {
              test: { rmse: 0.5, r2: 0.61 },
              cross_validation: { rmse: 0.48, rmse_std: 0.02 },
              performance_gap: { rmse: 0.02 },
              applicability_domain_fraction: 0.75,
              inference_molecules_per_second: 980
            },
            plots: { residual_distribution: [[-0.2, 3], [0, 12], [0.2, 4]] }
          },
          { algorithm: "svm", status: "failed", duration_seconds: 0.2, error: "training failed" }
        ]
      });
      throw new Error(`unexpected request: ${url}`);
    });
    vi.stubGlobal("fetch", fetchMock);
    const user = userEvent.setup();
    const onBack = vi.fn();

    render(<ResultsPanel job={comparisonJob} onBack={onBack} />);

    expect(await screen.findByRole("heading", { name: "P2X3 model comparison" })).toBeInTheDocument();
    expect(screen.getByText("extra_trees", { selector: ".result-summary-strip strong" })).toBeInTheDocument();
    expect(screen.getByText("0.4123 ± 0.0123")).toBeInTheDocument();
    expect(screen.getByText("0.7219")).toBeInTheDocument();
    expect(screen.getByText("82.0%")).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "extra_trees_actual_vs_predicted" })).toBeInTheDocument();
    expect(screen.queryByRole("columnheader", { name: "推論/秒" })).not.toBeInTheDocument();
    expect(screen.queryByRole("columnheader", { name: "時間" })).not.toBeInTheDocument();

    await user.selectOptions(screen.getByLabelText("分布を表示するモデル"), "random_forest");
    expect(screen.getByRole("img", { name: "random_forest_residual_distribution" })).toBeInTheDocument();
    expect(screen.queryByRole("img", { name: "extra_trees_actual_vs_predicted" })).not.toBeInTheDocument();

    await user.click(screen.getByRole("button", { name: "履歴へ戻る" }));
    expect(onBack).toHaveBeenCalledOnce();
  });
});
