import { cleanup, render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { afterEach, describe, expect, it, vi } from "vitest";

import { App } from "./App";

const emptyJobs = { items: [], total: 0, limit: 100, offset: 0 };
const resultDescriptorNames = [
  "MW", "RotB", "HBD", "HBA", "LogP", "TPSA", "QED", "SAScore",
  "FractionCSP3", "RingCount", "AromaticRingCount", "AliphaticRingCount",
  "HeavyAtomCount", "FormalCharge", "LipinskiViolations", "MolMR",
  "ParentSimilarity", "ScaffoldSimilarity"
];

function response(body: unknown, ok = true, status = 200) {
  return { ok, status, json: async () => body };
}

function installBaseFetch() {
  const fetchMock = vi.fn(async (input: RequestInfo | URL) => {
    const url = String(input);
    if (url === "/api/health") {
      return response({ status: "ok", database: "ok", version: "0.1.0" });
    }
    if (url.startsWith("/api/jobs?")) return response(emptyJobs);
    throw new Error(`unexpected request: ${url}`);
  });
  vi.stubGlobal("fetch", fetchMock);
  return fetchMock;
}

describe("App", () => {
  afterEach(() => {
    cleanup();
    window.localStorage.clear();
    vi.restoreAllMocks();
  });

  it("shows the API version and the job creation workspace", async () => {
    installBaseFetch();

    render(<App />);

    expect(await screen.findByText("API v0.1.0")).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "Molecular Generator" })).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "Parent 分子から探索を開始" })).toBeInTheDocument();
    expect(screen.getByLabelText("Molecular Name")).toBeInTheDocument();
    expect(screen.getByLabelText("Parent SMILES")).toHaveAttribute("rows", "2");
    const device = screen.getByLabelText("実行デバイス") as HTMLSelectElement;
    expect(device).toHaveValue("mps");
    expect(Array.from(device.options).map((option) => option.text)).toEqual(["CPU", "GPU", "MPS"]);
  });

  it("fills Parent SMILES from a PubChem molecular-name lookup", async () => {
    const fetchMock = installBaseFetch();
    fetchMock.mockImplementation(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url === "/api/health") return response({ status: "ok", database: "ok", version: "0.1.0" });
      if (url.startsWith("/api/jobs?")) return response(emptyJobs);
      if (url === "/api/parents/pubchem?name=Aspirin") {
        return response({
          cid: 2244,
          molecular_name: "Aspirin",
          isomeric_smiles: "CC(=O)OC1=CC=CC=C1C(=O)O"
        });
      }
      throw new Error(`unexpected request: ${url}`);
    });
    const user = userEvent.setup();
    render(<App />);

    await user.type(screen.getByLabelText("Molecular Name"), "Aspirin");
    await user.click(screen.getByRole("button", { name: "PubChemから取得" }));

    expect(screen.getByLabelText("Parent SMILES")).toHaveValue("CC(=O)OC1=CC=CC=C1C(=O)O");
    expect(await screen.findByText(/PubChem CID 2244/)).toBeInTheDocument();
  });

  it("shows tutorial samples and applies them to the creation form", async () => {
    installBaseFetch();
    const user = userEvent.setup();
    render(<App />);

    await user.click(screen.getByRole("button", { name: /00 Tutorial/ }));
    expect(screen.getByRole("heading", { name: "分子生成をはじめる" })).toBeInTheDocument();
    expect(screen.getByText("N#CC[C@H](C1CCCC1)[*:1]", { exact: true })).toBeInTheDocument();
    expect(screen.getByText("Oc1cncc(*)c1|*c1ccoc1", { exact: true })).toBeInTheDocument();
    expect(screen.getByText("RDKit")).toBeInTheDocument();

    await user.click(screen.getAllByRole("button", { name: "この例を使う" })[1]);
    expect(screen.getByLabelText("Model")).toHaveValue("linkinvent");
    expect(screen.getByLabelText("LinkInvent warheads")).toHaveValue("Oc1cncc(*)c1|*c1ccoc1");
  });

  it("shows a correctable backend validation error", async () => {
    const fetchMock = installBaseFetch();
    fetchMock.mockImplementation(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url === "/api/health") return response({ status: "ok", database: "ok", version: "0.1.0" });
      if (url.startsWith("/api/jobs?")) return response(emptyJobs);
      if (url === "/api/parents/validate") {
        return response(
          { error: { code: "invalid_parent_smiles", message: "SMILES を解釈できません。" } },
          false,
          422
        );
      }
      throw new Error(`unexpected request: ${url}`);
    });
    const user = userEvent.setup();
    render(<App />);

    await user.type(screen.getByLabelText("Parent SMILES"), "not smiles");
    await user.click(screen.getByRole("button", { name: "検証する" }));

    expect(await screen.findByRole("alert")).toHaveTextContent("SMILES を解釈できません。");
  });

  it("uses model-specific LibInvent and LinkInvent inputs", async () => {
    let validationPayload: Record<string, unknown> | null = null;
    vi.stubGlobal(
      "fetch",
      vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
        const url = String(input);
        if (url === "/api/health") {
          return response({ status: "ok", database: "ok", version: "0.1.0" });
        }
        if (url.startsWith("/api/jobs?")) return response(emptyJobs);
        if (url === "/api/parents/validate") {
          validationPayload = JSON.parse(String(init?.body));
          return response({
            input_smiles: "Oc1cncc(*)c1|*c1ccoc1",
            canonical_smiles: "*c1ccoc1|*c1cncc(O)c1",
            molblock: "mol",
            atom_count: 14,
            molecular_formula: "C9H7*2NO2"
          });
        }
        throw new Error(`unexpected request: ${url}`);
      })
    );
    const user = userEvent.setup();
    render(<App />);

    await user.selectOptions(screen.getByLabelText("Model"), "libinvent");
    expect(screen.getByLabelText("LibInvent scaffold")).toBeInTheDocument();
    expect(screen.getByText(/側鎖を発生させる位置に \*/)).toBeInTheDocument();

    await user.selectOptions(screen.getByLabelText("Model"), "linkinvent");
    const input = screen.getByLabelText("LinkInvent warheads");
    await user.type(input, "Oc1cncc(*)c1|*c1ccoc1");
    await user.click(screen.getByRole("button", { name: "検証する" }));

    expect(await screen.findByText("有効な入力")).toBeInTheDocument();
    expect(validationPayload).toEqual({
      smiles: "Oc1cncc(*)c1|*c1ccoc1",
      model_key: "linkinvent"
    });
  });

  it("creates a job and switches to its persistent history", async () => {
    const job = {
      id: "12345678-1234-5678-1234-567812345678",
      status: "queued",
      created_at: "2026-10-16T12:00:00Z",
      queued_at: "2026-10-16T12:00:00Z",
      started_at: null,
      finished_at: null,
      updated_at: "2026-10-16T12:00:00Z",
      comment: "first run",
      parent_smiles: "CCO",
      parent_smiles_canonical: "CCO",
      progress: { current: null, total: null, percent: null, message: "queued" },
      exit_code: null,
      error_summary: null
    };
    let created = false;
    let createPayload: Record<string, unknown> | null = null;
    vi.stubGlobal(
      "fetch",
      vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
        const url = String(input);
        if (url === "/api/health") return response({ status: "ok", database: "ok", version: "0.1.0" });
        if (url === "/api/parents/validate") {
          return response({
            input_smiles: "CCO",
            canonical_smiles: "CCO",
            molblock: "mol",
            atom_count: 3,
            molecular_formula: "C2H6O"
          });
        }
        if (url === "/api/jobs") {
          created = true;
          createPayload = JSON.parse(String(init?.body));
          return response(job);
        }
        if (url.startsWith("/api/jobs?")) {
          return response(created ? { items: [job], total: 1, limit: 100, offset: 0 } : emptyJobs);
        }
        if (url === `/api/jobs/${job.id}`) return response(job);
        if (url.startsWith(`/api/jobs/${job.id}/logs`)) {
          return response({ job_id: job.id, lines: [], truncated: false });
        }
        throw new Error(`unexpected request: ${url}`);
      })
    );
    const user = userEvent.setup();
    render(<App />);

    await user.type(screen.getByLabelText("Parent SMILES"), "CCO");
    await user.type(screen.getByLabelText("コメント"), "first run");
    await user.click(screen.getByRole("button", { name: "ジョブを開始" }));

    expect(await screen.findByRole("heading", { name: "ジョブ履歴" })).toBeInTheDocument();
    expect(screen.getAllByText("first run").length).toBeGreaterThan(0);
    expect(window.localStorage.getItem("reinvent4.selectedJobId")).toBe(job.id);
    expect(createPayload).toMatchObject({ run_conditions: { device: "mps" } });
  });

  it("opens completed results with filters and matched downloads", async () => {
    const job = {
      id: "12345678-1234-5678-1234-567812345678",
      status: "completed",
      created_at: "2026-10-16T12:00:00Z",
      queued_at: "2026-10-16T12:00:00Z",
      started_at: "2026-10-16T12:00:01Z",
      finished_at: "2026-10-16T12:00:02Z",
      updated_at: "2026-10-16T12:00:02Z",
      comment: "result run",
      parent_smiles: "CCO",
      parent_smiles_canonical: "CCO",
      progress: { current: 1, total: 1, percent: 100, message: "completed" },
      exit_code: 0,
      error_summary: null
    };
    const descriptors = resultDescriptorNames.map((name) => ({
      name, source: "RDKit", aliases: [], minimum: 0, maximum: 100, missing_count: 0
    }));
    const values = Object.fromEntries(resultDescriptorNames.map((name) => [name, 1]));
    vi.stubGlobal("fetch", vi.fn(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url === "/api/health") return response({ status: "ok", database: "ok", version: "0.1.0" });
      if (url.startsWith("/api/jobs?")) return response({ items: [job], total: 1, limit: 100, offset: 0 });
      if (url === `/api/jobs/${job.id}`) return response(job);
      if (url.startsWith(`/api/jobs/${job.id}/logs`)) return response({ job_id: job.id, lines: [], truncated: false });
      if (url.endsWith("/results/schema")) return response({ job_id: job.id, columns: ["SMILES"], descriptors, total_rows: 1, valid_rows: 1, invalid_smiles: 0, warnings: [], rdkit_version: "2024.09.6" });
      if (url.includes("/results/molecules")) return response({ job_id: job.id, items: [{ row_number: 1, smiles: "CCO", canonical_smiles: "CCO", structure_svg: "<svg xmlns='http://www.w3.org/2000/svg'/>", values: { ...values, MW: 46, LogP: -0.1, TPSA: 20 } }], total: 1, limit: 50, offset: 0, invalid_smiles: 0 });
      if (url.includes("/results/statistics")) return response({ job_id: job.id, total_rows: 1, filtered_rows: 1, invalid_smiles: 0, summaries: Object.fromEntries(descriptors.map(({ name }) => [name, { count: 1, missing_count: 0, minimum: 1, maximum: 1, mean: 1, median: 1, standard_deviation: 0, q1: 1, q3: 1 }])), histograms: Object.fromEntries(descriptors.map(({ name }) => [name, [{ start: 1, end: 1, count: 1 }]])), scatter_x: "MW", scatter_y: "LogP", scatter_sampled: false, scatter: [{ row_number: 1, x: 46, y: -0.1 }] });
      throw new Error(`unexpected request: ${url}`);
    }));
    const user = userEvent.setup();
    render(<App />);

    await user.click(screen.getByRole("button", { name: /進捗と履歴/ }));
    await user.click(await screen.findByText("result run"));
    await user.click(screen.getByRole("button", { name: "結果を開く" }));

    expect(await screen.findByRole("heading", { name: "result run" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "表" })).toHaveClass("active");
    expect(screen.getByRole("columnheader", { name: "Structure" })).toBeInTheDocument();
    expect(screen.queryByRole("columnheader", { name: "SMILES" })).not.toBeInTheDocument();
    const logpSort = screen.getByRole("button", { name: "LogPで並び替え" });
    await user.click(logpSort);
    expect(logpSort.closest("th")).toHaveAttribute("aria-sort", "ascending");
    await user.click(logpSort);
    expect(logpSort.closest("th")).toHaveAttribute("aria-sort", "descending");
    expect(screen.getByLabelText("MW 最小")).toBeInTheDocument();
    expect(screen.getByLabelText("MW 下限スライダー")).toBeInTheDocument();
    expect(screen.getByText("追加物性・類似度").closest("details")).not.toHaveAttribute("open");
    expect(screen.getByLabelText("QED 最小")).toBeInTheDocument();
    expect(screen.getAllByText("ParentSimilarity").length).toBeGreaterThan(0);
    expect(document.querySelectorAll(".axis-tick").length).toBeGreaterThan(0);
    expect(screen.getByLabelText("MW 横軸")).toBeInTheDocument();
    expect(screen.getByRole("link", { name: "CSV" })).toHaveAttribute("href", expect.stringContaining("download.csv"));
    expect((await screen.findAllByText("CCO")).length).toBeGreaterThan(0);
  });
});
