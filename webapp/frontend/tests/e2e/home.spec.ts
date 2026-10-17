import { expect, test } from "@playwright/test";

const descriptorNames = [
  "MW", "RotB", "HBD", "HBA", "LogP", "TPSA", "QED", "SAScore",
  "FractionCSP3", "RingCount", "AromaticRingCount", "AliphaticRingCount",
  "HeavyAtomCount", "FormalCharge", "LipinskiViolations", "MolMR",
  "ParentSimilarity", "ScaffoldSimilarity"
];

test("opens the tutorial and applies a LinkInvent sample", async ({ page }) => {
  await page.route("http://127.0.0.1:5173/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (url.pathname === "/api/health") return route.fulfill({ json: { status: "ok", database: "ok", version: "0.1.0" } });
    if (url.pathname === "/api/jobs") return route.fulfill({ json: { items: [], total: 0, limit: 100, offset: 0 } });
    return route.abort();
  });

  await page.goto("/");
  await expect(page.locator('link[rel="icon"]')).toHaveAttribute("href", "/favicon.svg");
  await page.getByRole("button", { name: /00 Tutorial/ }).click();
  await expect(page.getByRole("heading", { name: "分子生成をはじめる" })).toBeVisible();
  await expect(page.getByText("REINVENT4", { exact: true }).last()).toBeVisible();
  await page.getByRole("button", { name: "この例を使う" }).nth(1).click();
  await expect(page.getByLabel("Model")).toHaveValue("linkinvent");
  await expect(page.getByLabel("LinkInvent warheads")).toHaveValue("Oc1cncc(*)c1|*c1ccoc1");
});

test("creates a job and restores it in progress history", async ({ page }) => {
  const job = {
    id: "12345678-1234-5678-1234-567812345678",
    status: "queued",
    created_at: "2026-10-16T12:00:00Z",
    queued_at: "2026-10-16T12:00:00Z",
    started_at: null,
    finished_at: null,
    updated_at: "2026-10-16T12:00:00Z",
    comment: "E2E run",
    parent_smiles: "CCO",
    parent_smiles_canonical: "CCO",
    progress: { current: null, total: null, percent: null, message: "queued" },
    exit_code: null,
    error_summary: null
  };
  let created = false;

  await page.route("http://127.0.0.1:5173/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (url.pathname === "/api/health") {
      await route.fulfill({ json: { status: "ok", database: "ok", version: "0.1.0" } });
    } else if (url.pathname === "/api/parents/validate") {
      await route.fulfill({
        json: {
          input_smiles: "CCO",
          canonical_smiles: "CCO",
          molblock: "mol",
          atom_count: 3,
          molecular_formula: "C2H6O"
        }
      });
    } else if (url.pathname === "/api/parents/pubchem") {
      await route.fulfill({
        json: { cid: 2244, molecular_name: "Aspirin", isomeric_smiles: "CCO" }
      });
    } else if (url.pathname === "/api/jobs" && route.request().method() === "POST") {
      created = true;
      await route.fulfill({ status: 202, json: job });
    } else if (url.pathname === "/api/jobs") {
      await route.fulfill({ json: created ? { items: [job], total: 1, limit: 100, offset: 0 } : { items: [], total: 0, limit: 100, offset: 0 } });
    } else if (url.pathname === `/api/jobs/${job.id}/logs`) {
      await route.fulfill({ json: { job_id: job.id, lines: [], truncated: false } });
    } else if (url.pathname === `/api/jobs/${job.id}`) {
      await route.fulfill({ json: job });
    } else {
      await route.abort();
    }
  });

  await page.goto("/");
  await expect(page).toHaveTitle("Molecular Generator");
  await expect(page.getByRole("heading", { name: "Molecular Generator" })).toBeVisible();
  await expect(page.getByLabel("Parent SMILES")).toHaveAttribute("rows", "2");
  await page.getByLabel("Molecular Name").fill("Aspirin");
  await page.getByRole("button", { name: "PubChemから取得" }).click();
  await expect(page.getByLabel("Parent SMILES")).toHaveValue("CCO");
  await expect(page.getByText(/PubChem CID 2244/)).toBeVisible();
  await page.getByLabel("コメント").fill("E2E run");
  await page.getByRole("button", { name: "ジョブを開始" }).click();
  await expect(page.getByRole("heading", { name: "ジョブ履歴" })).toBeVisible();
  await expect(page.getByText("E2E run").first()).toBeVisible();

  await page.reload();
  await page.getByRole("button", { name: /進捗と履歴/ }).click();
  await expect(page.getByText("E2E run").first()).toBeVisible();
  await expect(page.getByText("計測なし")).toBeVisible();
});

test("loads the local Ketcher editor without an external structure service", async ({ page }) => {
  await page.route("http://127.0.0.1:5173/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (url.pathname === "/api/health") {
      await route.fulfill({ json: { status: "ok", database: "ok", version: "0.1.0" } });
    } else if (url.pathname === "/api/jobs") {
      await route.fulfill({ json: { items: [], total: 0, limit: 100, offset: 0 } });
    } else {
      await route.abort();
    }
  });

  await page.goto("/");
  await page.getByRole("button", { name: "構造エディタを開く" }).click();
  await expect(page.getByText("ready", { exact: true })).toBeVisible({ timeout: 20_000 });
  await expect(page.getByLabel("Ketcher 構造エディタ")).toBeVisible();
  await page.waitForTimeout(2_000);
  await expect(page.getByText("構造エディタを読み込めませんでした。", { exact: false })).toHaveCount(0);
});

test("opens a completed job result and applies a descriptor filter", async ({ page }) => {
  const id = "12345678-1234-5678-1234-567812345678";
  const job = { id, status: "completed", created_at: "2026-10-16T12:00:00Z", queued_at: "2026-10-16T12:00:00Z", started_at: "2026-10-16T12:00:01Z", finished_at: "2026-10-16T12:00:02Z", updated_at: "2026-10-16T12:00:02Z", comment: "analysis run", parent_smiles: "CCO", parent_smiles_canonical: "CCO", progress: { current: 2, total: 2, percent: 100, message: "completed" }, exit_code: 0, error_summary: null };
  const names = descriptorNames;
  const descriptors = names.map((name) => ({ name, source: "RDKit", aliases: [], minimum: 0, maximum: 100, missing_count: 0 }));
  const descriptorValues = Object.fromEntries(names.map((name) => [name, 1]));

  await page.route("http://127.0.0.1:5173/api/**", async (route) => {
    const url = new URL(route.request().url());
    if (url.pathname === "/api/health") return route.fulfill({ json: { status: "ok", database: "ok", version: "0.1.0" } });
    if (url.pathname === "/api/jobs") return route.fulfill({ json: { items: [job], total: 1, limit: 100, offset: 0 } });
    if (url.pathname === `/api/jobs/${id}`) return route.fulfill({ json: job });
    if (url.pathname === `/api/jobs/${id}/logs`) return route.fulfill({ json: { job_id: id, lines: [], truncated: false } });
    if (url.pathname.endsWith("/results/schema")) return route.fulfill({ json: { job_id: id, columns: ["SMILES"], descriptors, total_rows: 2, valid_rows: 2, invalid_smiles: 0, warnings: [], rdkit_version: "2024.09.6" } });
    if (url.pathname.endsWith("/results/molecules")) {
      const filtered = url.searchParams.get("mw_min") === "50";
      const values = filtered ? [{ row_number: 2, smiles: "CCCC", canonical_smiles: "CCCC", structure_svg: "<svg xmlns='http://www.w3.org/2000/svg'/>", values: { ...descriptorValues, MW: 58, LogP: 1.8, TPSA: 0 } }] : [{ row_number: 1, smiles: "CCO", canonical_smiles: "CCO", structure_svg: "<svg xmlns='http://www.w3.org/2000/svg'/>", values: { ...descriptorValues, MW: 46, LogP: -0.1, TPSA: 20 } }, { row_number: 2, smiles: "CCCC", canonical_smiles: "CCCC", structure_svg: "<svg xmlns='http://www.w3.org/2000/svg'/>", values: { ...descriptorValues, MW: 58, LogP: 1.8, TPSA: 0 } }];
      return route.fulfill({ json: { job_id: id, items: values, total: values.length, limit: 50, offset: 0, invalid_smiles: 0 } });
    }
    if (url.pathname.endsWith("/results/statistics")) return route.fulfill({ json: { job_id: id, total_rows: 2, filtered_rows: url.searchParams.get("mw_min") === "50" ? 1 : 2, invalid_smiles: 0, summaries: Object.fromEntries(names.map((name) => [name, { count: 1, missing_count: 0, minimum: 1, maximum: 1, mean: 1, median: 1, standard_deviation: 0, q1: 1, q3: 1 }])), histograms: Object.fromEntries(names.map((name) => [name, [{ start: 1, end: 1, count: 1 }]])), scatter_x: "MW", scatter_y: "LogP", scatter_sampled: false, scatter: [{ row_number: 2, x: 58, y: 1.8 }] } });
    return route.abort();
  });

  await page.goto("/");
  await page.getByRole("button", { name: /進捗と履歴/ }).click();
  await page.getByText("analysis run").first().click();
  await page.getByRole("button", { name: "結果を開く" }).click();
  await expect(page.getByRole("heading", { name: "analysis run" })).toBeVisible();
  await expect(page.getByRole("button", { name: "表" })).toHaveClass(/active/);
  await expect(page.locator(".result-table th:nth-child(2)")).toHaveCSS("position", "sticky");
  await expect(page.locator(".result-table th", { hasText: /^SMILES$/ })).toHaveCount(0);
  await page.getByRole("button", { name: "LogPで並び替え" }).click();
  await expect(page.locator(".result-table th", { has: page.getByRole("button", { name: "LogPで並び替え" }) })).toHaveAttribute("aria-sort", "ascending");
  await page.getByRole("button", { name: "LogPで並び替え" }).click();
  await expect(page.locator(".result-table th", { has: page.getByRole("button", { name: "LogPで並び替え" }) })).toHaveAttribute("aria-sort", "descending");
  await expect(page.getByLabel("MW 下限スライダー")).toBeVisible();
  await expect(page.getByLabel("QED 最小")).toBeHidden();
  await page.getByText("追加物性・類似度", { exact: true }).click();
  await expect(page.getByLabel("QED 最小")).toBeVisible();
  await page.getByLabel("MW 最小").fill("50");
  await page.getByRole("button", { name: "絞り込みを適用" }).click();
  await expect(page.locator(".result-table tbody tr")).toHaveCount(1);
  await expect(page.locator(".result-table tbody tr")).toContainText("58");
  await expect(page.locator(".axis-tick")).not.toHaveCount(0);
  await expect(page.getByLabel("MW 横軸")).toBeVisible();
  await expect(page.getByRole("link", { name: "CSV" })).toHaveAttribute("href", /mw_min=50/);
});
