import { useCallback, useEffect, useMemo, useState } from "react";

import {
  ApiError,
  cancelJob,
  createJob,
  getHealth,
  getJob,
  getJobLogs,
  listJobs,
  validateParent,
  type HealthResponse,
  type Job,
  type JobLogResponse,
  type ParentValidationResponse
} from "./api/client";
import { JobsPanel } from "./features/jobs/JobsPanel";
import { ParentInput } from "./features/parent-input/ParentInput";
import { ResultsPanel } from "./features/results/ResultsPanel";
import { TutorialPanel } from "./features/tutorial/TutorialPanel";

const POLL_INTERVAL_MS = 3_000;
const selectedJobStorageKey = "reinvent4.selectedJobId";
type ComputeDevice = "cpu" | "cuda" | "mps";

function errorMessage(error: unknown): string {
  if (error instanceof ApiError) return error.message;
  return "API に接続できません。サービスの起動状態を確認してください。";
}

export function App() {
  const [activeView, setActiveView] = useState<"tutorial" | "create" | "jobs" | "results">("create");
  const [health, setHealth] = useState<HealthResponse | null>(null);
  const [healthError, setHealthError] = useState(false);
  const [smiles, setSmiles] = useState("");
  const [comment, setComment] = useState("");
  const [modelKey, setModelKey] = useState("mol2mol_high");
  const [device, setDevice] = useState<ComputeDevice>("mps");
  const [numSmiles, setNumSmiles] = useState(100);
  const [validation, setValidation] = useState<ParentValidationResponse | null>(null);
  const [validationError, setValidationError] = useState<string | null>(null);
  const [createError, setCreateError] = useState<string | null>(null);
  const [validating, setValidating] = useState(false);
  const [creating, setCreating] = useState(false);
  const [jobs, setJobs] = useState<Job[]>([]);
  const [jobTotal, setJobTotal] = useState(0);
  const [jobsLoading, setJobsLoading] = useState(false);
  const [jobsError, setJobsError] = useState<string | null>(null);
  const [selectedJobId, setSelectedJobId] = useState<string | null>(() =>
    window.localStorage.getItem(selectedJobStorageKey)
  );
  const [selectedJob, setSelectedJob] = useState<Job | null>(null);
  const [logs, setLogs] = useState<JobLogResponse | null>(null);

  const refreshJobs = useCallback(async (showLoading = true) => {
    if (showLoading) setJobsLoading(true);
    try {
      const response = await listJobs();
      setJobs(response.items);
      setJobTotal(response.total);
      setJobsError(null);
      return response.items;
    } catch (error) {
      setJobsError(errorMessage(error));
      return [];
    } finally {
      if (showLoading) setJobsLoading(false);
    }
  }, []);

  const refreshSelected = useCallback(async (jobId: string) => {
    try {
      const [job, logResponse] = await Promise.all([getJob(jobId), getJobLogs(jobId)]);
      setSelectedJob(job);
      setLogs(logResponse);
      setJobs((current) => current.map((item) => (item.id === job.id ? job : item)));
    } catch (error) {
      setJobsError(errorMessage(error));
    }
  }, []);

  useEffect(() => {
    const controller = new AbortController();
    getHealth(controller.signal)
      .then(setHealth)
      .catch(() => setHealthError(true));
    void refreshJobs(false);
    return () => controller.abort();
  }, [refreshJobs]);

  useEffect(() => {
    if (!selectedJobId) return;
    window.localStorage.setItem(selectedJobStorageKey, selectedJobId);
    void refreshSelected(selectedJobId);
  }, [refreshSelected, selectedJobId]);

  useEffect(() => {
    const poll = window.setInterval(() => {
      void refreshJobs(false);
      if (selectedJobId) void refreshSelected(selectedJobId);
    }, POLL_INTERVAL_MS);
    return () => window.clearInterval(poll);
  }, [refreshJobs, refreshSelected, selectedJobId]);

  const activeJobs = useMemo(
    () => jobs.filter((job) => ["queued", "running", "cancelling"].includes(job.status)).length,
    [jobs]
  );

  const handleSmilesChange = (value: string) => {
    setSmiles(value);
    if (validation && value.trim() !== validation.canonical_smiles) setValidation(null);
    setValidationError(null);
  };

  const handleValidate = async () => {
    setValidating(true);
    setValidationError(null);
    try {
      const result = await validateParent(smiles, modelKey);
      setValidation(result);
      setSmiles(result.canonical_smiles);
    } catch (error) {
      setValidation(null);
      setValidationError(errorMessage(error));
    } finally {
      setValidating(false);
    }
  };

  const handleCreate = async () => {
    setCreating(true);
    setCreateError(null);
    try {
      const checked = await validateParent(smiles, modelKey);
      setValidation(checked);
      setSmiles(checked.canonical_smiles);
      const job = await createJob({
        parent_smiles: checked.canonical_smiles,
        comment,
        run_conditions: {
          model_key: modelKey,
          num_smiles: numSmiles,
          device,
          randomize_smiles: true,
          unique_molecules: true
        }
      });
      setSelectedJobId(job.id);
      setSelectedJob(job);
      setComment("");
      await refreshJobs(false);
      setActiveView("jobs");
    } catch (error) {
      setCreateError(errorMessage(error));
    } finally {
      setCreating(false);
    }
  };

  const handleSelect = (job: Job) => {
    setSelectedJobId(job.id);
    setSelectedJob(job);
    setLogs(null);
  };

  const handleCancel = async (job: Job) => {
    try {
      const updated = await cancelJob(job.id);
      setSelectedJob(updated);
      await refreshJobs(false);
    } catch (error) {
      setJobsError(errorMessage(error));
    }
  };

  const resultJob = selectedJob?.status === "completed" ? selectedJob : null;

  return (
    <main className="app-shell">
      <header className="topbar">
        <div className="brand">
          <div className="brand-mark" aria-hidden="true">R4</div>
          <div>
            <p className="eyebrow">REINVENT4</p>
            <h1>Molecular Generator</h1>
          </div>
        </div>
        <div className={`service-pill ${healthError ? "offline" : ""}`}>
          <span aria-hidden="true" />
          {health && <span role="status">API v{health.version}</span>}
          {!health && !healthError && <span role="status">接続確認中</span>}
          {healthError && <span role="alert">API offline</span>}
        </div>
      </header>

      <nav className="view-tabs" aria-label="Main navigation">
        <button
          type="button"
          className={activeView === "tutorial" ? "active" : ""}
          onClick={() => setActiveView("tutorial")}
        >
          <span>00</span> Tutorial
        </button>
        <button
          type="button"
          className={activeView === "create" ? "active" : ""}
          onClick={() => setActiveView("create")}
        >
          <span>01</span> ジョブ作成
        </button>
        <button
          type="button"
          className={activeView === "jobs" ? "active" : ""}
          onClick={() => setActiveView("jobs")}
        >
          <span>02</span> 進捗と履歴
          {activeJobs > 0 && <b>{activeJobs}</b>}
        </button>
        <button
          type="button"
          className={activeView === "results" ? "active" : ""}
          onClick={() => resultJob && setActiveView("results")}
          disabled={!resultJob}
        >
          <span>03</span> 結果解析
        </button>
      </nav>

      {activeView === "tutorial" ? (
        <section className="workspace-card tutorial-card">
          <TutorialPanel
            onUseSample={(nextModelKey, nextSmiles) => {
              setModelKey(nextModelKey);
              setSmiles(nextSmiles);
              setValidation(null);
              setValidationError(null);
              setCreateError(null);
              setActiveView("create");
            }}
          />
        </section>
      ) : activeView === "create" ? (
        <section className="workspace-card">
          <div className="section-heading">
            <div>
              <p className="section-kicker">NEW SAMPLING JOB</p>
              <h2>Parent 分子から探索を開始</h2>
              <p>SMILES を直接入力するか、Ketcher で構造を描画できます。</p>
            </div>
            <span className="step-indicator">1 / 2</span>
          </div>

          <ParentInput
            smiles={smiles}
            modelKey={modelKey}
            validation={validation}
            busy={validating}
            error={validationError}
            onSmilesChange={handleSmilesChange}
            onValidate={handleValidate}
          />

          <div className="divider" />

          <div className="conditions-grid">
            <div className="field-group">
              <label htmlFor="model">Model</label>
              <select
                id="model"
                value={modelKey}
                onChange={(event) => {
                  setModelKey(event.target.value);
                  setValidation(null);
                  setValidationError(null);
                }}
              >
                <option value="mol2mol_high">Mol2Mol high similarity</option>
                <option value="mol2mol_medium">Mol2Mol medium</option>
                <option value="mol2mol_mmp">Mol2Mol MMP</option>
                <option value="mol2mol_scaffold">Mol2Mol scaffold</option>
                <option value="libinvent">LibInvent</option>
                <option value="reinvent">REINVENT</option>
                <option value="linkinvent">LinkInvent</option>
              </select>
            </div>
            <div className="field-group">
              <label htmlFor="num-smiles">生成件数</label>
              <input
                id="num-smiles"
                type="number"
                min={1}
                max={10_000}
                value={numSmiles}
                onChange={(event) => setNumSmiles(Number(event.target.value))}
              />
            </div>
            <div className="field-group">
              <label htmlFor="device">実行デバイス</label>
              <select
                id="device"
                value={device}
                onChange={(event) => setDevice(event.target.value as ComputeDevice)}
              >
                <option value="cpu">CPU</option>
                <option value="cuda">GPU</option>
                <option value="mps">MPS</option>
              </select>
              <p className="field-hint">Apple SiliconではMPSを使用します。</p>
            </div>
            <div className="field-group wide">
              <label htmlFor="comment">コメント</label>
              <input
                id="comment"
                value={comment}
                maxLength={1_000}
                onChange={(event) => setComment(event.target.value)}
                placeholder="例: kinase A / first pass"
              />
              <p className="field-hint">表示名と安全な結果ファイル名に使用します。</p>
            </div>
          </div>

          {createError && <p className="notice error" role="alert">{createError}</p>}
          <div className="submit-row">
            <p>ジョブ開始後は、この画面を閉じても worker 上で処理が続きます。</p>
            <button
              className="button primary large"
              type="button"
              onClick={handleCreate}
              disabled={creating || !smiles.trim()}
            >
              {creating ? "登録中…" : "ジョブを開始"}
            </button>
          </div>
        </section>
      ) : activeView === "jobs" ? (
        <section className="workspace-card jobs-card">
          <JobsPanel
            jobs={jobs}
            total={jobTotal}
            selected={selectedJob}
            logs={logs}
            loading={jobsLoading}
            error={jobsError}
            onRefresh={() => void refreshJobs()}
            onSelect={handleSelect}
            onCancel={(job) => void handleCancel(job)}
            onOpenResults={(job) => {
              setSelectedJobId(job.id);
              setSelectedJob(job);
              setActiveView("results");
            }}
          />
        </section>
      ) : resultJob ? (
        <section className="workspace-card results-card">
          <ResultsPanel job={resultJob} onBack={() => setActiveView("jobs")} />
        </section>
      ) : (
        <section className="workspace-card"><div className="empty-state">完了したジョブを履歴から選択してください。</div></section>
      )}

      <footer>
        <span>Local-only workspace</span>
        <span>FastAPI · SQLite · independent worker</span>
      </footer>
    </main>
  );
}
