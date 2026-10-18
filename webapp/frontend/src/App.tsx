import { useCallback, useEffect, useMemo, useState } from "react";
import { useTranslation } from "react-i18next";

import "./i18n";

import {
  ApiError,
  cancelJob,
  createJob,
  getHealth,
  getJob,
  getJobLogs,
  deleteModel,
  listJobs,
  listModels,
  validateParent,
  type HealthResponse,
  type Job,
  type JobLogResponse,
  type ModelInfo,
  type ParentValidationResponse
} from "./api/client";
import { JobsPanel } from "./features/jobs/JobsPanel";
import { ParentInput } from "./features/parent-input/ParentInput";
import { ResultsPanel } from "./features/results/ResultsPanel";
import { QsarPanel } from "./features/qsar/QsarPanel";
import { RlPanel } from "./features/rl/RlPanel";
import { TutorialPanel } from "./features/tutorial/TutorialPanel";
import { TransferLearningPanel } from "./features/training/TransferLearningPanel";

const POLL_INTERVAL_MS = 3_000;
const selectedJobStorageKey = "reinvent4.selectedJobId";
type ComputeDevice = "cpu" | "cuda" | "mps";

function errorMessage(error: unknown): string {
  if (error instanceof ApiError) return error.message;
  return "Unable to connect to the API.";
}

export function App() {
  const { t, i18n } = useTranslation();
  const [activeView, setActiveView] = useState<"tutorial" | "create" | "training" | "qsar" | "rl" | "jobs" | "results">("create");
  const [health, setHealth] = useState<HealthResponse | null>(null);
  const [healthError, setHealthError] = useState(false);
  const [smiles, setSmiles] = useState("");
  const [comment, setComment] = useState("");
  const [modelId, setModelId] = useState("standard:reinvent");
  const [models, setModels] = useState<ModelInfo[]>([]);
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

  const refreshModels = useCallback(async () => {
    try {
      const response = await listModels();
      setModels(response.items);
      setJobsError(null);
      return response.items;
    } catch (error) {
      setJobsError(errorMessage(error));
      return [];
    }
  }, []);

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
    void refreshModels();
    return () => controller.abort();
  }, [refreshJobs, refreshModels]);

  useEffect(() => {
    if (!selectedJobId) return;
    window.localStorage.setItem(selectedJobStorageKey, selectedJobId);
    void refreshSelected(selectedJobId);
  }, [refreshSelected, selectedJobId]);

  useEffect(() => {
    const poll = window.setInterval(() => {
      void refreshJobs(false);
      void refreshModels();
      if (selectedJobId) void refreshSelected(selectedJobId);
    }, POLL_INTERVAL_MS);
    return () => window.clearInterval(poll);
  }, [refreshJobs, refreshModels, refreshSelected, selectedJobId]);

  const selectedModel = useMemo(
    () => models.find((model) => model.id === modelId) ?? models[0] ?? null,
    [modelId, models]
  );
  const parentModelKey = selectedModel?.base_model_key ?? "reinvent";
  const parentRequired = selectedModel ? selectedModel.model_type !== "reinvent" : false;
  const selectedModelDescription = selectedModel
    ? selectedModel.source === "standard"
      ? t(`generation.modelDescription.${selectedModel.base_model_key}`)
      : t(`generation.modelDescription.${selectedModel.source}`, { type: selectedModel.model_type.toUpperCase() })
    : t("generation.loadingModels");

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
      const result = await validateParent(smiles, parentModelKey);
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
      const checked = parentRequired ? await validateParent(smiles, parentModelKey) : null;
      if (checked) {
        setValidation(checked);
        setSmiles(checked.canonical_smiles);
      }
      const job = await createJob({
        parent_smiles: checked?.canonical_smiles ?? "",
        comment,
        run_conditions: {
          model_id: selectedModel?.id ?? modelId,
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

  const handleTrainingCreated = async (job: Job) => {
    setSelectedJobId(job.id);
    setSelectedJob(job);
    await refreshJobs(false);
    await refreshModels();
    setActiveView("jobs");
  };

  const handleDeleteModel = async (targetModelId: string) => {
    try {
      await deleteModel(targetModelId);
      if (modelId === targetModelId) setModelId("standard:reinvent");
      await refreshModels();
    } catch (error) {
      setJobsError(errorMessage(error));
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

  const resultJob = selectedJob?.status === "completed" && (
    !selectedJob.job_type
    || selectedJob.job_type === "generation"
    || selectedJob.job_type === "transfer_learning"
    || selectedJob.job_type === "reinforcement_learning"
    || selectedJob.job_type === "qsar_comparison"
  ) ? selectedJob : null;

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
        <div className="topbar-actions"><label className="language-select"><span>{t("language.label")}</span><select value={i18n.language.startsWith("en") ? "en" : "ja"} onChange={(event) => void i18n.changeLanguage(event.target.value)}><option value="ja">日本語</option><option value="en">English</option></select></label><div className={`service-pill ${healthError ? "offline" : ""}`}>
          <span aria-hidden="true" />
          {health && <span role="status">API v{health.version}</span>}
          {!health && !healthError && <span role="status">{t("service.checking")}</span>}
          {healthError && <span role="alert">API offline</span>}
        </div></div>
      </header>

      <nav className="view-tabs" aria-label="Main navigation">
        <button
          type="button"
          className={activeView === "tutorial" ? "active" : ""}
          onClick={() => setActiveView("tutorial")}
        >
          <span>00</span> {t("nav.tutorial")}
        </button>
        <button
          type="button"
          className={activeView === "create" ? "active" : ""}
          onClick={() => setActiveView("create")}
        >
          <span>01</span> {t("nav.generate")}
        </button>
        <button
          type="button"
          className={activeView === "training" ? "active" : ""}
          onClick={() => setActiveView("training")}
        >
          <span>02</span> {t("nav.training")}
        </button>
        <button
          type="button"
          className={activeView === "qsar" ? "active" : ""}
          onClick={() => setActiveView("qsar")}
        >
          <span>03</span> {t("nav.qsar")}
        </button>
        <button
          type="button"
          className={activeView === "rl" ? "active" : ""}
          onClick={() => setActiveView("rl")}
        >
          <span>04</span> {t("nav.rl")}
        </button>
        <button
          type="button"
          className={activeView === "jobs" ? "active" : ""}
          onClick={() => setActiveView("jobs")}
        >
          <span>05</span> {t("nav.jobs")}
          {activeJobs > 0 && <b>{activeJobs}</b>}
        </button>
        <button
          type="button"
          className={activeView === "results" ? "active" : ""}
          onClick={() => resultJob && setActiveView("results")}
          disabled={!resultJob}
        >
          <span>06</span> {t("nav.results")}
        </button>
      </nav>

      {activeView === "tutorial" ? (
        <section className="workspace-card tutorial-card">
          <TutorialPanel
            onUseSample={(nextModelKey, nextSmiles) => {
              setModelId(`standard:${nextModelKey}`);
              setSmiles(nextSmiles);
              setValidation(null);
              setValidationError(null);
              setCreateError(null);
              setActiveView("create");
            }}
            onOpenTraining={() => setActiveView("training")}
          />
        </section>
      ) : activeView === "create" ? (
        <section className="workspace-card">
          <div className="section-heading">
            <div>
              <p className="section-kicker">NEW SAMPLING JOB</p>
              <h2>{parentRequired ? t("generation.parentTitle") : t("generation.priorTitle")}</h2>
              <p>{parentRequired ? t("generation.parentRequired") : t("generation.parentOptional")}</p>
            </div>
            <span className="step-indicator">1 / 2</span>
          </div>

          <div className="model-selector-panel">
            <div className="model-selector-row">
              <div className="field-group">
                <label htmlFor="model">Model</label>
                <select
                  id="model"
                  value={selectedModel?.id ?? modelId}
                  onChange={(event) => {
                    setModelId(event.target.value);
                    setValidation(null);
                    setValidationError(null);
                  }}
                >
                  <optgroup label={t("generation.standardModels")}>
                    {models.filter((model) => model.source === "standard").map((model) => (
                      <option value={model.id} key={model.id}>{model.name}</option>
                    ))}
                  </optgroup>
                  {models.some((model) => model.source === "transfer_learning") && (
                    <optgroup label={t("generation.transferModels")}>
                      {models.filter((model) => model.source === "transfer_learning").map((model) => (
                        <option value={model.id} key={model.id}>{model.name}</option>
                      ))}
                    </optgroup>
                  )}
                  {models.some((model) => model.source === "reinforcement_learning") && (
                    <optgroup label={t("generation.rlModels")}>
                      {models.filter((model) => model.source === "reinforcement_learning").map((model) => (
                        <option value={model.id} key={model.id}>{model.name}</option>
                      ))}
                    </optgroup>
                  )}
                </select>
                <p className="model-description" role="status">
                  {selectedModelDescription}
                  {parentRequired ? t("generation.requiresParent") : t("generation.noParent")}
                </p>
              </div>
              <button
                className="button primary large model-start-button"
                type="button"
                onClick={handleCreate}
                disabled={creating || !selectedModel || (parentRequired && !smiles.trim())}
              >
                {creating ? t("generation.starting") : t("generation.start")}
              </button>
            </div>
          </div>

          {parentRequired ? (
            <ParentInput
              smiles={smiles}
              modelKey={parentModelKey}
              validation={validation}
              busy={validating}
              error={validationError}
              onSmilesChange={handleSmilesChange}
              onValidate={handleValidate}
            />
          ) : (
            <div className="parent-not-required">{t("generation.noParentDetail")}</div>
          )}

          <div className="divider" />

          <div className="conditions-grid">
            <div className="field-group">
              <label htmlFor="num-smiles">{t("generation.count")}</label>
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
              <label htmlFor="device">{t("generation.device")}</label>
              <select
                id="device"
                value={device}
                onChange={(event) => setDevice(event.target.value as ComputeDevice)}
              >
                <option value="cpu">CPU</option>
                <option value="cuda">GPU</option>
                <option value="mps">MPS</option>
              </select>
              <p className="field-hint">{t("generation.mpsHint")}</p>
            </div>
            <div className="field-group wide">
              <label htmlFor="comment">{t("generation.comment")}</label>
              <input
                id="comment"
                value={comment}
                maxLength={1_000}
                onChange={(event) => setComment(event.target.value)}
                placeholder={t("generation.commentPlaceholder")}
              />
              <p className="field-hint">{t("generation.commentHint")}</p>
            </div>
          </div>

          {createError && <p className="notice error" role="alert">{createError}</p>}
          <p className="job-persistence-note">{t("generation.persistence")}</p>
        </section>
      ) : activeView === "training" ? (
        <section className="workspace-card">
          <TransferLearningPanel
            models={models}
            onCreated={(job) => void handleTrainingCreated(job)}
            onDeleteModel={handleDeleteModel}
          />
        </section>
      ) : activeView === "qsar" ? (
        <section className="workspace-card">
          <QsarPanel onCreated={(job) => void handleTrainingCreated(job)} />
        </section>
      ) : activeView === "rl" ? (
        <section className="workspace-card">
          <RlPanel models={models} onCreated={(job) => void handleTrainingCreated(job)} />
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
        <section className="workspace-card"><div className="empty-state">{t("results.selectCompleted")}</div></section>
      )}

      <footer>
        <span>Local-only workspace</span>
        <span>FastAPI · SQLite · independent worker</span>
      </footer>
    </main>
  );
}
