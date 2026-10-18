import { jobStandardizationArtifactUrl, type Job, type JobLogResponse } from "../../api/client";
import { useEffect, useState } from "react";
import { useTranslation } from "react-i18next";

function formatDate(value: string | null, locale: string): string {
  if (!value) return "—";
  return new Intl.DateTimeFormat(locale, {
    dateStyle: "short",
    timeStyle: "medium"
  }).format(new Date(value));
}

function elapsed(job: Job, translate: (key: string, values?: Record<string, number>) => string): string {
  if (!job.started_at) return "—";
  const end = job.finished_at ? new Date(job.finished_at).getTime() : Date.now();
  const seconds = Math.max(0, Math.round((end - new Date(job.started_at).getTime()) / 1000));
  if (seconds < 60) return translate("jobs.seconds", { count: seconds });
  return translate("jobs.minutes", { minutes: Math.floor(seconds / 60), seconds: seconds % 60 });
}

function isGeneration(job: Job): boolean {
  return !job.job_type || job.job_type === "generation";
}

function hasAnalyzableResults(job: Job): boolean {
  return isGeneration(job)
    || job.job_type === "transfer_learning"
    || job.job_type === "reinforcement_learning"
    || job.job_type === "qsar_comparison";
}

interface Props {
  jobs: Job[];
  total: number;
  selected: Job | null;
  logs: JobLogResponse | null;
  loading: boolean;
  error: string | null;
  onRefresh: () => void;
  onSelect: (job: Job) => void;
  onCancel: (job: Job) => void;
  onOpenResults: (job: Job) => void;
}

export function JobsPanel({
  jobs,
  total,
  selected,
  logs,
  loading,
  error,
  onRefresh,
  onSelect,
  onCancel,
  onOpenResults
}: Props) {
  const { t, i18n } = useTranslation();
  const [copyState, setCopyState] = useState<"idle" | "copied" | "failed">("idle");
  const locale = i18n.language.startsWith("en") ? "en-US" : "ja-JP";
  const statusLabel = (job: Job) => t(`jobs.status.${job.status}`);
  const jobTypeLabel = (job: Job) => t(`jobs.type.${job.job_type || "generation"}`);
  const logText = logs?.lines.join("\n") ?? "";

  useEffect(() => {
    setCopyState("idle");
  }, [selected?.id, logText]);

  const copyLogs = async () => {
    if (!logText) return;
    try {
      await navigator.clipboard.writeText(logText);
      setCopyState("copied");
    } catch {
      setCopyState("failed");
    }
  };

  const copyLabel = copyState === "copied"
    ? t("jobs.logCopied")
    : copyState === "failed"
      ? t("jobs.logCopyFailed")
      : t("jobs.copyLogs");

  return (
    <div className="jobs-layout">
      <div className="jobs-list-panel">
        <div className="section-heading compact">
          <div>
            <p className="section-kicker">HISTORY</p>
            <h2>{t("jobs.title")}</h2>
          </div>
          <button className="button quiet" type="button" onClick={onRefresh} disabled={loading}>
            {loading ? t("jobs.refreshing") : t("jobs.refresh")}
          </button>
        </div>
        <p className="muted">{t("jobs.saved", { count: total })}</p>
        {error && <p className="notice error" role="alert">{error}</p>}
        {jobs.length === 0 ? (
          <div className="empty-state">{t("jobs.empty")}</div>
        ) : (
          <div className="job-list">
            {jobs.map((job) => (
              <button
                className={`job-row ${selected?.id === job.id ? "selected" : ""}`}
                type="button"
                key={job.id}
                onClick={() => onSelect(job)}
              >
                <span className={`status-badge ${job.status}`}>{statusLabel(job)}</span>
                <span className="job-row-main">
                  <strong>{job.job_type === "transfer_learning" ? "TL · " : job.job_type === "qsar_training" ? "QSAR · " : job.job_type === "qsar_comparison" ? "QSAR · " : job.job_type === "reinforcement_learning" ? "RL · " : ""}{job.comment || t("jobs.noComment")}</strong>
                  <small>{isGeneration(job) ? (job.parent_smiles_canonical || job.parent_smiles || t("jobs.noParent")) : jobTypeLabel(job)}</small>
                </span>
                <time>{formatDate(job.created_at, locale)}</time>
              </button>
            ))}
          </div>
        )}
      </div>

      <aside className="job-detail">
        {!selected ? (
          <div className="empty-state">{t("jobs.select")}</div>
        ) : (
          <>
            <div className="detail-title">
              <span className={`status-badge ${selected.status}`}>{statusLabel(selected)}</span>
              <h3>{selected.comment || t("jobs.noComment")}</h3>
              <code>{selected.id}</code>
            </div>
            <dl className="detail-grid">
              <div><dt>Phase</dt><dd>{selected.progress.message || selected.status}</dd></div>
              <div><dt>{t("jobs.type")}</dt><dd>{jobTypeLabel(selected)}</dd></div>
              <div><dt>{t("jobs.started")}</dt><dd>{formatDate(selected.started_at, locale)}</dd></div>
              <div><dt>{t("jobs.finished")}</dt><dd>{formatDate(selected.finished_at, locale)}</dd></div>
              <div><dt>{t("jobs.elapsed")}</dt><dd>{elapsed(selected, t)}</dd></div>
              <div><dt>{t("jobs.exitCode")}</dt><dd>{selected.exit_code ?? "—"}</dd></div>
              <div>
                <dt>{t("jobs.progress")}</dt>
                <dd>{selected.progress.percent == null ? t("jobs.notMeasured") : `${selected.progress.percent}%`}</dd>
              </div>
              {selected.job_type === "transfer_learning" && (
                <>
                  <div><dt>Epoch</dt><dd>{selected.progress.current ?? 0} / {selected.progress.total ?? "—"}</dd></div>
                  <div><dt>Training Loss</dt><dd>{selected.metrics.training_loss?.toFixed(4) ?? "—"}</dd></div>
                  <div><dt>Validation Loss</dt><dd>{selected.metrics.validation_loss?.toFixed(4) ?? "—"}</dd></div>
                </>
              )}
              {selected.job_type === "qsar_training" && selected.status === "completed" && (
                <>
                  <div><dt>Test R²</dt><dd>{selected.metrics.test?.r2?.toFixed(4) ?? "—"}</dd></div>
                  <div><dt>Test RMSE</dt><dd>{selected.metrics.test?.rmse?.toFixed(4) ?? "—"}</dd></div>
                  <div><dt>Test ROC-AUC</dt><dd>{selected.metrics.test?.roc_auc?.toFixed(4) ?? "—"}</dd></div>
                  <div><dt>Test F1</dt><dd>{selected.metrics.test?.f1?.toFixed(4) ?? "—"}</dd></div>
                </>
              )}
              {selected.job_type === "reinforcement_learning" && (
                <>
                  <div><dt>Step</dt><dd>{selected.progress.current ?? 0} / {selected.progress.total ?? "—"}</dd></div>
                  <div><dt>Mean Score</dt><dd>{selected.metrics.mean_total_score?.toFixed(4) ?? "—"}</dd></div>
                  <div><dt>Valid SMILES</dt><dd>{selected.metrics.fraction_valid_smiles == null ? "—" : `${(selected.metrics.fraction_valid_smiles * 100).toFixed(1)}%`}</dd></div>
                </>
              )}
            </dl>
            {selected.error_summary && (
              <p className="notice error" role="alert">{selected.error_summary}</p>
            )}
            {selected.status === "completed" && hasAnalyzableResults(selected) && (
              <button className="button primary" type="button" onClick={() => onOpenResults(selected)}>
                {t("jobs.openResults")}
              </button>
            )}
            {["transfer_learning", "qsar_training", "qsar_comparison", "reinforcement_learning"].includes(selected.job_type) && (
              <div className="artifact-links">
                <a className="button quiet" href={jobStandardizationArtifactUrl(selected.id, "accepted")}>標準化 採用CSV</a>
                <a className="button quiet" href={jobStandardizationArtifactUrl(selected.id, "rejected")}>標準化 除外CSV</a>
              </div>
            )}
            {(selected.status === "queued" || selected.status === "running" || selected.status === "cancelling") && (
              <button
                className="button danger"
                type="button"
                onClick={() => onCancel(selected)}
                disabled={selected.status === "cancelling"}
              >
                {selected.status === "cancelling" ? t("jobs.cancelling") : t("jobs.cancel")}
              </button>
            )}
            <div className="log-panel">
              <div className="log-title">
                <div className="log-title-copy">
                  <strong>{t("jobs.logTail")}</strong>
                  {logs?.truncated && <span>{t("jobs.logTruncated")}</span>}
                </div>
                <button
                  className={`log-copy-button ${copyState}`}
                  type="button"
                  onClick={() => void copyLogs()}
                  disabled={!logText}
                >
                  {copyLabel}
                </button>
              </div>
              <pre>{logText || t("jobs.noLogs")}</pre>
            </div>
          </>
        )}
      </aside>
    </div>
  );
}
