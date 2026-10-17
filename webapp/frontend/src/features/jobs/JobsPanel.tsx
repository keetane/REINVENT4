import type { Job, JobLogResponse } from "../../api/client";

const statusLabels: Record<Job["status"], string> = {
  queued: "待機中",
  running: "実行中",
  completed: "完了",
  failed: "失敗",
  cancelling: "中断処理中",
  cancelled: "中断済み",
  interrupted: "要確認"
};

function formatDate(value: string | null): string {
  if (!value) return "—";
  return new Intl.DateTimeFormat("ja-JP", {
    dateStyle: "short",
    timeStyle: "medium"
  }).format(new Date(value));
}

function elapsed(job: Job): string {
  if (!job.started_at) return "—";
  const end = job.finished_at ? new Date(job.finished_at).getTime() : Date.now();
  const seconds = Math.max(0, Math.round((end - new Date(job.started_at).getTime()) / 1000));
  if (seconds < 60) return `${seconds}秒`;
  return `${Math.floor(seconds / 60)}分 ${seconds % 60}秒`;
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
  return (
    <div className="jobs-layout">
      <div className="jobs-list-panel">
        <div className="section-heading compact">
          <div>
            <p className="section-kicker">HISTORY</p>
            <h2>ジョブ履歴</h2>
          </div>
          <button className="button quiet" type="button" onClick={onRefresh} disabled={loading}>
            {loading ? "更新中…" : "更新"}
          </button>
        </div>
        <p className="muted">保存済み {total} 件 · 3秒ごとに自動更新</p>
        {error && <p className="notice error" role="alert">{error}</p>}
        {jobs.length === 0 ? (
          <div className="empty-state">まだジョブはありません。</div>
        ) : (
          <div className="job-list">
            {jobs.map((job) => (
              <button
                className={`job-row ${selected?.id === job.id ? "selected" : ""}`}
                type="button"
                key={job.id}
                onClick={() => onSelect(job)}
              >
                <span className={`status-badge ${job.status}`}>{statusLabels[job.status]}</span>
                <span className="job-row-main">
                  <strong>{job.comment || "コメントなし"}</strong>
                  <small>{job.parent_smiles_canonical || job.parent_smiles}</small>
                </span>
                <time>{formatDate(job.created_at)}</time>
              </button>
            ))}
          </div>
        )}
      </div>

      <aside className="job-detail">
        {!selected ? (
          <div className="empty-state">ジョブを選ぶと詳細とログを確認できます。</div>
        ) : (
          <>
            <div className="detail-title">
              <span className={`status-badge ${selected.status}`}>{statusLabels[selected.status]}</span>
              <h3>{selected.comment || "コメントなし"}</h3>
              <code>{selected.id}</code>
            </div>
            <dl className="detail-grid">
              <div><dt>Phase</dt><dd>{selected.progress.message || selected.status}</dd></div>
              <div><dt>開始</dt><dd>{formatDate(selected.started_at)}</dd></div>
              <div><dt>終了</dt><dd>{formatDate(selected.finished_at)}</dd></div>
              <div><dt>実行時間</dt><dd>{elapsed(selected)}</dd></div>
              <div><dt>終了コード</dt><dd>{selected.exit_code ?? "—"}</dd></div>
              <div>
                <dt>進捗</dt>
                <dd>{selected.progress.percent == null ? "計測なし" : `${selected.progress.percent}%`}</dd>
              </div>
            </dl>
            {selected.error_summary && (
              <p className="notice error" role="alert">{selected.error_summary}</p>
            )}
            {selected.status === "completed" && (
              <button className="button primary" type="button" onClick={() => onOpenResults(selected)}>
                結果を開く
              </button>
            )}
            {(selected.status === "queued" || selected.status === "running" || selected.status === "cancelling") && (
              <button
                className="button danger"
                type="button"
                onClick={() => onCancel(selected)}
                disabled={selected.status === "cancelling"}
              >
                {selected.status === "cancelling" ? "中断処理中…" : "このジョブを中断"}
              </button>
            )}
            <div className="log-panel">
              <div className="log-title">
                <strong>ログ末尾</strong>
                {logs?.truncated && <span>一部のみ表示</span>}
              </div>
              <pre>{logs?.lines.length ? logs.lines.join("\n") : "ログはまだありません。"}</pre>
            </div>
          </>
        )}
      </aside>
    </div>
  );
}
