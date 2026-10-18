import { useEffect, useMemo, useState } from "react";
import { useTranslation } from "react-i18next";

import {
  ApiError,
  getQsarComparison,
  getResultMolecules,
  getResultSchema,
  getResultStatistics,
  resultDownloadUrl,
  type DescriptorSchema,
  type Job,
  type MoleculeListResponse,
  type QsarComparisonResult,
  type ResultQuery,
  type ResultSchemaResponse,
  type StatisticsResponse
} from "../../api/client";

const CORE_DESCRIPTORS = ["MW", "RotB", "HBD", "HBA", "LogP", "TPSA"] as const;
const ADDED_DESCRIPTORS = [
  "QED",
  "SAScore",
  "FractionCSP3",
  "RingCount",
  "AromaticRingCount",
  "AliphaticRingCount",
  "HeavyAtomCount",
  "FormalCharge",
  "LipinskiViolations",
  "MolMR",
  "ParentSimilarity",
  "ScaffoldSimilarity"
] as const;
const FILTER_DESCRIPTORS = [...CORE_DESCRIPTORS, ...ADDED_DESCRIPTORS] as const;
const QUERY_NAMES: Record<(typeof FILTER_DESCRIPTORS)[number], string> = {
  MW: "mw",
  RotB: "rotb",
  HBD: "hbd",
  HBA: "hba",
  LogP: "logp",
  TPSA: "tpsa",
  QED: "qed",
  SAScore: "sa_score",
  FractionCSP3: "fraction_csp3",
  RingCount: "ring_count",
  AromaticRingCount: "aromatic_ring_count",
  AliphaticRingCount: "aliphatic_ring_count",
  HeavyAtomCount: "heavy_atom_count",
  FormalCharge: "formal_charge",
  LipinskiViolations: "lipinski_violations",
  MolMR: "molmr",
  ParentSimilarity: "parent_similarity",
  ScaffoldSimilarity: "scaffold_similarity"
};
const PAGE_SIZE = 50;

type Bounds = Record<string, { minimum: string; maximum: string }>;

function blankBounds(): Bounds {
  return Object.fromEntries(
    FILTER_DESCRIPTORS.map((name) => [name, { minimum: "", maximum: "" }])
  );
}

function storedBounds(jobId: string): Bounds {
  try {
    const parsed = JSON.parse(window.sessionStorage.getItem(`reinvent4.results.${jobId}`) ?? "null") as Bounds | null;
    if (parsed && FILTER_DESCRIPTORS.every((name) => typeof parsed[name]?.minimum === "string" && typeof parsed[name]?.maximum === "string")) {
      return parsed;
    }
  } catch {
    // Ignore malformed browser state and return to the unfiltered view.
  }
  return blankBounds();
}

function message(error: unknown, fallback: string): string {
  return error instanceof ApiError ? error.message : fallback;
}

function filtersToQuery(bounds: Bounds): ResultQuery {
  const query: ResultQuery = {};
  FILTER_DESCRIPTORS.forEach((descriptor) => {
    const key = QUERY_NAMES[descriptor];
    query[`${key}_min`] = bounds[descriptor].minimum;
    query[`${key}_max`] = bounds[descriptor].maximum;
  });
  return query;
}

function formatNumber(value: number | null | undefined, digits = 2): string {
  if (value === null || value === undefined) return "—";
  return new Intl.NumberFormat("ja-JP", { maximumFractionDigits: digits }).format(value);
}

function clamp(value: number, minimum: number, maximum: number): number {
  return Math.min(maximum, Math.max(minimum, value));
}

function sliderValue(value: number): string {
  return String(Number(value.toFixed(6)));
}

function axisTicks(minimum: number, maximum: number, count = 5): number[] {
  if (minimum === maximum) return [minimum];
  return Array.from({ length: count }, (_, index) => minimum + ((maximum - minimum) * index) / (count - 1));
}

function DescriptorFilter({
  descriptor,
  bounds,
  limits,
  onChange
}: {
  descriptor: (typeof FILTER_DESCRIPTORS)[number];
  bounds: { minimum: string; maximum: string };
  limits?: DescriptorSchema;
  onChange: (next: { minimum: string; maximum: string }) => void;
}) {
  const { t } = useTranslation();
  const rangeMinimum = limits?.minimum;
  const rangeMaximum = limits?.maximum;
  const hasRange = typeof rangeMinimum === "number" && typeof rangeMaximum === "number" && rangeMaximum > rangeMinimum;
  const minimumValue = typeof rangeMinimum === "number" ? rangeMinimum : 0;
  const maximumValue = typeof rangeMaximum === "number" ? rangeMaximum : minimumValue;
  const lower = hasRange
    ? clamp(bounds.minimum === "" ? minimumValue : Number(bounds.minimum), minimumValue, maximumValue)
    : 0;
  const upper = hasRange
    ? clamp(bounds.maximum === "" ? maximumValue : Number(bounds.maximum), minimumValue, maximumValue)
    : 0;
  const safeLower = Number.isFinite(lower) ? Math.min(lower, upper) : minimumValue;
  const safeUpper = Number.isFinite(upper) ? Math.max(upper, lower) : maximumValue;
  const width = hasRange ? maximumValue - minimumValue : 1;
  const start = ((safeLower - minimumValue) / width) * 100;
  const end = ((safeUpper - minimumValue) / width) * 100;
  const step = hasRange ? Math.max(width / 200, 0.0001) : 1;

  return (
    <fieldset className="descriptor-filter">
      <legend>{descriptor}</legend>
      <div className="filter-number-row">
        <input
          aria-label={t("results.minimum", { descriptor })}
          inputMode="decimal"
          placeholder="min"
          value={bounds.minimum}
          onChange={(event) => onChange({ ...bounds, minimum: event.target.value })}
        />
        <span>–</span>
        <input
          aria-label={t("results.maximum", { descriptor })}
          inputMode="decimal"
          placeholder="max"
          value={bounds.maximum}
          onChange={(event) => onChange({ ...bounds, maximum: event.target.value })}
        />
      </div>
      {hasRange ? (
        <div
          className="range-slider"
          style={{ background: `linear-gradient(to right, #294138 0%, #294138 ${start}%, #55d69d ${start}%, #55d69d ${end}%, #294138 ${end}%, #294138 100%)` }}
        >
          <input
            type="range"
            aria-label={t("results.lowerSlider", { descriptor })}
            min={minimumValue}
            max={maximumValue}
            step={step}
            value={safeLower}
            onChange={(event) => onChange({ ...bounds, minimum: sliderValue(Math.min(Number(event.target.value), safeUpper)) })}
          />
          <input
            type="range"
            aria-label={t("results.upperSlider", { descriptor })}
            min={minimumValue}
            max={maximumValue}
            step={step}
            value={safeUpper}
            onChange={(event) => onChange({ ...bounds, maximum: sliderValue(Math.max(Number(event.target.value), safeLower)) })}
          />
        </div>
      ) : (
        <div className="range-static">{rangeMinimum === null || rangeMinimum === undefined ? t("results.rangeLoading") : formatNumber(rangeMinimum)}</div>
      )}
    </fieldset>
  );
}

function ScatterPlot({ statistics }: { statistics: StatisticsResponse }) {
  const { t } = useTranslation();
  const width = 560;
  const height = 290;
  const left = 58;
  const right = 18;
  const top = 18;
  const bottom = 46;
  const points = statistics.scatter;
  if (!points.length) return <div className="chart-empty">{t("results.noPoints")}</div>;
  const xs = points.map((point) => point.x);
  const ys = points.map((point) => point.y);
  const xMin = Math.min(...xs);
  const xMax = Math.max(...xs);
  const yMin = Math.min(...ys);
  const yMax = Math.max(...ys);
  const x = (value: number) => left + ((value - xMin) / (xMax - xMin || 1)) * (width - left - right);
  const y = (value: number) => height - bottom - ((value - yMin) / (yMax - yMin || 1)) * (height - top - bottom);
  const xTicks = axisTicks(xMin, xMax);
  const yTicks = axisTicks(yMin, yMax);

  return (
    <svg className="scatter-plot" viewBox={`0 0 ${width} ${height}`} role="img" aria-label={t("results.scatterLabel", { x: statistics.scatter_x, y: statistics.scatter_y })}>
      <line className="axis-line" x1={left} y1={height - bottom} x2={width - right} y2={height - bottom} />
      <line className="axis-line" x1={left} y1={top} x2={left} y2={height - bottom} />
      {xTicks.map((tick) => (
        <g className="axis-tick" key={`x-${tick}`}>
          <line x1={x(tick)} y1={height - bottom} x2={x(tick)} y2={height - bottom + 5} />
          <text x={x(tick)} y={height - bottom + 18} textAnchor="middle">{formatNumber(tick)}</text>
        </g>
      ))}
      {yTicks.map((tick) => (
        <g className="axis-tick" key={`y-${tick}`}>
          <line x1={left - 5} y1={y(tick)} x2={left} y2={y(tick)} />
          <text x={left - 9} y={y(tick) + 4} textAnchor="end">{formatNumber(tick)}</text>
        </g>
      ))}
      {points.map((point) => (
        <circle key={point.row_number} cx={x(point.x)} cy={y(point.y)} r="3.5">
          <title>{`#${point.row_number}: ${formatNumber(point.x)} / ${formatNumber(point.y)}`}</title>
        </circle>
      ))}
      <text className="axis-label" x={(left + width - right) / 2} y={height - 5} textAnchor="middle">{statistics.scatter_x}</text>
      <text className="axis-label" x="14" y={height / 2} textAnchor="middle" transform={`rotate(-90 14 ${height / 2})`}>{statistics.scatter_y}</text>
    </svg>
  );
}

interface Props {
  job: Job;
  onBack: () => void;
}

type LearningPoint = { index: number; values: Record<string, number> };

function finiteNumber(value: unknown): number | null {
  return typeof value === "number" && Number.isFinite(value) ? value : null;
}

function learningPoints(job: Job): LearningPoint[] {
  const rawHistory = Array.isArray(job.metrics.history) ? job.metrics.history : [];
  const points = rawHistory.flatMap((entry, position) => {
    const index = finiteNumber(entry[job.job_type === "transfer_learning" ? "epoch" : "step"]);
    if (index === null) return [];
    const values: Record<string, number> = {};
    const candidates: Record<string, unknown> = job.job_type === "transfer_learning"
      ? { training_loss: entry.training_loss, validation_loss: entry.validation_loss }
      : {
          mean_total_score: entry.mean_total_score,
          fraction_valid_smiles: entry.fraction_valid_smiles,
          ...(typeof entry.components === "object" && entry.components ? entry.components as Record<string, unknown> : {})
        };
    Object.entries(candidates).forEach(([key, value]) => {
      const numeric = finiteNumber(value);
      if (numeric !== null) values[key] = numeric;
    });
    return [{ index, values, position }];
  });
  if (points.length) return points.sort((a, b) => a.index - b.index || a.position - b.position);
  const index = job.progress.current ?? job.progress.total ?? 1;
  const values: Record<string, number> = {};
  const candidates = job.job_type === "transfer_learning"
    ? { training_loss: job.metrics.training_loss, validation_loss: job.metrics.validation_loss }
    : {
        mean_total_score: job.metrics.mean_total_score,
        fraction_valid_smiles: job.metrics.fraction_valid_smiles,
        ...(job.metrics.components ?? {})
      };
  Object.entries(candidates).forEach(([key, value]) => {
    const numeric = finiteNumber(value);
    if (numeric !== null) values[key] = numeric;
  });
  return Object.keys(values).length ? [{ index, values }] : [];
}

const SERIES_COLORS = ["#55d69d", "#f1c75b", "#66a8ff", "#ed82b1", "#b79cff", "#ff8e66"];

function LearningChart({ points, series, xLabel }: { points: LearningPoint[]; series: string[]; xLabel: string }) {
  const width = 760;
  const height = 330;
  const left = 64;
  const right = 24;
  const top = 24;
  const bottom = 54;
  const allValues = points.flatMap((point) => series.flatMap((name) => point.values[name] === undefined ? [] : [point.values[name]]));
  if (!points.length || !allValues.length) return <div className="chart-empty">データがありません</div>;
  const xMin = Math.min(...points.map((point) => point.index));
  const xMax = Math.max(...points.map((point) => point.index));
  const yMinRaw = Math.min(...allValues);
  const yMaxRaw = Math.max(...allValues);
  const padding = (yMaxRaw - yMinRaw || Math.max(Math.abs(yMaxRaw), 1)) * 0.08;
  const yMin = yMinRaw - padding;
  const yMax = yMaxRaw + padding;
  const x = (value: number) => left + ((value - xMin) / (xMax - xMin || 1)) * (width - left - right);
  const y = (value: number) => height - bottom - ((value - yMin) / (yMax - yMin || 1)) * (height - top - bottom);
  const pathFor = (name: string) => points.flatMap((point) => point.values[name] === undefined ? [] : [point])
    .map((point, index) => `${index ? "L" : "M"}${x(point.index)},${y(point.values[name])}`).join(" ");
  return (
    <div className="learning-chart-wrap">
      <svg className="learning-chart" viewBox={`0 0 ${width} ${height}`} role="img">
        {axisTicks(yMin, yMax).map((tick) => <g key={tick}><line x1={left} y1={y(tick)} x2={width - right} y2={y(tick)} /><text x={left - 8} y={y(tick) + 4} textAnchor="end">{formatNumber(tick, 4)}</text></g>)}
        <line className="axis-line" x1={left} y1={height - bottom} x2={width - right} y2={height - bottom} />
        {series.map((name, index) => <path key={name} d={pathFor(name)} fill="none" stroke={SERIES_COLORS[index % SERIES_COLORS.length]} strokeWidth="2.5" />)}
        <text className="axis-label" x={(left + width - right) / 2} y={height - 12} textAnchor="middle">{xLabel}</text>
      </svg>
      <div className="chart-legend">{series.map((name, index) => <span key={name}><i style={{ background: SERIES_COLORS[index % SERIES_COLORS.length] }} />{name}</span>)}</div>
    </div>
  );
}

function slope(values: number[]): number {
  if (values.length < 2) return 0;
  return (values.at(-1)! - values[0]) / (values.length - 1);
}

function LearningResultsPanel({ job, onBack }: Props) {
  const { t } = useTranslation();
  const points = learningPoints(job);
  const isTl = job.job_type === "transfer_learning";
  const primarySeries = isTl
    ? ["training_loss", "validation_loss"].filter((name) => points.some((point) => point.values[name] !== undefined))
    : ["mean_total_score", "fraction_valid_smiles"].filter((name) => points.some((point) => point.values[name] !== undefined));
  const componentSeries = isTl ? [] : Array.from(new Set(points.flatMap((point) => Object.keys(point.values))))
    .filter((name) => !primarySeries.includes(name) && name !== "total_score").slice(0, 6);
  const recent = points.slice(-Math.min(5, points.length));
  let diagnostic = t("learning.insufficient");
  let diagnosticClass = "warning";
  if (points.length >= 3 && isTl) {
    const train = recent.flatMap((point) => point.values.training_loss === undefined ? [] : [point.values.training_loss]);
    const validation = recent.flatMap((point) => point.values.validation_loss === undefined ? [] : [point.values.validation_loss]);
    if (train.length >= 2 && validation.length >= 2 && slope(train) < 0 && slope(validation) > 0) {
      diagnostic = t("learning.overfitRisk");
      diagnosticClass = "error";
    } else if (validation.length >= 2 && Math.abs(slope(validation)) <= Math.max(Math.abs(validation.at(-1)!) * 0.01, 0.0001)) {
      diagnostic = t("learning.converged");
      diagnosticClass = "success";
    } else {
      diagnostic = t("learning.progressing");
    }
  } else if (points.length >= 3 && !isTl) {
    const scores = recent.flatMap((point) => point.values.mean_total_score === undefined ? [] : [point.values.mean_total_score]);
    const validity = recent.flatMap((point) => point.values.fraction_valid_smiles === undefined ? [] : [point.values.fraction_valid_smiles]);
    if (validity.length && validity.at(-1)! < 0.8) {
      diagnostic = t("learning.validityRisk");
      diagnosticClass = "error";
    } else if (scores.length >= 2 && Math.abs(slope(scores)) <= Math.max(Math.abs(scores.at(-1)!) * 0.01, 0.001)) {
      diagnostic = t("learning.converged");
      diagnosticClass = "success";
    } else {
      diagnostic = t("learning.progressing");
    }
  }
  return (
    <div className="results-panel learning-results">
      <div className="section-heading result-heading"><div><p className="section-kicker">LEARNING ANALYSIS</p><h2>{job.comment || (isTl ? "Transfer Learning" : t("jobs.type.reinforcement_learning"))}</h2><p>{isTl ? t("learning.tlDescription") : t("learning.rlDescription")}</p></div><button className="button quiet" type="button" onClick={onBack}>{t("results.back")}</button></div>
      <div className="result-summary-strip learning-summary">
        <div><strong>{points.length}</strong><span>{isTl ? "Epochs" : "Steps"}</span></div>
        <div><strong>{formatNumber(isTl ? job.metrics.training_loss : job.metrics.mean_total_score, 4)}</strong><span>{isTl ? "Training Loss" : "Mean Score"}</span></div>
        <div><strong>{formatNumber(isTl ? job.metrics.validation_loss : job.metrics.fraction_valid_smiles, 4)}</strong><span>{isTl ? "Validation Loss" : "Valid SMILES"}</span></div>
      </div>
      <p className={`notice ${diagnosticClass}`}>{diagnostic}</p>
      <p className="muted learning-caveat">{t("learning.heuristic")}</p>
      <section className="analytics-card learning-analysis-card"><h3>{isTl ? t("learning.lossCurve") : t("learning.scoreCurve")}</h3><LearningChart points={points} series={primarySeries} xLabel={isTl ? "Epoch" : "Step"} /></section>
      {componentSeries.length > 0 && <section className="analytics-card learning-analysis-card"><h3>{t("learning.components")}</h3><LearningChart points={points} series={componentSeries} xLabel="Step" /></section>}
      <section className="analytics-card learning-analysis-card"><h3>{t("learning.history")}</h3><div className="learning-table-wrap"><table className="learning-table"><thead><tr><th>{isTl ? "Epoch" : "Step"}</th>{[...primarySeries, ...componentSeries].map((name) => <th key={name}>{name}</th>)}</tr></thead><tbody>{points.map((point) => <tr key={point.index}><td>{point.index}</td>{[...primarySeries, ...componentSeries].map((name) => <td key={name}>{formatNumber(point.values[name], 5)}</td>)}</tr>)}</tbody></table></div></section>
    </div>
  );
}

type ComparisonSort = "metric" | "algorithm";

function sortComparisonResults(comparison: QsarComparisonResult, sort: ComparisonSort) {
  const results = [...comparison.results];
  if (sort === "algorithm") return results.sort((a, b) => a.algorithm.localeCompare(b.algorithm));
  const lowerIsBetter = ["rmse", "mae"].includes(comparison.primary_metric);
  return results.sort((a, b) => {
    const aValue = finiteNumber(a.metrics?.cross_validation?.[comparison.primary_metric]);
    const bValue = finiteNumber(b.metrics?.cross_validation?.[comparison.primary_metric]);
    if (aValue === null) return bValue === null ? 0 : 1;
    if (bValue === null) return -1;
    return lowerIsBetter ? aValue - bValue : bValue - aValue;
  });
}

function ComparisonPlot({ name, points }: { name: string; points: number[][] }) {
  const finitePoints = points.filter((point) => point.length >= 2 && Number.isFinite(point[0]) && Number.isFinite(point[1]));
  if (!finitePoints.length) return null;
  const xs = finitePoints.map((point) => point[0]);
  const ys = finitePoints.map((point) => point[1]);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const mapX = (value: number) => 12 + ((value - minX) / (maxX - minX || 1)) * 176;
  const mapY = (value: number) => 88 - ((value - minY) / (maxY - minY || 1)) * 76;
  return (
    <figure>
      <figcaption>{name.replaceAll("_", " ")}</figcaption>
      <svg viewBox="0 0 200 100" role="img" aria-label={name}>
        <line x1="12" y1="88" x2="188" y2="88" />
        <line x1="12" y1="12" x2="12" y2="88" />
        {finitePoints.map((point, index) => <circle key={index} cx={mapX(point[0])} cy={mapY(point[1])} r="2.5" />)}
      </svg>
    </figure>
  );
}

function ComparisonCharts({ comparison }: { comparison: QsarComparisonResult }) {
  const { t } = useTranslation();
  const completed = comparison.results.filter((result) => result.status === "completed");
  const [selectedAlgorithm, setSelectedAlgorithm] = useState(comparison.recommended_algorithm);
  useEffect(() => setSelectedAlgorithm(comparison.recommended_algorithm), [comparison]);
  const values = completed.flatMap((result) => {
    const value = finiteNumber(result.metrics?.cross_validation?.[comparison.primary_metric]);
    return value === null ? [] : [{ name: result.algorithm, value }];
  });
  const maximum = Math.max(...values.map((item) => Math.abs(item.value)), 1);
  const selected = completed.find((result) => result.algorithm === selectedAlgorithm) ?? completed[0];
  const folds = selected?.metrics?.cross_validation?.[`${comparison.primary_metric}_values`];
  return (
    <div className="comparison-charts">
      <label className="comparison-model-selector">{t("qsar.distributionModel")}
        <select value={selected?.algorithm ?? ""} onChange={(event) => setSelectedAlgorithm(event.target.value)}>
          {completed.map((result) => <option key={result.algorithm} value={result.algorithm}>{result.algorithm}{result.algorithm === comparison.recommended_algorithm ? ` (${t("qsar.recommendedShort")})` : ""}</option>)}
        </select>
      </label>
      <figure>
        <figcaption>Cross-validation {comparison.primary_metric}</figcaption>
        <div className="comparison-bars">
          {values.map((item) => <div key={item.name}><span>{item.name}</span><i style={{ width: `${Math.max(2, Math.abs(item.value) / maximum * 100)}%` }} /><b>{item.value.toFixed(4)}</b></div>)}
        </div>
      </figure>
      {Array.isArray(folds) && selected && <ComparisonPlot name={`${selected.algorithm}_cv_${comparison.primary_metric}`} points={folds.map((value, index) => [index + 1, value])} />}
      {selected?.plots && Object.entries(selected.plots).map(([name, points]) => <ComparisonPlot key={name} name={`${selected.algorithm}_${name}`} points={points} />)}
    </div>
  );
}

function QsarComparisonResultsPanel({ job, onBack }: Props) {
  const { t } = useTranslation();
  const [comparison, setComparison] = useState<QsarComparisonResult | null>(null);
  const [sort, setSort] = useState<ComparisonSort>("metric");
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let active = true;
    setComparison(null);
    setError(null);
    getQsarComparison(job.id)
      .then((result) => { if (active) setComparison(result); })
      .catch((loadError) => { if (active) setError(message(loadError, t("qsar.comparisonLoadError"))); });
    return () => { active = false; };
  }, [job.id, t]);

  const rows = comparison ? sortComparisonResults(comparison, sort) : [];
  return (
    <div className="results-panel qsar-comparison-results">
      <div className="section-heading result-heading">
        <div>
          <p className="section-kicker">QSAR MODEL COMPARISON</p>
          <h2>{job.comment || t("qsar.results")}</h2>
          <p>{t("qsar.comparisonAnalysisHint")}</p>
        </div>
        <button className="button quiet" type="button" onClick={onBack}>{t("results.back")}</button>
      </div>

      {!comparison && !error && <div className="empty-state">{t("qsar.comparisonLoading")}</div>}
      {error && <p className="notice error" role="alert">{error}</p>}
      {comparison && <>
        <div className="result-summary-strip">
          <div><strong>{comparison.recommended_algorithm || "—"}</strong><span>{t("qsar.recommendedModel")}</span></div>
          <div><strong>{comparison.primary_metric}</strong><span>{t("qsar.primary")}</span></div>
          <div><strong>{comparison.successful_count}</strong><span>{t("qsar.successfulModels")}</span></div>
          <div><strong>{comparison.failed_count}</strong><span>{t("qsar.failedModels")}</span></div>
        </div>
        <p className="notice success">{t("qsar.recommended", { algorithm: comparison.recommended_algorithm, metric: comparison.primary_metric })}</p>
        {comparison.selection_basis && <p className="muted">{t("qsar.selectionBasis")}: {comparison.selection_basis}</p>}
        <section className="model-management comparison-analysis-section">
          <div className="section-heading compact">
            <div><h3>{t("qsar.results")}</h3></div>
            <label className="comparison-sort">{t("qsar.sortOrder")}
              <select value={sort} onChange={(event) => setSort(event.target.value as ComparisonSort)}>
                <option value="metric">{t("qsar.primary")}</option>
                <option value="algorithm">{t("qsar.modelName")}</option>
              </select>
            </label>
          </div>
          <div className="qsar-preview comparison-result-table-wrap">
            <table>
              <thead><tr><th>{t("qsar.model")}</th><th>{t("qsar.status")}</th><th>{t("qsar.cvMean")}</th><th>Test</th><th>{t("qsar.testR2")}</th><th>{t("qsar.gap")}</th><th>{t("qsar.inAd")}</th></tr></thead>
              <tbody>{rows.map((result) => {
                const cv = finiteNumber(result.metrics?.cross_validation?.[comparison.primary_metric]);
                const std = finiteNumber(result.metrics?.cross_validation?.[`${comparison.primary_metric}_std`]);
                const test = finiteNumber(result.metrics?.test?.[comparison.primary_metric]);
                const testR2 = finiteNumber(result.metrics?.test?.r2);
                const gap = finiteNumber(result.metrics?.performance_gap?.[comparison.primary_metric]);
                const ad = finiteNumber(result.metrics?.applicability_domain_fraction);
                return <tr key={result.algorithm}><td>{result.algorithm}{result.algorithm === comparison.recommended_algorithm ? " ★" : ""}</td><td>{result.status === "completed" ? t("qsar.success") : t("qsar.failure", { error: result.error ?? "" })}</td><td>{cv === null ? "—" : `${cv.toFixed(4)} ± ${std === null ? "—" : std.toFixed(4)}`}</td><td>{test === null ? "—" : test.toFixed(4)}</td><td>{testR2 === null ? "—" : testR2.toFixed(4)}</td><td>{gap === null ? "—" : gap.toFixed(4)}</td><td>{ad === null ? "—" : `${(ad * 100).toFixed(1)}%`}</td></tr>;
              })}</tbody>
            </table>
          </div>
          <ComparisonCharts comparison={comparison} />
        </section>
      </>}
    </div>
  );
}

function GenerationResultsPanel({ job, onBack }: Props) {
  const { t } = useTranslation();
  const [schema, setSchema] = useState<ResultSchemaResponse | null>(null);
  const [molecules, setMolecules] = useState<MoleculeListResponse | null>(null);
  const [statistics, setStatistics] = useState<StatisticsResponse | null>(null);
  const [draftBounds, setDraftBounds] = useState<Bounds>(() => storedBounds(job.id));
  const [appliedBounds, setAppliedBounds] = useState<Bounds>(() => storedBounds(job.id));
  const [sortBy, setSortBy] = useState("MW");
  const [sortDirection, setSortDirection] = useState("asc");
  const [scatterX, setScatterX] = useState("MW");
  const [scatterY, setScatterY] = useState("LogP");
  const [offset, setOffset] = useState(0);
  const [display, setDisplay] = useState<"cards" | "table">("table");
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const baseQuery = useMemo<ResultQuery>(() => ({
    ...filtersToQuery(appliedBounds),
    sort_by: sortBy,
    sort_direction: sortDirection
  }), [appliedBounds, sortBy, sortDirection]);

  useEffect(() => {
    const controller = new AbortController();
    setLoading(true);
    Promise.all([
      getResultSchema(job.id, controller.signal),
      getResultMolecules(job.id, { ...baseQuery, limit: PAGE_SIZE, offset }, controller.signal),
      getResultStatistics(job.id, { ...baseQuery, scatter_x: scatterX, scatter_y: scatterY }, controller.signal)
    ])
      .then(([nextSchema, nextMolecules, nextStatistics]) => {
        setSchema(nextSchema);
        setMolecules(nextMolecules);
        setStatistics(nextStatistics);
        setError(null);
      })
      .catch((loadError) => {
        if (!controller.signal.aborted) setError(message(loadError, t("results.loadError")));
      })
      .finally(() => {
        if (!controller.signal.aborted) setLoading(false);
      });
    return () => controller.abort();
  }, [baseQuery, job.id, offset, scatterX, scatterY, t]);

  const availableDescriptors = schema?.descriptors.map((item) => item.name) ?? [...FILTER_DESCRIPTORS];
  const descriptorSchemas = useMemo(
    () => new Map(schema?.descriptors.map((item) => [item.name, item]) ?? []),
    [schema]
  );
  const totalPages = Math.max(1, Math.ceil((molecules?.total ?? 0) / PAGE_SIZE));
  const currentPage = Math.floor(offset / PAGE_SIZE) + 1;

  const applyFilters = () => {
    const next = structuredClone(draftBounds);
    window.sessionStorage.setItem(`reinvent4.results.${job.id}`, JSON.stringify(next));
    setAppliedBounds(next);
    setOffset(0);
  };

  const clearFilters = () => {
    const empty = blankBounds();
    setDraftBounds(empty);
    setAppliedBounds(blankBounds());
    window.sessionStorage.removeItem(`reinvent4.results.${job.id}`);
    setOffset(0);
  };

  const sortFromTable = (descriptor: string) => {
    if (sortBy === descriptor) {
      setSortDirection((current) => current === "asc" ? "desc" : "asc");
    } else {
      setSortBy(descriptor);
      setSortDirection("asc");
    }
    setOffset(0);
  };

  return (
    <div className="results-panel">
      <div className="section-heading result-heading">
        <div>
          <p className="section-kicker">RESULT EXPLORER</p>
          <h2>{job.comment || t("results.generated")}</h2>
          <p>{job.parent_smiles_canonical || job.parent_smiles}</p>
        </div>
        <button className="button quiet" type="button" onClick={onBack}>{t("results.back")}</button>
      </div>

      <div className="result-summary-strip">
        <div><strong>{schema?.valid_rows ?? "—"}</strong><span>{t("results.valid")}</span></div>
        <div><strong>{molecules?.total ?? "—"}</strong><span>{t("results.filtered")}</span></div>
        <div><strong>{schema?.invalid_smiles ?? "—"}</strong><span>{t("results.invalid")}</span></div>
        <div><strong>{schema ? `RDKit ${schema.rdkit_version}` : "—"}</strong><span>{t("results.environment")}</span></div>
      </div>

      {schema?.warnings.map((warning) => <p className="notice warning" key={warning}>{warning}</p>)}
      {error && <p className="notice error" role="alert">{error}</p>}

      <section className="filter-panel" aria-label={t("results.filtersLabel")}>
        <div className="filter-title"><h3>{t("results.filters")}</h3><span>{t("results.inclusive")}</span></div>
        <p className="filter-group-label">{t("results.basic")}</p>
        <div className="filter-grid">
          {CORE_DESCRIPTORS.map((descriptor) => (
            <DescriptorFilter
              key={descriptor}
              descriptor={descriptor}
              bounds={draftBounds[descriptor]}
              limits={descriptorSchemas.get(descriptor)}
              onChange={(next) => setDraftBounds((current) => ({ ...current, [descriptor]: next }))}
            />
          ))}
        </div>
        <details className="advanced-filters">
          <summary>{t("results.additional")}</summary>
          <div className="filter-grid">
            {ADDED_DESCRIPTORS.map((descriptor) => (
              <DescriptorFilter
                key={descriptor}
                descriptor={descriptor}
                bounds={draftBounds[descriptor]}
                limits={descriptorSchemas.get(descriptor)}
                onChange={(next) => setDraftBounds((current) => ({ ...current, [descriptor]: next }))}
              />
            ))}
          </div>
        </details>
        <div className="filter-actions">
          <label>{t("results.sort")}
            <select value={sortBy} onChange={(event) => { setSortBy(event.target.value); setOffset(0); }}>
              {availableDescriptors.map((name) => <option value={name} key={name}>{name}</option>)}
            </select>
          </label>
          <label>{t("results.order")}
            <select value={sortDirection} onChange={(event) => { setSortDirection(event.target.value); setOffset(0); }}>
              <option value="asc">{t("results.ascending")}</option>
              <option value="desc">{t("results.descending")}</option>
            </select>
          </label>
          <button className="button quiet" type="button" onClick={clearFilters}>{t("results.clear")}</button>
          <button className="button primary" type="button" onClick={applyFilters}>{t("results.apply")}</button>
        </div>
      </section>

      <div className="result-toolbar">
        <div className="segmented" aria-label={t("results.viewMode")}>
          <button className={display === "cards" ? "active" : ""} type="button" onClick={() => setDisplay("cards")}>{t("results.cards")}</button>
          <button className={display === "table" ? "active" : ""} type="button" onClick={() => setDisplay("table")}>{t("results.table")}</button>
        </div>
        <div className="download-actions">
          <a className="button secondary" href={resultDownloadUrl(job.id, "csv", baseQuery)}>CSV</a>
          <a className="button secondary" href={resultDownloadUrl(job.id, "sdf", baseQuery)}>SDF</a>
        </div>
      </div>

      {loading && !molecules ? <div className="empty-state">{t("results.loading")}</div> : null}
      {molecules && molecules.items.length === 0 ? <div className="empty-state">{t("results.empty")}</div> : null}
      {molecules && molecules.items.length > 0 && display === "cards" ? (
        <div className="molecule-cards">
          {molecules.items.map((item) => (
            <article className="molecule-card" key={item.row_number}>
              <img src={`data:image/svg+xml;charset=utf-8,${encodeURIComponent(item.structure_svg)}`} alt={t("results.structureAlt", { row: item.row_number })} />
              <div className="molecule-meta">
                <strong>#{item.row_number}</strong>
                <code title={item.canonical_smiles}>{item.canonical_smiles}</code>
                <dl>{FILTER_DESCRIPTORS.map((name) => <div key={name}><dt>{name}</dt><dd>{formatNumber(item.values[name])}</dd></div>)}</dl>
              </div>
            </article>
          ))}
        </div>
      ) : null}

      {molecules && molecules.items.length > 0 && display === "table" ? (
        <div className="result-table-wrap">
          <table className="result-table">
            <thead><tr><th>#</th><th>Structure</th>{FILTER_DESCRIPTORS.map((name) => (
              <th
                key={name}
                aria-sort={sortBy === name ? (sortDirection === "asc" ? "ascending" : "descending") : "none"}
              >
                <button className="sort-header" type="button" aria-label={t("results.sortBy", { name })} onClick={() => sortFromTable(name)}>
                  <span>{name}</span><span aria-hidden="true">{sortBy === name ? (sortDirection === "asc" ? "↑" : "↓") : "↕"}</span>
                </button>
              </th>
            ))}</tr></thead>
            <tbody>{molecules.items.map((item) => (
              <tr key={item.row_number}>
                <td>{item.row_number}</td>
                <td><img src={`data:image/svg+xml;charset=utf-8,${encodeURIComponent(item.structure_svg)}`} alt="" /></td>
                {FILTER_DESCRIPTORS.map((name) => <td key={name}>{formatNumber(item.values[name])}</td>)}
              </tr>
            ))}</tbody>
          </table>
        </div>
      ) : null}

      {molecules && molecules.total > 0 ? (
        <nav className="pagination" aria-label={t("results.pages")}>
          <button className="button quiet" type="button" disabled={offset === 0} onClick={() => setOffset(Math.max(0, offset - PAGE_SIZE))}>{t("results.previous")}</button>
          <span>{currentPage} / {totalPages}</span>
          <button className="button quiet" type="button" disabled={currentPage >= totalPages} onClick={() => setOffset(offset + PAGE_SIZE)}>{t("results.next")}</button>
        </nav>
      ) : null}

      {statistics ? (
        <section className="analytics-grid">
          <div className="analytics-card scatter-card">
            <div className="analytics-title">
              <h3>{t("results.relationship")}</h3>
              <div>
                <select aria-label={t("results.xAxis")} value={scatterX} onChange={(event) => setScatterX(event.target.value)}>{availableDescriptors.map((name) => <option key={name}>{name}</option>)}</select>
                <span>×</span>
                <select aria-label={t("results.yAxis")} value={scatterY} onChange={(event) => setScatterY(event.target.value)}>{availableDescriptors.map((name) => <option key={name}>{name}</option>)}</select>
              </div>
            </div>
            <ScatterPlot statistics={statistics} />
            {statistics.scatter_sampled && <p className="muted">{t("results.sampled")}</p>}
          </div>
          <div className="analytics-card">
            <h3>{t("results.distribution")}</h3>
            <div className="histogram-list">
              {FILTER_DESCRIPTORS.map((name) => {
                const bins = statistics.histograms[name] ?? [];
                const peak = Math.max(1, ...bins.map((bin) => bin.count));
                return (
                  <div className="mini-histogram" key={name}>
                    <strong>{name}</strong>
                    <div className="histogram-bars">{bins.map((bin, index) => <span key={`${bin.start}-${index}`} style={{ height: `${Math.max(3, (bin.count / peak) * 48)}px` }} title={`${formatNumber(bin.start)}–${formatNumber(bin.end)}: ${bin.count}`} />)}</div>
                    <div className="histogram-axis" aria-label={t("results.histogramAxis", { name })}><span>{formatNumber(bins.at(0)?.start)}</span><span>{formatNumber(bins.at(-1)?.end)}</span></div>
                    <small>{t("results.median")} {formatNumber(statistics.summaries[name]?.median)}</small>
                  </div>
                );
              })}
            </div>
          </div>
        </section>
      ) : null}
    </div>
  );
}

export function ResultsPanel(props: Props) {
  if (props.job.job_type === "qsar_comparison") {
    return <QsarComparisonResultsPanel {...props} />;
  }
  if (props.job.job_type === "transfer_learning" || props.job.job_type === "reinforcement_learning") {
    return <LearningResultsPanel {...props} />;
  }
  return <GenerationResultsPanel {...props} />;
}
