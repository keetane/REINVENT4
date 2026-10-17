import { useEffect, useMemo, useState } from "react";

import {
  ApiError,
  getResultMolecules,
  getResultSchema,
  getResultStatistics,
  resultDownloadUrl,
  type DescriptorSchema,
  type Job,
  type MoleculeListResponse,
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

function message(error: unknown): string {
  return error instanceof ApiError ? error.message : "結果データを読み込めませんでした。";
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
          aria-label={`${descriptor} 最小`}
          inputMode="decimal"
          placeholder="min"
          value={bounds.minimum}
          onChange={(event) => onChange({ ...bounds, minimum: event.target.value })}
        />
        <span>–</span>
        <input
          aria-label={`${descriptor} 最大`}
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
            aria-label={`${descriptor} 下限スライダー`}
            min={minimumValue}
            max={maximumValue}
            step={step}
            value={safeLower}
            onChange={(event) => onChange({ ...bounds, minimum: sliderValue(Math.min(Number(event.target.value), safeUpper)) })}
          />
          <input
            type="range"
            aria-label={`${descriptor} 上限スライダー`}
            min={minimumValue}
            max={maximumValue}
            step={step}
            value={safeUpper}
            onChange={(event) => onChange({ ...bounds, maximum: sliderValue(Math.max(Number(event.target.value), safeLower)) })}
          />
        </div>
      ) : (
        <div className="range-static">{rangeMinimum === null || rangeMinimum === undefined ? "範囲を計算中" : formatNumber(rangeMinimum)}</div>
      )}
    </fieldset>
  );
}

function ScatterPlot({ statistics }: { statistics: StatisticsResponse }) {
  const width = 560;
  const height = 290;
  const left = 58;
  const right = 18;
  const top = 18;
  const bottom = 46;
  const points = statistics.scatter;
  if (!points.length) return <div className="chart-empty">表示できる点がありません。</div>;
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
    <svg className="scatter-plot" viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${statistics.scatter_x} と ${statistics.scatter_y} の散布図`}>
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

export function ResultsPanel({ job, onBack }: Props) {
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
        if (!controller.signal.aborted) setError(message(loadError));
      })
      .finally(() => {
        if (!controller.signal.aborted) setLoading(false);
      });
    return () => controller.abort();
  }, [baseQuery, job.id, offset, scatterX, scatterY]);

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
          <h2>{job.comment || "生成結果"}</h2>
          <p>{job.parent_smiles_canonical || job.parent_smiles}</p>
        </div>
        <button className="button quiet" type="button" onClick={onBack}>履歴へ戻る</button>
      </div>

      <div className="result-summary-strip">
        <div><strong>{schema?.valid_rows ?? "—"}</strong><span>有効分子</span></div>
        <div><strong>{molecules?.total ?? "—"}</strong><span>絞り込み後</span></div>
        <div><strong>{schema?.invalid_smiles ?? "—"}</strong><span>無効 SMILES</span></div>
        <div><strong>{schema ? `RDKit ${schema.rdkit_version}` : "—"}</strong><span>計算環境</span></div>
      </div>

      {schema?.warnings.map((warning) => <p className="notice warning" key={warning}>{warning}</p>)}
      {error && <p className="notice error" role="alert">{error}</p>}

      <section className="filter-panel" aria-label="結果フィルター">
        <div className="filter-title"><h3>物性フィルター</h3><span>範囲は両端を含みます</span></div>
        <p className="filter-group-label">基本物性</p>
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
          <summary>追加物性・類似度</summary>
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
          <label>並び替え
            <select value={sortBy} onChange={(event) => { setSortBy(event.target.value); setOffset(0); }}>
              {availableDescriptors.map((name) => <option value={name} key={name}>{name}</option>)}
            </select>
          </label>
          <label>順序
            <select value={sortDirection} onChange={(event) => { setSortDirection(event.target.value); setOffset(0); }}>
              <option value="asc">昇順</option>
              <option value="desc">降順</option>
            </select>
          </label>
          <button className="button quiet" type="button" onClick={clearFilters}>クリア</button>
          <button className="button primary" type="button" onClick={applyFilters}>絞り込みを適用</button>
        </div>
      </section>

      <div className="result-toolbar">
        <div className="segmented" aria-label="表示形式">
          <button className={display === "cards" ? "active" : ""} type="button" onClick={() => setDisplay("cards")}>カード</button>
          <button className={display === "table" ? "active" : ""} type="button" onClick={() => setDisplay("table")}>表</button>
        </div>
        <div className="download-actions">
          <a className="button secondary" href={resultDownloadUrl(job.id, "csv", baseQuery)}>CSV</a>
          <a className="button secondary" href={resultDownloadUrl(job.id, "sdf", baseQuery)}>SDF</a>
        </div>
      </div>

      {loading && !molecules ? <div className="empty-state">結果を読み込み中…</div> : null}
      {molecules && molecules.items.length === 0 ? <div className="empty-state">条件に一致する分子はありません。</div> : null}
      {molecules && molecules.items.length > 0 && display === "cards" ? (
        <div className="molecule-cards">
          {molecules.items.map((item) => (
            <article className="molecule-card" key={item.row_number}>
              <img src={`data:image/svg+xml;charset=utf-8,${encodeURIComponent(item.structure_svg)}`} alt={`分子 ${item.row_number} の構造`} />
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
                <button className="sort-header" type="button" aria-label={`${name}で並び替え`} onClick={() => sortFromTable(name)}>
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
        <nav className="pagination" aria-label="結果ページ">
          <button className="button quiet" type="button" disabled={offset === 0} onClick={() => setOffset(Math.max(0, offset - PAGE_SIZE))}>前へ</button>
          <span>{currentPage} / {totalPages}</span>
          <button className="button quiet" type="button" disabled={currentPage >= totalPages} onClick={() => setOffset(offset + PAGE_SIZE)}>次へ</button>
        </nav>
      ) : null}

      {statistics ? (
        <section className="analytics-grid">
          <div className="analytics-card scatter-card">
            <div className="analytics-title">
              <h3>物性の関係</h3>
              <div>
                <select aria-label="散布図 X 軸" value={scatterX} onChange={(event) => setScatterX(event.target.value)}>{availableDescriptors.map((name) => <option key={name}>{name}</option>)}</select>
                <span>×</span>
                <select aria-label="散布図 Y 軸" value={scatterY} onChange={(event) => setScatterY(event.target.value)}>{availableDescriptors.map((name) => <option key={name}>{name}</option>)}</select>
              </div>
            </div>
            <ScatterPlot statistics={statistics} />
            {statistics.scatter_sampled && <p className="muted">再現可能な等間隔サンプル 500 点を表示しています。</p>}
          </div>
          <div className="analytics-card">
            <h3>分布</h3>
            <div className="histogram-list">
              {FILTER_DESCRIPTORS.map((name) => {
                const bins = statistics.histograms[name] ?? [];
                const peak = Math.max(1, ...bins.map((bin) => bin.count));
                return (
                  <div className="mini-histogram" key={name}>
                    <strong>{name}</strong>
                    <div className="histogram-bars">{bins.map((bin, index) => <span key={`${bin.start}-${index}`} style={{ height: `${Math.max(3, (bin.count / peak) * 48)}px` }} title={`${formatNumber(bin.start)}–${formatNumber(bin.end)}: ${bin.count}`} />)}</div>
                    <div className="histogram-axis" aria-label={`${name} 横軸`}><span>{formatNumber(bins.at(0)?.start)}</span><span>{formatNumber(bins.at(-1)?.end)}</span></div>
                    <small>中央値 {formatNumber(statistics.summaries[name]?.median)}</small>
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
