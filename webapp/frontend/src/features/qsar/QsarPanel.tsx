import { useEffect, useMemo, useState } from "react";
import { useTranslation } from "react-i18next";

import {
  ApiError,
  createQsarComparisonJob,
  createQsarJob,
  deleteQsarModel,
  getQsarCapabilities,
  getQsarComparison,
  getQsarEvaluation,
  inspectQsarData,
  listQsarModels,
  predictQsar,
  qsarArtifactUrl,
  registerQsarComparisonModel,
  validateParent,
  type Job,
  type ParentValidationResponse,
  type QsarCapabilities,
  type QsarComparisonResult,
  type QsarDataInspection,
  type QsarDataPayload,
  type QsarModel,
  type QsarPredictionRow
} from "../../api/client";
import {
  isDataUploadTooLarge,
  suggestQsarTargetColumn,
  suggestSmilesColumn
} from "../../lib/uploads";
import { StandardizationNotice } from "../standardization/StandardizationNotice";
import { ParentInput } from "../parent-input/ParentInput";

interface Props {
  onCreated: (job: Job) => void;
}

type HyperParameters = Record<string, string | number>;

function message(error: unknown, fallback: string): string {
  return error instanceof ApiError ? error.message : fallback;
}

function metric(model: QsarModel, key: string): string {
  const test = model.metrics.test as Record<string, unknown> | undefined;
  const value = test?.[key];
  return typeof value === "number" ? value.toFixed(3) : "—";
}

export function QsarPanel({ onCreated }: Props) {
  const { t } = useTranslation();
  const [fileName, setFileName] = useState("");
  const [content, setContent] = useState("");
  const [inspection, setInspection] = useState<QsarDataInspection | null>(null);
  const [capabilities, setCapabilities] = useState<QsarCapabilities | null>(null);
  const [models, setModels] = useState<QsarModel[]>([]);
  const [smilesColumn, setSmilesColumn] = useState("");
  const [targetColumn, setTargetColumn] = useState("");
  const [taskType, setTaskType] = useState<"regression" | "classification">("regression");
  const [targetTransform, setTargetTransform] = useState<"none" | "p_activity">("none");
  const [concentrationUnit, setConcentrationUnit] = useState<QsarDataPayload["concentration_unit"]>("nM");
  const [activityType, setActivityType] = useState<QsarDataPayload["activity_type"]>("IC50");
  const [duplicateStrategy, setDuplicateStrategy] = useState<QsarDataPayload["duplicate_strategy"]>("exclude_conflicts");
  const [algorithm, setAlgorithm] = useState("random_forest");
  const [mode, setMode] = useState<"single" | "comparison">("single");
  const [selectedAlgorithms, setSelectedAlgorithms] = useState<string[]>(["random_forest", "extra_trees", "gradient_boosting"]);
  const [primaryMetric, setPrimaryMetric] = useState("rmse");
  const [comparisonJobId, setComparisonJobId] = useState<string | null>(null);
  const [comparison, setComparison] = useState<QsarComparisonResult | null>(null);
  const [comparisonSort, setComparisonSort] = useState<"metric" | "algorithm">("metric");
  const [fingerprint, setFingerprint] = useState("morgan");
  const [featureMode, setFeatureMode] = useState("fingerprint");
  const [splitMethod, setSplitMethod] = useState("scaffold");
  const [trainFraction, setTrainFraction] = useState(0.8);
  const [folds, setFolds] = useState(3);
  const [adThreshold, setAdThreshold] = useState(0.4);
  const [comment, setComment] = useState("");
  const [advanced, setAdvanced] = useState(false);
  const [hyperparameters, setHyperparameters] = useState<HyperParameters>({});
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [evaluation, setEvaluation] = useState<{ modelId: string; plots: Record<string, number[][]> } | null>(null);
  const [predictionModelId, setPredictionModelId] = useState<string | null>(null);
  const [predictionSmiles, setPredictionSmiles] = useState("");
  const [predictionValidation, setPredictionValidation] = useState<ParentValidationResponse | null>(null);
  const [predictionResult, setPredictionResult] = useState<QsarPredictionRow | null>(null);
  const [predictionBusy, setPredictionBusy] = useState(false);
  const [predictionError, setPredictionError] = useState<string | null>(null);

  const dataPayload = useMemo<QsarDataPayload>(() => ({
    file_name: fileName,
    content,
    smiles_column: smilesColumn || null,
    target_column: targetColumn || null,
    task_type: taskType,
    target_transform: taskType === "classification" ? "none" : targetTransform,
    concentration_unit: concentrationUnit,
    activity_type: activityType,
    duplicate_strategy: taskType === "classification" && duplicateStrategy === "mean" ? "exclude_conflicts" : duplicateStrategy
  }), [fileName, content, smilesColumn, targetColumn, taskType, targetTransform, concentrationUnit, activityType, duplicateStrategy]);

  const refreshModels = async () => {
    const response = await listQsarModels();
    setModels(response.items);
  };

  const openPrediction = (modelId: string) => {
    setPredictionModelId((current) => current === modelId ? null : modelId);
    setPredictionSmiles("");
    setPredictionValidation(null);
    setPredictionResult(null);
    setPredictionError(null);
  };

  const validatePredictionSmiles = async () => {
    if (!predictionSmiles.trim()) return;
    setPredictionBusy(true);
    setPredictionError(null);
    try {
      const checked = await validateParent(predictionSmiles, "reinvent");
      setPredictionSmiles(checked.canonical_smiles);
      setPredictionValidation(checked);
      setPredictionResult(null);
    } catch (reason) {
      setPredictionValidation(null);
      setPredictionResult(null);
      setPredictionError(message(reason, t("qsar.predictionError")));
    } finally {
      setPredictionBusy(false);
    }
  };

  const calculatePrediction = async () => {
    if (!predictionModelId || !predictionSmiles.trim()) return;
    setPredictionBusy(true);
    setPredictionError(null);
    try {
      const checked = await validateParent(predictionSmiles, "reinvent");
      setPredictionSmiles(checked.canonical_smiles);
      setPredictionValidation(checked);
      const response = await predictQsar(predictionModelId, [checked.canonical_smiles]);
      setPredictionResult(response.predictions[0] ?? null);
      if (!response.predictions[0]) setPredictionError(t("qsar.predictionError"));
    } catch (reason) {
      setPredictionResult(null);
      setPredictionError(message(reason, t("qsar.predictionError")));
    } finally {
      setPredictionBusy(false);
    }
  };

  useEffect(() => {
    void Promise.all([getQsarCapabilities().then(setCapabilities), refreshModels()]).catch((reason) => setError(message(reason, t("qsar.error"))));
  }, [t]);

  useEffect(() => {
    const selected = capabilities?.algorithms.find((item) => item.key === algorithm);
    if (selected && !selected.task_types.includes(taskType)) {
      const fallback = capabilities?.algorithms.find((item) => item.available && item.task_types.includes(taskType));
      if (fallback) {
        setAlgorithm(fallback.key);
        setHyperparameters({});
      }
    }
  }, [algorithm, capabilities, taskType]);

  useEffect(() => {
    setPrimaryMetric(taskType === "regression" ? "rmse" : "roc_auc");
    setSelectedAlgorithms((current) => current.filter((key) => capabilities?.algorithms.find((item) => item.key === key)?.task_types.includes(taskType)));
  }, [capabilities, taskType]);

  useEffect(() => {
    if (!comparisonJobId || comparison) return;
    const poll = () => void getQsarComparison(comparisonJobId).then(setComparison).catch((reason) => {
      if (!(reason instanceof ApiError) || reason.status !== 409) setError(message(reason, t("qsar.error")));
    });
    poll();
    const timer = window.setInterval(poll, 3000);
    return () => window.clearInterval(timer);
  }, [comparison, comparisonJobId, t]);

  useEffect(() => {
    if (!content) return;
    const timer = window.setTimeout(() => {
      void inspectQsarData(dataPayload).then(setInspection).catch((reason) => setError(message(reason, t("qsar.error"))));
    }, 250);
    return () => window.clearTimeout(timer);
  }, [dataPayload, content, t]);

  const onFile = async (file: File | undefined) => {
    if (!file) return;
    setError(null);
    if (isDataUploadTooLarge(file)) {
      setFileName("");
      setContent("");
      setSmilesColumn("");
      setTargetColumn("");
      setInspection(null);
      setError(t("upload.tooLarge"));
      return;
    }
    const text = await file.text();
    setFileName(file.name);
    setContent(text);
    const initial = await inspectQsarData({ ...dataPayload, file_name: file.name, content: text, smiles_column: null, target_column: null });
    setInspection(initial);
    const suggestedSmiles = suggestSmilesColumn(initial.columns);
    const suggestedTarget = suggestQsarTargetColumn(initial.columns, suggestedSmiles);
    setSmilesColumn(suggestedSmiles);
    setTargetColumn(suggestedTarget);
  };

  const start = async () => {
    setBusy(true);
    setError(null);
    try {
      const submittedHyperparameters = algorithm === "mlp" && typeof hyperparameters.hidden_layer_sizes === "string"
        ? { ...hyperparameters, hidden_layer_sizes: hyperparameters.hidden_layer_sizes.split(",").map((item) => Number(item.trim())).filter((item) => Number.isInteger(item) && item > 0) }
        : hyperparameters;
      const common = {
        ...dataPayload,
        algorithm,
        comment,
        fingerprint,
        feature_mode: featureMode,
        morgan_radius: 2,
        bit_length: 2048,
        use_chirality: true,
        split_method: splitMethod,
        train_fraction: trainFraction,
        cross_validation_folds: folds,
        ad_threshold: adThreshold,
        device: "cpu",
        advanced,
        hyperparameters: advanced ? submittedHyperparameters : {}
      };
      const job = mode === "single"
        ? await createQsarJob(common)
        : await createQsarComparisonJob({ ...common, algorithms: selectedAlgorithms, primary_metric: primaryMetric });
      if (mode === "comparison") {
        setComparisonJobId(job.id);
        setComparison(null);
      }
      onCreated(job);
    } catch (reason) {
      setError(message(reason, t("qsar.error")));
    } finally {
      setBusy(false);
    }
  };

  const downloadRejected = () => {
    if (!inspection) return;
    const rows = ["Row,SMILES,Target,Reason", ...inspection.rejected_rows.map((row) =>
      [row.row_number, row.smiles, row.target, row.reason].map((value) => `"${String(value).replaceAll('"', '""')}"`).join(",")
    )];
    const url = URL.createObjectURL(new Blob([rows.join("\n")], { type: "text/csv;charset=utf-8" }));
    const link = document.createElement("a");
    link.href = url;
    link.download = "qsar-rejected-rows.csv";
    link.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="qsar-panel">
      <div className="section-heading">
        <div><p className="section-kicker">QSAR MODELING</p><h2>{t("qsar.title")}</h2><p>{t("qsar.intro")}</p></div>
      </div>

      <StandardizationNotice version={inspection?.standardization_version} />

      <div className="qsar-grid">
        <div className="field-group wide">
          <label htmlFor="qsar-file">{t("qsar.data")}</label>
          <input id="qsar-file" type="file" accept=".csv,text/csv" onChange={(event) => void onFile(event.target.files?.[0])} />
          <p className="field-hint">{t("qsar.dataHint")}</p>
        </div>
        <div className="field-group"><label htmlFor="qsar-smiles">{t("qsar.smilesColumn")}</label><select id="qsar-smiles" value={smilesColumn} onChange={(event) => setSmilesColumn(event.target.value)}>{inspection?.columns.map((column) => <option key={column}>{column}</option>)}</select></div>
        <div className="field-group"><label htmlFor="qsar-target">{t("qsar.targetColumn")}</label><select id="qsar-target" value={targetColumn} onChange={(event) => setTargetColumn(event.target.value)}>{inspection?.columns.map((column) => <option key={column}>{column}</option>)}</select></div>
        <div className="field-group"><label htmlFor="qsar-task">{t("qsar.task")}</label><select id="qsar-task" value={taskType} onChange={(event) => setTaskType(event.target.value as typeof taskType)}><option value="regression">{t("qsar.regression")}</option><option value="classification">{t("qsar.classification")}</option></select><p className="field-hint">{t("qsar.taskHint")}</p></div>
        <div className="field-group"><label htmlFor="qsar-duplicate">{t("qsar.duplicates")}</label><select id="qsar-duplicate" value={duplicateStrategy} onChange={(event) => setDuplicateStrategy(event.target.value as typeof duplicateStrategy)}><option value="exclude_conflicts">{t("qsar.excludeConflicts")}</option><option value="keep_first">{t("qsar.keepFirst")}</option>{taskType === "regression" && <option value="mean">{t("qsar.useMean")}</option>}</select></div>
        {taskType === "regression" && <><div className="field-group"><label htmlFor="qsar-transform">{t("qsar.targetTransform")}</label><select id="qsar-transform" value={targetTransform} onChange={(event) => setTargetTransform(event.target.value as typeof targetTransform)}><option value="none">{t("qsar.noTransform")}</option><option value="p_activity">{t("qsar.toPActivity")}</option></select></div>{targetTransform === "p_activity" && <><div className="field-group"><label htmlFor="qsar-activity">{t("qsar.activityType")}</label><select id="qsar-activity" value={activityType} onChange={(event) => setActivityType(event.target.value as typeof activityType)}>{["IC50", "EC50", "Ki", "Kd", "custom"].map((value) => <option key={value}>{value}</option>)}</select></div><div className="field-group"><label htmlFor="qsar-unit">{t("qsar.unit")}</label><select id="qsar-unit" value={concentrationUnit} onChange={(event) => setConcentrationUnit(event.target.value as typeof concentrationUnit)}>{["M", "mM", "uM", "μM", "nM", "pM"].map((unit) => <option key={unit}>{unit}</option>)}</select></div></>}</>}
      </div>

      {inspection && <div className="qsar-summary">
        <div><span>{t("qsar.allRows")}</span><strong>{inspection.total_rows}</strong></div><div><span>標準化済み</span><strong>{inspection.standardized_count}</strong></div><div><span>{t("qsar.valid")}</span><strong>{inspection.final_accepted_count}</strong></div><div><span>{t("qsar.invalidSmiles")}</span><strong>{inspection.invalid_smiles_count}</strong></div><div><span>{t("qsar.targetErrors")}</span><strong>{inspection.missing_target_count + inspection.invalid_target_count}</strong></div><div><span>{t("qsar.duplicateCount")}</span><strong>{inspection.duplicate_count}</strong></div><div><span>{t("qsar.conflicts")}</span><strong>{inspection.conflicting_duplicate_count}</strong></div><div><span>モデル非互換</span><strong>0</strong></div>
        <div><span>フラグメント抽出</span><strong>{inspection.fragment_extracted_count}</strong></div><div><span>中和</span><strong>{inspection.neutralized_count}</strong></div><div><span>永久電荷保持</span><strong>{inspection.charged_structure_retained_count}</strong></div><div><span>MW≤200除外</span><strong>{inspection.low_molecular_weight_excluded_count}</strong></div>
        {inspection.rejected_rows.length > 0 && <button type="button" className="button quiet" onClick={downloadRejected}>{t("qsar.rejectedCsv")}</button>}
      </div>}
      {inspection && inspection.rejected_rows.length > 0 && <details className="advanced-settings"><summary>{t("qsar.rejectedList", { count: inspection.rejected_rows.length })}</summary><div className="qsar-preview"><table><thead><tr><th>{t("qsar.row")}</th><th>SMILES</th><th>{t("qsar.target")}</th><th>{t("qsar.reason")}</th></tr></thead><tbody>{inspection.rejected_rows.slice(0, 200).map((row) => <tr key={`${row.row_number}-${row.reason}`}><td>{row.row_number}</td><td>{row.smiles}</td><td>{row.target}</td><td>{row.reason}</td></tr>)}</tbody></table></div></details>}

      {inspection?.preview.length ? <div className="qsar-preview"><table><thead><tr>{inspection.columns.map((column) => <th key={column}>{column}</th>)}</tr></thead><tbody>{inspection.preview.map((row, index) => <tr key={index}>{inspection.columns.map((column) => <td key={column}>{row[column]}</td>)}</tr>)}</tbody></table></div> : null}

      <div className="divider" />
      <div className="qsar-grid">
        <div className="field-group wide"><label htmlFor="qsar-mode">{t("qsar.mode")}</label><select id="qsar-mode" value={mode} onChange={(event) => setMode(event.target.value as typeof mode)}><option value="single">Single Model</option><option value="comparison">Model Comparison</option></select><p className="field-hint">{t("qsar.comparisonHint")}</p></div>
        {mode === "comparison" && <div className="field-group wide"><label>{t("qsar.comparisonAlgorithms")}</label><div className="comparison-presets"><button type="button" className="button quiet" onClick={() => setSelectedAlgorithms(capabilities?.algorithms.filter((item) => item.available && item.task_types.includes(taskType) && !["chemprop_v2", "qsartuna"].includes(item.key)).map((item) => item.key) ?? [])}>{t("qsar.selectAll")}</button><button type="button" className="button quiet" onClick={() => setSelectedAlgorithms(["random_forest", "extra_trees", "gradient_boosting", "lightgbm", "xgboost", "catboost"].filter((key) => capabilities?.algorithms.find((item) => item.key === key)?.available))}>{t("qsar.treeModels")}</button><button type="button" className="button quiet" onClick={() => setSelectedAlgorithms(["svm", "kernel_ridge", "gaussian_process"].filter((key) => capabilities?.algorithms.find((item) => item.key === key)?.available && capabilities?.algorithms.find((item) => item.key === key)?.task_types.includes(taskType))) }>{t("qsar.kernelModels")}</button><button type="button" className="button quiet" onClick={() => setSelectedAlgorithms(["mlp"].filter((key) => capabilities?.algorithms.find((item) => item.key === key)?.available))}>{t("qsar.neuralModels")}</button><button type="button" className="button quiet" onClick={() => setSelectedAlgorithms([])}>{t("qsar.clearSelection")}</button></div><div className="algorithm-checks">{capabilities?.algorithms.filter((item) => item.available && item.task_types.includes(taskType) && !["chemprop_v2", "qsartuna"].includes(item.key)).map((item) => <label key={item.key}><input type="checkbox" checked={selectedAlgorithms.includes(item.key)} onChange={(event) => setSelectedAlgorithms((current) => event.target.checked ? [...current, item.key] : current.filter((key) => key !== item.key))} />{item.label}</label>)}</div></div>}
        {mode === "comparison" && <div className="field-group"><label htmlFor="qsar-primary">{t("qsar.primaryMetric")}</label><select id="qsar-primary" value={primaryMetric} onChange={(event) => setPrimaryMetric(event.target.value)}>{taskType === "regression" ? <><option value="rmse">RMSE {t("qsar.minimum")}</option><option value="mae">MAE {t("qsar.minimum")}</option><option value="r2">R² {t("qsar.maximum")}</option></> : <><option value="roc_auc">ROC-AUC {t("qsar.maximum")}</option><option value="pr_auc">PR-AUC {t("qsar.maximum")}</option></>}</select></div>}
        <div className="field-group wide"><label htmlFor="qsar-algorithm">{t("qsar.algorithm")}</label><select id="qsar-algorithm" value={algorithm} onChange={(event) => { setAlgorithm(event.target.value); setHyperparameters({}); }}>{capabilities?.algorithms.map((item) => { const supported = item.task_types.includes(taskType); return <option key={item.key} value={item.key} disabled={!item.available || !supported}>{item.label}{!item.available ? t("common.unavailable") : !supported ? t("qsar.taskOnly", { task: taskType === "regression" ? t("qsar.classification") : t("qsar.regression") }) : ""}</option>; })}</select>{capabilities?.algorithms.filter((item) => !item.available).map((item) => <p className="field-hint" key={item.key}>{item.label}: {item.reason}</p>)}</div>
        <div className="field-group"><label htmlFor="qsar-feature">{t("qsar.feature")}</label><select id="qsar-feature" value={featureMode} disabled={algorithm === "chemprop_v2" || algorithm === "qsartuna"} onChange={(event) => setFeatureMode(event.target.value)}><option value="fingerprint">Fingerprint</option><option value="descriptors">{t("qsar.descriptors")}</option><option value="both">{t("qsar.both")}</option></select>{algorithm === "chemprop_v2" && <p className="field-hint">{t("qsar.chempropHint")}</p>}{algorithm === "qsartuna" && <p className="field-hint">{t("qsar.qsartunaHint")}</p>}</div>
        <div className="field-group"><label htmlFor="qsar-fingerprint">Fingerprint</label><select id="qsar-fingerprint" value={fingerprint} disabled={algorithm === "chemprop_v2" || algorithm === "qsartuna"} onChange={(event) => setFingerprint(event.target.value)}><option value="morgan">Morgan</option><option value="rdkit">RDKit</option><option value="maccs">MACCS</option></select></div>
        <div className="field-group"><label htmlFor="qsar-split">{t("qsar.split")}</label><select id="qsar-split" value={splitMethod} onChange={(event) => setSplitMethod(event.target.value)}><option value="scaffold">Scaffold Split</option><option value="random">Random Split</option></select><p className="field-hint">{t("qsar.scaffoldHint")}</p></div>
        <div className="field-group"><label htmlFor="qsar-ratio">{t("qsar.trainingRatio")}</label><input id="qsar-ratio" type="number" min="0.5" max="0.95" step="0.05" value={trainFraction} onChange={(event) => setTrainFraction(Number(event.target.value))} /></div>
        <div className="field-group"><label htmlFor="qsar-folds">{t("qsar.folds")}</label><input id="qsar-folds" type="number" min="2" max="10" value={folds} onChange={(event) => setFolds(Number(event.target.value))} /></div>
        <div className="field-group"><label htmlFor="qsar-ad">{t("qsar.adThreshold")}</label><input id="qsar-ad" type="number" min="0" max="1" step="0.05" value={adThreshold} onChange={(event) => setAdThreshold(Number(event.target.value))} /><p className="field-hint">{t("qsar.adHint")}</p></div>
        <div className="field-group wide"><label htmlFor="qsar-comment">{t("qsar.modelComment")}</label><input id="qsar-comment" value={comment} onChange={(event) => setComment(event.target.value.replace(/[^A-Za-z0-9_-]/g, ""))} placeholder="kinase_activity" maxLength={50} /><p className="field-hint">{t("qsar.modelCommentHint")}</p></div>
      </div>
      <details className="advanced-settings" open={advanced} onToggle={(event) => setAdvanced(event.currentTarget.open)}><summary>Advanced Settings</summary><div className="qsar-grid">
        {(algorithm === "random_forest" || algorithm === "extra_trees") && <><div className="field-group"><label>n_estimators</label><input type="number" min="1" max="10000" value={hyperparameters.n_estimators ?? 300} onChange={(event) => setHyperparameters((value) => ({ ...value, n_estimators: Number(event.target.value) }))} /></div><div className="field-group"><label>max_depth</label><input type="number" min="1" max="1000" value={hyperparameters.max_depth ?? 20} onChange={(event) => setHyperparameters((value) => ({ ...value, max_depth: Number(event.target.value) }))} /></div><div className="field-group"><label>min_samples_leaf</label><input type="number" min="1" value={hyperparameters.min_samples_leaf ?? 1} onChange={(event) => setHyperparameters((value) => ({ ...value, min_samples_leaf: Number(event.target.value) }))} /></div>{algorithm === "random_forest" && <div className="field-group"><label>min_samples_split</label><input type="number" min="2" value={hyperparameters.min_samples_split ?? 2} onChange={(event) => setHyperparameters((value) => ({ ...value, min_samples_split: Number(event.target.value) }))} /></div>}</>}
        {algorithm === "svm" && <><div className="field-group"><label>C</label><input type="number" min="0.000001" step="0.1" value={hyperparameters.C ?? 1} onChange={(event) => setHyperparameters((value) => ({ ...value, C: Number(event.target.value) }))} /></div><div className="field-group"><label>kernel</label><select value={hyperparameters.kernel ?? "rbf"} onChange={(event) => setHyperparameters((value) => ({ ...value, kernel: event.target.value }))}><option>rbf</option><option>linear</option><option>poly</option><option>sigmoid</option></select></div><div className="field-group"><label>gamma</label><select value={hyperparameters.gamma ?? "scale"} onChange={(event) => setHyperparameters((value) => ({ ...value, gamma: event.target.value }))}><option>scale</option><option>auto</option></select></div>{taskType === "regression" && <div className="field-group"><label>epsilon</label><input type="number" min="0" step="0.01" value={hyperparameters.epsilon ?? 0.1} onChange={(event) => setHyperparameters((value) => ({ ...value, epsilon: Number(event.target.value) }))} /></div>}</>}
        {algorithm === "gradient_boosting" && <><div className="field-group"><label>n_estimators</label><input type="number" min="1" value={hyperparameters.n_estimators ?? 200} onChange={(event) => setHyperparameters((value) => ({ ...value, n_estimators: Number(event.target.value) }))} /></div><div className="field-group"><label>learning_rate</label><input type="number" min="0.000001" max="1" step="0.01" value={hyperparameters.learning_rate ?? 0.05} onChange={(event) => setHyperparameters((value) => ({ ...value, learning_rate: Number(event.target.value) }))} /></div><div className="field-group"><label>max_depth</label><input type="number" min="1" value={hyperparameters.max_depth ?? 3} onChange={(event) => setHyperparameters((value) => ({ ...value, max_depth: Number(event.target.value) }))} /></div></>}
        {algorithm === "lightgbm" && <><NumberParameter name="n_estimators" value={hyperparameters.n_estimators ?? 300} update={setHyperparameters} min={1} /><NumberParameter name="learning_rate" value={hyperparameters.learning_rate ?? 0.05} update={setHyperparameters} min={0.000001} max={1} step={0.01} /><NumberParameter name="max_depth" value={hyperparameters.max_depth ?? -1} update={setHyperparameters} min={-1} /><NumberParameter name="num_leaves" value={hyperparameters.num_leaves ?? 31} update={setHyperparameters} min={2} /><NumberParameter name="min_child_samples" value={hyperparameters.min_child_samples ?? 20} update={setHyperparameters} min={1} /><NumberParameter name="subsample" value={hyperparameters.subsample ?? 1} update={setHyperparameters} min={0.01} max={1} step={0.05} /><NumberParameter name="colsample_bytree" value={hyperparameters.colsample_bytree ?? 1} update={setHyperparameters} min={0.01} max={1} step={0.05} /></>}
        {algorithm === "xgboost" && <><NumberParameter name="n_estimators" value={hyperparameters.n_estimators ?? 300} update={setHyperparameters} min={1} /><NumberParameter name="learning_rate" value={hyperparameters.learning_rate ?? 0.05} update={setHyperparameters} min={0.000001} max={1} step={0.01} /><NumberParameter name="max_depth" value={hyperparameters.max_depth ?? 6} update={setHyperparameters} min={1} /><NumberParameter name="min_child_weight" value={hyperparameters.min_child_weight ?? 1} update={setHyperparameters} min={0} /><NumberParameter name="subsample" value={hyperparameters.subsample ?? 1} update={setHyperparameters} min={0.01} max={1} step={0.05} /><NumberParameter name="colsample_bytree" value={hyperparameters.colsample_bytree ?? 1} update={setHyperparameters} min={0.01} max={1} step={0.05} /><NumberParameter name="reg_alpha" value={hyperparameters.reg_alpha ?? 0} update={setHyperparameters} min={0} /><NumberParameter name="reg_lambda" value={hyperparameters.reg_lambda ?? 1} update={setHyperparameters} min={0} /></>}
        {algorithm === "catboost" && <><NumberParameter name="iterations" value={hyperparameters.iterations ?? 300} update={setHyperparameters} min={1} /><NumberParameter name="learning_rate" value={hyperparameters.learning_rate ?? 0.05} update={setHyperparameters} min={0.000001} max={1} step={0.01} /><NumberParameter name="depth" value={hyperparameters.depth ?? 6} update={setHyperparameters} min={1} max={16} /><NumberParameter name="l2_leaf_reg" value={hyperparameters.l2_leaf_reg ?? 3} update={setHyperparameters} min={0} /></>}
        {(algorithm === "kernel_ridge" || algorithm === "gaussian_process") && <><div className="field-group"><label>kernel</label><select value={hyperparameters.kernel ?? "rbf"} onChange={(event) => setHyperparameters((value) => ({ ...value, kernel: event.target.value }))}><option value="rbf">RBF</option>{algorithm === "kernel_ridge" && <><option value="linear">Linear</option><option value="poly">Polynomial</option><option value="laplacian">Laplacian</option><option value="sigmoid">Sigmoid</option></>} {algorithm === "gaussian_process" && <option value="matern">Matérn</option>}</select></div>{algorithm === "kernel_ridge" && <><NumberParameter name="alpha" value={hyperparameters.alpha ?? 1} update={setHyperparameters} min={0} step={0.1} /><NumberParameter name="gamma" value={hyperparameters.gamma ?? 0.01} update={setHyperparameters} min={0.000001} step={0.01} /></>}</>}
        {(algorithm === "ridge" || algorithm === "elastic_net") && <><NumberParameter name="alpha" value={hyperparameters.alpha ?? 1} update={setHyperparameters} min={0} step={0.1} />{algorithm === "elastic_net" && <NumberParameter name="l1_ratio" value={hyperparameters.l1_ratio ?? 0.5} update={setHyperparameters} min={0} max={1} step={0.05} />}</>}
        {algorithm === "logistic_regression" && <><NumberParameter name="C" value={hyperparameters.C ?? 1} update={setHyperparameters} min={0.000001} step={0.1} /><NumberParameter name="max_iter" value={hyperparameters.max_iter ?? 5000} update={setHyperparameters} min={1} /></>}
        {algorithm === "mlp" && <><div className="field-group"><label>hidden_layer_sizes</label><input value={hyperparameters.hidden_layer_sizes ?? "128,64"} onChange={(event) => setHyperparameters((value) => ({ ...value, hidden_layer_sizes: event.target.value }))} /><p className="field-hint">{t("qsar.hiddenLayersHint")}</p></div><NumberParameter name="alpha" value={hyperparameters.alpha ?? 0.0001} update={setHyperparameters} min={0} step={0.0001} /><NumberParameter name="learning_rate_init" value={hyperparameters.learning_rate_init ?? 0.001} update={setHyperparameters} min={0.000001} max={1} step={0.001} /><NumberParameter name="max_iter" value={hyperparameters.max_iter ?? 500} update={setHyperparameters} min={1} /></>}
        {algorithm === "chemprop_v2" && <><NumberParameter name="epochs" value={hyperparameters.epochs ?? 30} update={setHyperparameters} min={1} /><NumberParameter name="batch_size" value={hyperparameters.batch_size ?? 64} update={setHyperparameters} min={1} /><NumberParameter name="learning_rate" value={hyperparameters.learning_rate ?? 0.0001} update={setHyperparameters} min={0.000001} max={1} step={0.0001} /><NumberParameter name="patience" value={hyperparameters.patience ?? 10} update={setHyperparameters} min={1} /></>}
        {algorithm === "qsartuna" && <><NumberParameter name="n_trials" value={hyperparameters.n_trials ?? 30} update={setHyperparameters} min={1} /><NumberParameter name="n_startup_trials" value={hyperparameters.n_startup_trials ?? 10} update={setHyperparameters} min={0} /><NumberParameter name="timeout_seconds" value={hyperparameters.timeout_seconds ?? 3600} update={setHyperparameters} min={1} /></>}
      </div></details>
      {error && <p className="notice error" role="alert">{error}</p>}
      <button className="button primary large" type="button" onClick={() => void start()} disabled={busy || !content || !smilesColumn || !targetColumn || !comment || !inspection?.valid_count || (mode === "comparison" && selectedAlgorithms.length < 2)}>{busy ? t("qsar.busy") : mode === "comparison" ? t("qsar.compare") : t("qsar.start")}</button>

      {comparison && <section className="model-management comparison-analysis-section">
        <div className="section-heading compact">
          <div><p className="section-kicker">MODEL COMPARISON</p><h3>{t("qsar.results")}</h3><p>{t("qsar.recommended", { algorithm: comparison.recommended_algorithm, metric: comparison.primary_metric })}</p></div>
          <label className="comparison-sort">{t("qsar.sortOrder")}<select value={comparisonSort} onChange={(event) => setComparisonSort(event.target.value as typeof comparisonSort)}><option value="metric">{t("qsar.primary")}</option><option value="algorithm">{t("qsar.modelName")}</option></select></label>
        </div>
        <div className="qsar-preview comparison-result-table-wrap"><table>
          <thead><tr><th>{t("qsar.model")}</th><th>{t("qsar.status")}</th><th>{t("qsar.cvMean")}</th><th>Test</th><th>{t("qsar.testR2")}</th><th>{t("qsar.gap")}</th><th>{t("qsar.inAd")}</th><th>{t("qsar.register")}</th></tr></thead>
          <tbody>{sortComparison(comparison, comparisonSort).map((result) => {
            const cv = result.metrics?.cross_validation?.[comparison.primary_metric];
            const std = result.metrics?.cross_validation?.[`${comparison.primary_metric}_std`];
            const test = result.metrics?.test?.[comparison.primary_metric];
            const testR2 = result.metrics?.test?.r2;
            const gap = result.metrics?.performance_gap?.[comparison.primary_metric];
            return <tr key={result.algorithm}>
              <td>{result.algorithm}{result.algorithm === comparison.recommended_algorithm ? " ★" : ""}</td>
              <td>{result.status === "completed" ? t("qsar.success") : t("qsar.failure", { error: result.error ?? "" })}</td>
              <td>{typeof cv === "number" ? `${cv.toFixed(4)} ± ${typeof std === "number" ? std.toFixed(4) : "—"}` : "—"}</td>
              <td>{typeof test === "number" ? test.toFixed(4) : "—"}</td>
              <td>{typeof testR2 === "number" ? testR2.toFixed(4) : "—"}</td>
              <td>{typeof gap === "number" ? gap.toFixed(4) : "—"}</td>
              <td>{typeof result.metrics?.applicability_domain_fraction === "number" ? `${(result.metrics.applicability_domain_fraction * 100).toFixed(1)}%` : "—"}</td>
              <td>{result.status === "completed" && comparisonJobId && <button type="button" className="button quiet" onClick={() => void registerQsarComparisonModel(comparisonJobId, result.algorithm, `${comment}_${result.algorithm}`.slice(0, 50)).then(refreshModels).catch((reason) => setError(message(reason, t("qsar.error"))))}>{t("qsar.registerModel")}</button>}</td>
            </tr>;
          })}</tbody>
        </table></div>
        <ComparisonCharts comparison={comparison} />
      </section>}

      <section className="model-management">
        <div className="section-heading compact"><div><p className="section-kicker">MODEL REGISTRY</p><h3>{t("qsar.registry")}</h3></div><button type="button" className="button quiet" onClick={() => void refreshModels()}>{t("jobs.refresh")}</button></div>
        {models.length === 0 ? <div className="empty-state">{t("qsar.noModels")}</div> : <div className="qsar-model-list">{models.map((model) => <article className={`qsar-model-card ${predictionModelId === model.id ? "prediction-open" : ""}`} key={model.id}><div><strong>{model.name}</strong><span>{model.algorithm} · {model.task_type} · {model.target_name} ({model.target_unit}) · Train {model.training_data_count} / Test {model.test_data_count}</span><span>{model.task_type === "regression" ? `R² ${metric(model, "r2")} · RMSE ${metric(model, "rmse")} · MAE ${metric(model, "mae")}` : `ROC-AUC ${metric(model, "roc_auc")} · F1 ${metric(model, "f1")}`}</span>{evaluation?.modelId === model.id && <div className="qsar-plots">{Object.entries(evaluation.plots).map(([name, points]) => <MiniPlot key={name} name={name} points={points} />)}{Object.keys(evaluation.plots).length === 0 && <span>{t("qsar.noPlots")}</span>}</div>}</div><div className="qsar-model-actions"><button className="button primary" type="button" onClick={() => openPrediction(model.id)}>{predictionModelId === model.id ? t("qsar.closeEvaluation") : t("qsar.evaluate")}</button><button className="button quiet" type="button" onClick={() => void getQsarEvaluation(model.id).then((result) => setEvaluation({ modelId: model.id, plots: result.plots ?? {} })).catch((reason) => setError(message(reason, t("qsar.error"))))}>{t("qsar.details")}</button><a className="button quiet" href={qsarArtifactUrl(model.id, "evaluation")}>{t("qsar.evaluationJson")}</a><a className="button quiet" href={qsarArtifactUrl(model.id, "predictions")}>{t("qsar.predictions")}</a><a className="button quiet" href={qsarArtifactUrl(model.id, "model")}>{t("qsar.model")}</a><button className="button danger" type="button" disabled={model.in_use} onClick={() => void deleteQsarModel(model.id).then(refreshModels).catch((reason) => setError(message(reason, t("qsar.error"))))}>{t("common.delete")}</button></div>{predictionModelId === model.id && <section className="qsar-prediction-panel"><div><p className="section-kicker">SINGLE MOLECULE PREDICTION</p><h4>{t("qsar.evaluateTitle", { name: model.name })}</h4><p className="muted">{t("qsar.evaluateHint", { target: model.target_name, unit: model.target_unit })}</p></div><ParentInput smiles={predictionSmiles} modelKey="reinvent" validation={predictionValidation} busy={predictionBusy} error={predictionError} onSmilesChange={(value) => { setPredictionSmiles(value); setPredictionValidation(null); setPredictionResult(null); setPredictionError(null); }} onValidate={() => void validatePredictionSmiles()} inputLabel={t("qsar.evaluationSmiles")} inputHint={t("qsar.evaluationSmilesHint")} /><button className="button primary large" type="button" onClick={() => void calculatePrediction()} disabled={predictionBusy || !predictionSmiles.trim()}>{predictionBusy ? t("qsar.predicting") : t("qsar.calculatePrediction")}</button>{predictionResult && <div className={`qsar-prediction-result ${predictionResult.valid ? "valid" : "invalid"}`} role="status"><div><span>{t("qsar.predictedValue")}</span><strong>{predictionResult.prediction == null ? "—" : predictionResult.prediction.toFixed(5)}</strong><small>{model.task_type === "classification" ? t("qsar.classProbability") : model.target_unit}</small></div><div><span>{t("qsar.applicabilityDomain")}</span><strong>{predictionResult.in_ad ? "In" : "Out"}</strong><small>{t("qsar.nearestSimilarity", { value: predictionResult.nearest_tanimoto == null ? "—" : predictionResult.nearest_tanimoto.toFixed(3), threshold: model.ad_threshold.toFixed(3) })}</small></div><div><span>{t("qsar.standardizedStructure")}</span><code>{predictionResult.standardized_smiles ?? "—"}</code><small>MW {predictionResult.molecular_weight?.toFixed(2) ?? "—"} · Charge {predictionResult.formal_charge ?? "—"}</small></div>{predictionResult.error && <p className="notice error">{predictionResult.error}</p>}</div>}</section>}</article>)}</div>}
      </section>
    </div>
  );
}

function NumberParameter({ name, value, update, min, max, step = 1 }: { name: string; value: string | number; update: React.Dispatch<React.SetStateAction<HyperParameters>>; min?: number; max?: number; step?: number }) {
  return <div className="field-group"><label>{name}</label><input type="number" value={value} min={min} max={max} step={step} onChange={(event) => update((current) => ({ ...current, [name]: Number(event.target.value) }))} /></div>;
}

function MiniPlot({ name, points }: { name: string; points: number[][] }) {
  if (!points.length) return null;
  const xs = points.map((point) => point[0]);
  const ys = points.map((point) => point[1]);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const mapX = (value: number) => 12 + (value - minX) / (maxX - minX || 1) * 176;
  const mapY = (value: number) => 88 - (value - minY) / (maxY - minY || 1) * 76;
  return <figure><figcaption>{name.replaceAll("_", " ")}</figcaption><svg viewBox="0 0 200 100" role="img" aria-label={name}><line x1="12" y1="88" x2="188" y2="88" /><line x1="12" y1="12" x2="12" y2="88" />{points.map((point, index) => <circle key={index} cx={mapX(point[0])} cy={mapY(point[1])} r="2.5" />)}</svg></figure>;
}

function sortComparison(comparison: QsarComparisonResult, sort: "metric" | "algorithm") {
  const results = [...comparison.results];
  if (sort === "algorithm") return results.sort((a, b) => a.algorithm.localeCompare(b.algorithm));
  const direction = ["rmse", "mae"].includes(comparison.primary_metric) ? 1 : -1;
  return results.sort((a, b) => {
    const av = a.metrics?.cross_validation?.[comparison.primary_metric];
    const bv = b.metrics?.cross_validation?.[comparison.primary_metric];
    return direction * ((typeof av === "number" ? av : Number.POSITIVE_INFINITY) - (typeof bv === "number" ? bv : Number.POSITIVE_INFINITY));
  });
}

function ComparisonCharts({ comparison }: { comparison: QsarComparisonResult }) {
  const { t } = useTranslation();
  const completed = comparison.results.filter((result) => result.status === "completed");
  const [selectedAlgorithm, setSelectedAlgorithm] = useState(comparison.recommended_algorithm);
  useEffect(() => setSelectedAlgorithm(comparison.recommended_algorithm), [comparison]);
  const values = completed.map((result) => ({ name: result.algorithm, value: Number(result.metrics?.cross_validation?.[comparison.primary_metric]) }));
  const finite = values.filter((item) => Number.isFinite(item.value));
  const max = Math.max(...finite.map((item) => Math.abs(item.value)), 1);
  const selected = completed.find((result) => result.algorithm === selectedAlgorithm) ?? completed[0];
  const folds = selected?.metrics?.cross_validation?.[`${comparison.primary_metric}_values`];
  return <div className="comparison-charts">
    <label className="comparison-model-selector">{t("qsar.distributionModel")}<select value={selected?.algorithm ?? ""} onChange={(event) => setSelectedAlgorithm(event.target.value)}>{completed.map((result) => <option key={result.algorithm} value={result.algorithm}>{result.algorithm}{result.algorithm === comparison.recommended_algorithm ? ` (${t("qsar.recommendedShort")})` : ""}</option>)}</select></label>
    <figure><figcaption>Cross-validation {comparison.primary_metric}</figcaption><div className="comparison-bars">{finite.map((item) => <div key={item.name}><span>{item.name}</span><i style={{ width: `${Math.max(2, Math.abs(item.value) / max * 100)}%` }} /><b>{item.value.toFixed(4)}</b></div>)}</div></figure>
    {Array.isArray(folds) && selected && <MiniPlot name={`${selected.algorithm}_cv_${comparison.primary_metric}`} points={folds.map((value, index) => [index + 1, value])} />}
    {selected?.plots && Object.entries(selected.plots).map(([name, points]) => <MiniPlot key={name} name={`${selected.algorithm}_${name}`} points={points} />)}
  </div>;
}
