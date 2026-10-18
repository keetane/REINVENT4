import { useEffect, useMemo, useState } from "react";
import { useTranslation } from "react-i18next";

import {
  ApiError, createRlJob, getRlCapabilities, listQsarModels, previewRlScoring,
  type Job, type ModelInfo, type QsarModel, type RlPreviewRow, type RlPropertyComponent, type RlQsarComponent
} from "../../api/client";
import { StandardizationNotice } from "../standardization/StandardizationNotice";

interface Props { models: ModelInfo[]; onCreated: (job: Job) => void; }

type QsarGoal = "maximize" | "minimize" | "range" | "target";
type QsarSetting = RlQsarComponent & { goal: QsarGoal; target: number; tolerance: number };

const propertyOptions = [
  ["QED", "QED"], ["MolecularWeight", "MW"], ["SlogP", "SlogP"], ["TPSA", "TPSA"],
  ["HBondDonors", "HBD"], ["HBondAcceptors", "HBA"], ["NumRotBond", "Rotatable Bonds"], ["SAScore", "SA Score"]
] as const;

function errorMessage(error: unknown, fallback: string) {
  return error instanceof ApiError ? error.message : fallback;
}

export function RlPanel({ models, onCreated }: Props) {
  const { t } = useTranslation();
  const generationModels = useMemo(() => models.filter((model) =>
    ["reinvent", "mol2mol"].includes(model.model_type) &&
    (model.source !== "standard" || ["reinvent", "mol2mol_scaffold_generic"].includes(model.base_model_key ?? ""))
  ), [models]);
  const [generationModelId, setGenerationModelId] = useState("standard:reinvent");
  const [qsarModels, setQsarModels] = useState<QsarModel[]>([]);
  const [selectedQsar, setSelectedQsar] = useState<string[]>([]);
  const [qsarSettings, setQsarSettings] = useState<Record<string, QsarSetting>>({});
  const [properties, setProperties] = useState<string[]>(["QED"]);
  const [comment, setComment] = useState("");
  const [device, setDevice] = useState<"cpu" | "cuda" | "mps">("mps");
  const [devices, setDevices] = useState<{ key: "cpu" | "cuda" | "mps"; label: string; available: boolean }[]>([
    { key: "cpu", label: "CPU", available: true }
  ]);
  const [inputSmiles, setInputSmiles] = useState("");
  const [batchSize, setBatchSize] = useState(64);
  const [sigma, setSigma] = useState(128);
  const [learningRate, setLearningRate] = useState(0.0001);
  const [minimumSteps, setMinimumSteps] = useState(25);
  const [maximumSteps, setMaximumSteps] = useState(300);
  const [targetScore, setTargetScore] = useState(0.8);
  const [diversity, setDiversity] = useState(true);
  const [previewSmiles, setPreviewSmiles] = useState("CCO\nCCN");
  const [preview, setPreview] = useState<RlPreviewRow[]>([]);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const selectedGeneration = generationModels.find((model) => model.id === generationModelId);
  const mol2mol = selectedGeneration?.model_type === "mol2mol";

  useEffect(() => {
    void Promise.all([
      listQsarModels().then((response) => setQsarModels(response.items)),
      getRlCapabilities().then((response) => {
        setDevices(response.devices);
        if (!response.devices.some((item) => item.key === "mps" && item.available)) {
          setDevice(response.devices.find((item) => item.available)?.key ?? "cpu");
        }
      })
    ]).catch((reason) => setError(errorMessage(reason, t("rl.error"))));
  }, [t]);

  const defaultSetting = (modelId: string): QsarSetting => ({
    model_id: modelId, weight: 1, transform: "sigmoid", low: 0, high: 10, k: 0.5,
    class_index: 1, ad_behavior: "penalize", goal: "maximize", target: 5, tolerance: 1
  });
  const settingFor = (modelId: string) => qsarSettings[modelId] ?? defaultSetting(modelId);
  const qsarComponents: RlQsarComponent[] = selectedQsar.map((modelId) => {
    const setting = settingFor(modelId);
    const transform = setting.goal === "maximize" ? "sigmoid" : setting.goal === "minimize" ? "reverse_sigmoid" : "double_sigmoid";
    const low = setting.goal === "target" ? setting.target - setting.tolerance : setting.low;
    const high = setting.goal === "target" ? setting.target + setting.tolerance : setting.high;
    return { model_id: modelId, weight: setting.weight, transform, low, high, k: setting.k, class_index: setting.class_index, ad_behavior: setting.ad_behavior };
  });
  const propertyComponents: RlPropertyComponent[] = properties.map((component) => ({
    component, weight: 1, transform: component === "QED" ? "none" : "double_sigmoid",
    low: component === "MolecularWeight" ? 200 : 0,
    high: component === "MolecularWeight" ? 500 : component === "TPSA" ? 140 : 10,
    k: 0.5, reference_smiles: [], smarts: []
  }));

  const toggleQsar = (modelId: string) => {
    setSelectedQsar((current) => current.includes(modelId) ? current.filter((id) => id !== modelId) : [...current, modelId]);
    setQsarSettings((current) => current[modelId] ? current : { ...current, [modelId]: defaultSetting(modelId) });
  };

  const updateQsar = (modelId: string, patch: Partial<QsarSetting>) => {
    setQsarSettings((current) => ({ ...current, [modelId]: { ...(current[modelId] ?? defaultSetting(modelId)), ...patch } }));
  };

  const start = async () => {
    setBusy(true); setError(null);
    try {
      const job = await createRlJob({
        generation_model_id: generationModelId, comment, device, input_smiles_content: inputSmiles,
        qsar_components: qsarComponents,
        property_components: propertyComponents,
        aggregation: "geometric_mean", batch_size: batchSize, sigma, learning_rate: learningRate,
        minimum_steps: minimumSteps, maximum_steps: maximumSteps, target_mean_score: targetScore,
        diversity_filter: diversity, diversity_bucket_size: 25, diversity_min_score: 0.4
      });
      onCreated(job);
    } catch (reason) { setError(errorMessage(reason, t("rl.error"))); } finally { setBusy(false); }
  };

  const testScoring = async () => {
    setBusy(true); setError(null);
    try {
      const smiles = previewSmiles.split(/\r?\n/).map((value) => value.trim()).filter(Boolean);
      const result = await previewRlScoring({ smiles, qsar_components: qsarComponents, property_components: propertyComponents, aggregation: "geometric_mean" });
      setPreview(result.rows);
    } catch (reason) { setError(errorMessage(reason, t("rl.error"))); } finally { setBusy(false); }
  };

  return <div className="rl-panel">
    <div className="section-heading"><div><p className="section-kicker">REINFORCEMENT LEARNING</p><h2>{t("rl.title")}</h2><p>{t("rl.intro")}</p></div></div>
    <div className="rl-grid">
      <div className="field-group wide"><label htmlFor="rl-model">{t("rl.agent")}</label><select id="rl-model" value={generationModelId} onChange={(event) => setGenerationModelId(event.target.value)}>{generationModels.map((model) => <option key={model.id} value={model.id}>{model.name} · {model.source}</option>)}</select><p className="field-hint">{t("rl.agentHint")}</p></div>
      <div className="field-group"><label htmlFor="rl-device">{t("rl.device")}</label><select id="rl-device" value={device} onChange={(event) => setDevice(event.target.value as typeof device)}>{devices.map((item) => <option key={item.key} value={item.key} disabled={!item.available}>{item.label}{item.available ? "" : t("common.unavailable")}</option>)}</select></div>
      <div className="field-group"><label htmlFor="rl-comment">{t("rl.comment")}</label><input id="rl-comment" value={comment} onChange={(event) => setComment(event.target.value.replace(/[^A-Za-z0-9_-]/g, ""))} placeholder="targetA_rl" /></div>
      {mol2mol && <div className="field-group wide"><label htmlFor="rl-input">{t("rl.mol2molInput")}</label><textarea id="rl-input" rows={5} value={inputSmiles} onChange={(event) => setInputSmiles(event.target.value)} placeholder={t("rl.onePerLine")} /></div>}
    </div>

    <StandardizationNotice />

    <section className="rl-score-section"><h3>{t("rl.scoring")}</h3><p className="field-hint">{t("rl.scoringHint")}</p>
      <div className="rl-check-grid">{qsarModels.map((model) => <label key={model.id}><input type="checkbox" checked={selectedQsar.includes(model.id)} onChange={() => toggleQsar(model.id)} /><span><strong>QSAR: {model.target_name}</strong><small>{model.name} · {model.algorithm} · {model.task_type}</small></span></label>)}</div>
      {qsarModels.length === 0 && <div className="empty-state">{t("rl.noQsar")}</div>}
      {selectedQsar.map((modelId) => { const model = qsarModels.find((item) => item.id === modelId); const setting = settingFor(modelId); return <article className="qsar-score-config" key={modelId}><h4>QSAR: {model?.target_name ?? model?.name}</h4><div className="rl-grid"><div className="field-group"><label>{t("rl.goal")}</label><select value={setting.goal} onChange={(event) => updateQsar(modelId, { goal: event.target.value as QsarGoal })}><option value="maximize">{t("rl.maximize")}</option><option value="minimize">{t("rl.minimize")}</option><option value="range">{t("rl.range")}</option><option value="target">{t("rl.target")}</option></select></div><div className="field-group"><label>Weight</label><input type="number" min="0.01" max="100" step="0.1" value={setting.weight} onChange={(event) => updateQsar(modelId, { weight: Number(event.target.value) })} /></div>{setting.goal === "target" ? <><div className="field-group"><label>{t("rl.target")}</label><input type="number" value={setting.target} onChange={(event) => updateQsar(modelId, { target: Number(event.target.value) })} /></div><div className="field-group"><label>{t("rl.tolerance")}</label><input type="number" min="0.000001" value={setting.tolerance} onChange={(event) => updateQsar(modelId, { tolerance: Number(event.target.value) })} /></div></> : <><div className="field-group"><label>Low</label><input type="number" value={setting.low} onChange={(event) => updateQsar(modelId, { low: Number(event.target.value) })} /></div><div className="field-group"><label>High</label><input type="number" value={setting.high} onChange={(event) => updateQsar(modelId, { high: Number(event.target.value) })} /></div></>}<div className="field-group"><label>{t("rl.adBehavior")}</label><select value={setting.ad_behavior} onChange={(event) => updateQsar(modelId, { ad_behavior: event.target.value as "penalize" | "allow" })}><option value="penalize">{t("rl.penalize")}</option><option value="allow">{t("rl.allow")}</option></select></div></div></article>; })}
      <h4>{t("rl.properties")}</h4><div className="rl-check-grid">{propertyOptions.map(([key, label]) => <label key={key}><input type="checkbox" checked={properties.includes(key)} onChange={() => setProperties((current) => current.includes(key) ? current.filter((item) => item !== key) : [...current, key])} /><span><strong>{label}</strong></span></label>)}</div>
    </section>

    <details className="advanced-settings"><summary>Advanced Settings</summary><div className="rl-grid">
      <div className="field-group"><label>Batch size</label><input type="number" min="1" value={batchSize} onChange={(event) => setBatchSize(Number(event.target.value))} /></div>
      <div className="field-group"><label>Sigma</label><input type="number" min="1" value={sigma} onChange={(event) => setSigma(Number(event.target.value))} /></div>
      <div className="field-group"><label>Learning rate</label><input type="number" min="0.00000001" step="0.0001" value={learningRate} onChange={(event) => setLearningRate(Number(event.target.value))} /></div>
      <div className="field-group"><label>{t("rl.minimumSteps")}</label><input type="number" min="0" value={minimumSteps} onChange={(event) => setMinimumSteps(Number(event.target.value))} /></div>
      <div className="field-group"><label>{t("rl.maximumSteps")}</label><input type="number" min="1" value={maximumSteps} onChange={(event) => setMaximumSteps(Number(event.target.value))} /></div>
      <div className="field-group"><label>{t("rl.targetMean")}</label><input type="number" min="0" max="1" step="0.05" value={targetScore} onChange={(event) => setTargetScore(Number(event.target.value))} /></div>
      <label className="check-row"><input type="checkbox" checked={diversity} onChange={(event) => setDiversity(event.target.checked)} /> Diversity Filter</label>
    </div></details>

    {(selectedQsar.length > 0 || properties.length > 0) && <section className="rl-preview"><h3>{t("rl.preview")}</h3><textarea rows={4} value={previewSmiles} onChange={(event) => setPreviewSmiles(event.target.value)} /><button className="button quiet" type="button" onClick={() => void testScoring()} disabled={busy}>{t("rl.previewButton")}</button>{preview.length > 0 && <div className="qsar-preview"><table><thead><tr><th>Generated</th><th>Standardized</th><th>MW / Charge</th><th>Scoring Function</th><th>{t("rl.prediction")}</th><th>AD</th><th>{t("rl.transformed")}</th><th>Total Score / Reason</th></tr></thead><tbody>{preview.flatMap((row) => row.components.map((component, index) => <tr key={`${row.smiles}-${component.kind}-${component.name}`}><td>{index === 0 ? row.generated_smiles : ""}</td><td>{index === 0 ? row.standardized_smiles ?? "—" : ""}</td><td>{index === 0 ? `${row.molecular_weight?.toFixed(2) ?? "—"} / ${row.formal_charge ?? "—"}` : ""}</td><td>{component.kind === "qsar" ? `QSAR: ${component.name}` : `RDKit: ${component.name}`}</td><td>{typeof component.prediction === "number" ? component.prediction.toFixed(4) : t("common.error")}</td><td>{component.in_ad === null ? "—" : component.in_ad ? "In" : "Out"}</td><td>{component.score.toFixed(4)}</td><td>{index === 0 ? `${row.total_score.toFixed(4)}${row.exclusion_reason ? ` / ${row.exclusion_reason}` : ""}` : ""}</td></tr>))}</tbody></table></div>}</section>}
    {error && <p className="notice error" role="alert">{error}</p>}
    <button className="button primary large" type="button" onClick={() => void start()} disabled={busy || !comment || (!selectedQsar.length && !properties.length) || (mol2mol && !inputSmiles.trim())}>{busy ? t("rl.busy") : t("rl.start")}</button>
  </div>;
}
