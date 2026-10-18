import { useMemo, useState } from "react";
import { useTranslation } from "react-i18next";

import {
  ApiError,
  createTrainingJob,
  inspectTrainingData,
  type Job,
  type ModelInfo,
  type TrainingDataInspection
} from "../../api/client";
import { isDataUploadTooLarge, suggestSmilesColumn } from "../../lib/uploads";
import { StandardizationNotice } from "../standardization/StandardizationNotice";

type Device = "cpu" | "cuda" | "mps";

interface Props {
  models: ModelInfo[];
  onCreated: (job: Job) => void;
  onDeleteModel: (modelId: string) => Promise<void>;
}

function errorMessage(error: unknown, fallback: string): string {
  return error instanceof ApiError ? error.message : fallback;
}

export function TransferLearningPanel({ models, onCreated, onDeleteModel }: Props) {
  const { t, i18n } = useTranslation();
  const [baseModel, setBaseModel] = useState<"reinvent" | "mol2mol">("reinvent");
  const [fileName, setFileName] = useState("");
  const [content, setContent] = useState("");
  const [smilesColumn, setSmilesColumn] = useState<string | null>(null);
  const [inspection, setInspection] = useState<TrainingDataInspection | null>(null);
  const [comment, setComment] = useState("");
  const [trainPercent, setTrainPercent] = useState(80);
  const [epochs, setEpochs] = useState(3);
  const [batchSize, setBatchSize] = useState(50);
  const [learningRate, setLearningRate] = useState(0.0001);
  const [tanimotoLower, setTanimotoLower] = useState(0.7);
  const [tanimotoUpper, setTanimotoUpper] = useState(1.0);
  const [device, setDevice] = useState<Device>("mps");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const trainedModels = useMemo(
    () => models.filter((model) => model.source === "transfer_learning"),
    [models]
  );

  const inspect = async (name: string, text: string, column: string | null) => {
    const result = await inspectTrainingData({
      file_name: name,
      content: text,
      smiles_column: column,
      base_model: baseModel
    });
    setInspection(result);
    return result;
  };

  const handleFile = async (file: File | undefined) => {
    if (!file) return;
    if (isDataUploadTooLarge(file)) {
      setFileName("");
      setContent("");
      setSmilesColumn(null);
      setInspection(null);
      setError(t("upload.tooLarge"));
      return;
    }
    setBusy(true);
    setError(null);
    setFileName(file.name);
    const text = await file.text();
    setContent(text);
    try {
      const first = await inspect(file.name, text, null);
      if (first.columns.length) {
        const preferred = suggestSmilesColumn(first.columns);
        setSmilesColumn(preferred);
        await inspect(file.name, text, preferred);
      } else {
        setSmilesColumn(null);
      }
    } catch (inspectionError) {
      setInspection(null);
      setError(errorMessage(inspectionError, t("tl.error")));
    } finally {
      setBusy(false);
    }
  };

  const handleColumn = async (column: string) => {
    setSmilesColumn(column);
    setBusy(true);
    setError(null);
    try {
      await inspect(fileName, content, column);
    } catch (inspectionError) {
      setInspection(null);
      setError(errorMessage(inspectionError, t("tl.error")));
    } finally {
      setBusy(false);
    }
  };

  const handleSubmit = async () => {
    setBusy(true);
    setError(null);
    try {
      const job = await createTrainingJob({
        file_name: fileName,
        content,
        smiles_column: smilesColumn,
        base_model: baseModel,
        comment,
        train_fraction: trainPercent / 100,
        num_epochs: epochs,
        batch_size: batchSize,
        learning_rate: learningRate,
        tanimoto_lower: tanimotoLower,
        tanimoto_upper: tanimotoUpper,
        device
      });
      onCreated(job);
    } catch (submitError) {
      setError(errorMessage(submitError, t("tl.error")));
    } finally {
      setBusy(false);
    }
  };

  const validComment = /^[A-Za-z0-9_-]{1,50}$/.test(comment);
  const canSubmit = Boolean(
    inspection && inspection.valid_count >= 2 && fileName && validComment && !busy
  );

  return (
    <div className="training-panel">
      <div className="section-heading">
        <div>
          <p className="section-kicker">TRANSFER LEARNING</p>
          <h2>{t("tl.title")}</h2>
          <p>{t("tl.intro")}</p>
        </div>
        <span className="step-indicator">TL</span>
      </div>

      <div className="training-grid">
        <div className="field-group">
          <label htmlFor="tl-base-model">Base Model</label>
          <select id="tl-base-model" value={baseModel} onChange={(event) => {
            const value = event.target.value as typeof baseModel;
            setBaseModel(value);
            if (content) void inspectTrainingData({ file_name: fileName, content, smiles_column: smilesColumn, base_model: value }).then(setInspection).catch((reason) => setError(errorMessage(reason, t("tl.error"))));
          }}>
            <option value="reinvent">Reinvent</option>
            <option value="mol2mol">Mol2Mol</option>
          </select>
          <p className="field-hint">{t("tl.baseHint")}</p>
        </div>
        <div className="field-group">
          <label htmlFor="tl-device">{t("tl.device")}</label>
          <select id="tl-device" value={device} onChange={(event) => setDevice(event.target.value as Device)}>
            <option value="cpu">CPU</option>
            <option value="cuda">GPU</option>
            <option value="mps">MPS</option>
          </select>
          <p className="field-hint">{t("tl.deviceHint")}</p>
        </div>
        <div className="field-group wide training-file">
          <label htmlFor="training-file">{t("tl.data")}</label>
          <input id="training-file" type="file" accept=".smi,.smiles,.txt,.csv,text/plain,text/csv" onChange={(event) => void handleFile(event.target.files?.[0])} />
          <p className="field-hint">{t("tl.dataHint")}</p>
        </div>
        {inspection?.columns.length ? (
          <div className="field-group">
            <label htmlFor="smiles-column">{t("tl.smilesColumn")}</label>
            <select id="smiles-column" value={smilesColumn ?? ""} onChange={(event) => void handleColumn(event.target.value)}>
              {inspection.columns.map((column) => <option key={column} value={column}>{column}</option>)}
            </select>
            <p className="field-hint">{t("tl.smilesHint")}</p>
          </div>
        ) : null}
      </div>

      <StandardizationNotice version={inspection?.standardization_version} />

      {inspection && inspection.total_rows > 0 && (
        <div className="training-summary" role="status">
          <div><span>入力分子</span><strong>{inspection.total_rows}</strong></div>
          <div><span>標準化済み</span><strong>{inspection.standardized_count}</strong></div>
          <div><span>{t("tl.valid")}</span><strong>{inspection.valid_count}</strong></div>
          <div><span>{t("tl.invalid")}</span><strong>{inspection.invalid_count}</strong></div>
          <div><span>{t("tl.duplicates")}</span><strong>{inspection.duplicate_count}</strong></div>
          <div><span>フラグメント抽出</span><strong>{inspection.fragment_extracted_count}</strong></div>
          <div><span>中和</span><strong>{inspection.neutralized_count}</strong></div>
          <div><span>永久電荷保持</span><strong>{inspection.charged_structure_retained_count}</strong></div>
          <div><span>MW≤200除外</span><strong>{inspection.low_molecular_weight_excluded_count}</strong></div>
          <div><span>モデル非互換</span><strong>{inspection.model_compatibility_error_count}</strong></div>
          <div><span>{t("tl.split")}</span><strong>{trainPercent}:{100 - trainPercent}</strong></div>
        </div>
      )}

      <div className="training-grid settings">
        <div className="field-group wide">
          <label htmlFor="model-comment">{t("tl.comment")}</label>
          <input id="model-comment" value={comment} onChange={(event) => setComment(event.target.value)} placeholder={t("tl.commentPlaceholder")} />
          <p className="field-hint">{t("tl.commentHint")}</p>
        </div>
        <div className="field-group">
          <label htmlFor="train-percent">{t("tl.trainRatio")}</label>
          <input id="train-percent" type="number" min={50} max={95} value={trainPercent} onChange={(event) => setTrainPercent(Number(event.target.value))} />
          <p className="field-hint">{t("tl.trainRatioHint")}</p>
        </div>
        <div className="field-group">
          <label htmlFor="epochs">Epoch</label>
          <input id="epochs" type="number" min={1} max={1000} value={epochs} onChange={(event) => setEpochs(Number(event.target.value))} />
          <p className="field-hint">{t("tl.epochHint")}</p>
        </div>
        <div className="field-group">
          <label htmlFor="batch-size">Batch size</label>
          <input id="batch-size" type="number" min={1} max={4096} value={batchSize} onChange={(event) => setBatchSize(Number(event.target.value))} />
          <p className="field-hint">{t("tl.batchHint")}</p>
        </div>
        <div className="field-group">
          <label htmlFor="learning-rate">Learning rate</label>
          <input id="learning-rate" type="number" min={0.00000001} max={1} step={0.00001} value={learningRate} onChange={(event) => setLearningRate(Number(event.target.value))} />
          <p className="field-hint">{t("tl.lrHint")}</p>
        </div>
        {baseModel === "mol2mol" && (
          <>
            <div className="field-group">
              <label htmlFor="tanimoto-lower">{t("tl.lower")}</label>
              <input id="tanimoto-lower" type="number" min={0} max={1} step={0.05} value={tanimotoLower} onChange={(event) => setTanimotoLower(Number(event.target.value))} />
              <p className="field-hint">{t("tl.lowerHint")}</p>
            </div>
            <div className="field-group">
              <label htmlFor="tanimoto-upper">{t("tl.upper")}</label>
              <input id="tanimoto-upper" type="number" min={0} max={1} step={0.05} value={tanimotoUpper} onChange={(event) => setTanimotoUpper(Number(event.target.value))} />
              <p className="field-hint">{t("tl.upperHint")}</p>
            </div>
          </>
        )}
      </div>

      {error && <p className="notice error" role="alert">{error}</p>}
      {!validComment && comment && <p className="notice warning">{t("tl.invalidComment")}</p>}
      <div className="submit-row">
        <p>{t("tl.progressHint")}</p>
        <button className="button primary large" type="button" disabled={!canSubmit} onClick={() => void handleSubmit()}>
          {busy ? t("tl.starting") : t("tl.start")}
        </button>
      </div>

      <div className="model-management">
        <div className="section-heading compact">
          <div><p className="section-kicker">MODEL REGISTRY</p><h3>{t("tl.registry")}</h3></div>
        </div>
        {trainedModels.length === 0 ? (
          <div className="empty-state">{t("tl.empty")}</div>
        ) : (
          <div className="model-list">
            {trainedModels.map((model) => (
              <article className="model-row" key={model.id}>
                <div>
                  <strong>{model.name}</strong>
                  <span>{model.model_type.toUpperCase()} · {model.training_data_count} · {model.created_at ? new Date(model.created_at).toLocaleString(i18n.language.startsWith("en") ? "en-US" : "ja-JP") : "—"}</span>
                </div>
                <button className="button danger" type="button" onClick={() => {
                  if (window.confirm(t("tl.deleteConfirm", { name: model.name }))) void onDeleteModel(model.id);
                }}>{t("common.delete")}</button>
              </article>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
