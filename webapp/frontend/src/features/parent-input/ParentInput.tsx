import { lazy, Suspense, useEffect, useState } from "react";
import { useTranslation } from "react-i18next";

import { ApiError, lookupPubChem, type ParentValidationResponse } from "../../api/client";
import { ErrorBoundary } from "../../components/ErrorBoundary";

const KetcherEditor = lazy(() => import("./KetcherEditor"));

interface Props {
  smiles: string;
  modelKey: string;
  validation: ParentValidationResponse | null;
  busy: boolean;
  error: string | null;
  onSmilesChange: (smiles: string) => void;
  onValidate: () => void;
  inputLabel?: string;
  inputHint?: string;
}

export function ParentInput({
  smiles,
  modelKey,
  validation,
  busy,
  error,
  onSmilesChange,
  onValidate,
  inputLabel: inputLabelOverride,
  inputHint: inputHintOverride
}: Props) {
  const { t } = useTranslation();
  const [showEditor, setShowEditor] = useState(false);
  const [editorStatus, setEditorStatus] = useState<"loading" | "ready" | "failed">("loading");
  const [editorMessage, setEditorMessage] = useState<string | null>(null);
  const [syncVersion, setSyncVersion] = useState(0);
  const [molecularName, setMolecularName] = useState("");
  const [pubchemBusy, setPubchemBusy] = useState(false);
  const [pubchemError, setPubchemError] = useState<string | null>(null);
  const [pubchemResult, setPubchemResult] = useState<string | null>(null);
  const isLibInvent = modelKey === "libinvent";
  const isLinkInvent = modelKey === "linkinvent";

  const inputLabel = inputLabelOverride ?? (isLinkInvent
    ? "LinkInvent warheads"
    : isLibInvent
      ? "LibInvent scaffold"
      : "Parent SMILES");
  const placeholder = isLinkInvent
    ? t("parent.linkPlaceholder")
    : isLibInvent
      ? t("parent.libPlaceholder")
      : t("parent.placeholder");
  const inputHint = inputHintOverride ?? (isLinkInvent
    ? t("parent.linkHint")
    : isLibInvent
      ? t("parent.libHint")
      : t("parent.hint"));

  useEffect(() => {
    if (showEditor && validation) setSyncVersion((value) => value + 1);
  }, [showEditor, validation]);

  useEffect(() => {
    if (showEditor && editorStatus === "ready" && smiles.trim()) {
      setSyncVersion((value) => value + 1);
    }
  }, [modelKey]);

  useEffect(() => {
    if (!showEditor || editorStatus !== "loading") return;
    const timeout = window.setTimeout(() => {
      setEditorStatus("failed");
      setEditorMessage(t("parent.editorTimeout"));
    }, 15_000);
    return () => window.clearTimeout(timeout);
  }, [editorStatus, showEditor, t]);

  const handlePubChemLookup = async () => {
    setPubchemBusy(true);
    setPubchemError(null);
    setPubchemResult(null);
    try {
      const result = await lookupPubChem(molecularName.trim());
      onSmilesChange(result.isomeric_smiles);
      setMolecularName(result.molecular_name);
      setPubchemResult(`PubChem CID ${result.cid} · ${result.molecular_name}`);
    } catch (lookupError) {
      setPubchemError(
        lookupError instanceof ApiError
          ? lookupError.message
          : t("parent.pubchemError")
      );
    } finally {
      setPubchemBusy(false);
    }
  };

  return (
    <div className="parent-input">
      <div className="field-group pubchem-field">
        <label htmlFor="molecular-name">Molecular Name</label>
        <div className="pubchem-lookup-row">
          <input
            id="molecular-name"
            value={molecularName}
            onChange={(event) => {
              setMolecularName(event.target.value);
              setPubchemError(null);
              setPubchemResult(null);
            }}
            onKeyDown={(event) => {
              if (event.key === "Enter" && molecularName.trim() && !pubchemBusy) {
                event.preventDefault();
                void handlePubChemLookup();
              }
            }}
            placeholder={t("parent.namePlaceholder")}
            autoComplete="off"
          />
          <button
            className="button secondary"
            type="button"
            onClick={() => void handlePubChemLookup()}
            disabled={pubchemBusy || !molecularName.trim()}
          >
            {pubchemBusy ? t("parent.searching") : t("parent.lookup")}
          </button>
        </div>
        <p className="field-hint">{t("parent.pubchemHint")}</p>
        {pubchemError && <p className="notice error compact" role="alert">{pubchemError}</p>}
        {pubchemResult && <p className="pubchem-result" role="status">✓ {t("parent.applied", { value: pubchemResult })}</p>}
      </div>

      <div className="field-group">
        <label htmlFor="parent-smiles">{inputLabel}</label>
        <textarea
          id="parent-smiles"
          value={smiles}
          onChange={(event) => onSmilesChange(event.target.value)}
          placeholder={placeholder}
          rows={2}
          className="smiles-input"
          spellCheck={false}
        />
        <p className="field-hint">
          {inputHint}
        </p>
      </div>

      <div className="inline-actions">
        <button className="button primary" type="button" onClick={onValidate} disabled={busy}>
          {busy ? t("parent.validating") : t("parent.validate")}
        </button>
        <button
          className="button secondary"
          type="button"
          onClick={() => setShowEditor((value) => !value)}
        >
          {showEditor ? t("parent.closeEditor") : t("parent.openEditor")}
        </button>
        {showEditor && (
          <button
            className="button quiet"
            type="button"
            onClick={() => setSyncVersion((value) => value + 1)}
            disabled={!smiles.trim() || editorStatus !== "ready"}
          >
            {t("parent.sync")}
          </button>
        )}
      </div>

      {error && <p className="notice error" role="alert">{error}</p>}
      {validation && (
        <div className="validation-card" role="status">
          <span className="validation-mark">✓</span>
          <div>
            <strong>{t("parent.valid")}</strong>
            <p>{validation.canonical_smiles}</p>
            <small>{validation.molecular_formula} · {validation.atom_count} atoms</small>
          </div>
        </div>
      )}

      {showEditor && (
        <div className="editor-shell" aria-label={t("parent.editorLabel")}>
          <div className="editor-caption">
            <span>Ketcher standalone</span>
            <span className={`editor-state ${editorStatus}`}>{editorStatus}</span>
          </div>
          {editorMessage && <p className="notice error" role="alert">{editorMessage}</p>}
          <ErrorBoundary
            fallback={
              <div className="editor-fallback" role="alert">
                {t("parent.editorUnavailable")}
              </div>
            }
            onError={() => setEditorStatus("failed")}
          >
            <Suspense fallback={<div className="editor-fallback">{t("parent.editorLoading")}</div>}>
              <KetcherEditor
                value={smiles}
                linkInvent={isLinkInvent}
                syncVersion={syncVersion}
                onChange={onSmilesChange}
                onError={(message) => {
                  setEditorStatus("failed");
                  setEditorMessage(message);
                }}
                onReady={() => {
                  setEditorStatus("ready");
                  if (smiles.trim()) setSyncVersion((value) => value + 1);
                }}
              />
            </Suspense>
          </ErrorBoundary>
        </div>
      )}
    </div>
  );
}
