import { lazy, Suspense, useEffect, useState } from "react";

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
}

export function ParentInput({
  smiles,
  modelKey,
  validation,
  busy,
  error,
  onSmilesChange,
  onValidate
}: Props) {
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

  const inputLabel = isLinkInvent
    ? "LinkInvent warheads"
    : isLibInvent
      ? "LibInvent scaffold"
      : "Parent SMILES";
  const placeholder = isLinkInvent
    ? "例: Oc1cncc(*)c1|*c1ccoc1"
    : isLibInvent
      ? "例: c1ccc([*])cc1"
      : "例: CCO または N[C@@H](C)C(=O)O";
  const inputHint = isLinkInvent
    ? "2つの固定側鎖にそれぞれ1個の * を付け、| で区切ってください。"
    : isLibInvent
      ? "固定する骨格を入力し、側鎖を発生させる位置に * を付けてください。"
      : "単一分子を入力してください。反応式と複数 fragment は対象外です。";

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
      setEditorMessage(
        "構造エディタの初期化が完了しませんでした。SMILES 入力はそのまま使用できます。"
      );
    }, 15_000);
    return () => window.clearTimeout(timeout);
  }, [editorStatus, showEditor]);

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
          : "PubChem に接続できませんでした。時間をおいて再試行してください。"
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
            placeholder="例: Aspirin, Imatinib"
            autoComplete="off"
          />
          <button
            className="button secondary"
            type="button"
            onClick={() => void handlePubChemLookup()}
            disabled={pubchemBusy || !molecularName.trim()}
          >
            {pubchemBusy ? "検索中…" : "PubChemから取得"}
          </button>
        </div>
        <p className="field-hint">分子名からisomeric SMILESを検索し、下の入力欄へ反映します。</p>
        {pubchemError && <p className="notice error compact" role="alert">{pubchemError}</p>}
        {pubchemResult && <p className="pubchem-result" role="status">✓ {pubchemResult} を反映しました。</p>}
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
          {busy ? "検証中…" : "検証する"}
        </button>
        <button
          className="button secondary"
          type="button"
          onClick={() => setShowEditor((value) => !value)}
        >
          {showEditor ? "構造エディタを閉じる" : "構造エディタを開く"}
        </button>
        {showEditor && (
          <button
            className="button quiet"
            type="button"
            onClick={() => setSyncVersion((value) => value + 1)}
            disabled={!smiles.trim() || editorStatus !== "ready"}
          >
            SMILES を描画へ反映
          </button>
        )}
      </div>

      {error && <p className="notice error" role="alert">{error}</p>}
      {validation && (
        <div className="validation-card" role="status">
          <span className="validation-mark">✓</span>
          <div>
            <strong>有効な入力</strong>
            <p>{validation.canonical_smiles}</p>
            <small>{validation.molecular_formula} · {validation.atom_count} atoms</small>
          </div>
        </div>
      )}

      {showEditor && (
        <div className="editor-shell" aria-label="Ketcher 構造エディタ">
          <div className="editor-caption">
            <span>Ketcher standalone</span>
            <span className={`editor-state ${editorStatus}`}>{editorStatus}</span>
          </div>
          {editorMessage && <p className="notice error" role="alert">{editorMessage}</p>}
          <ErrorBoundary
            fallback={
              <div className="editor-fallback" role="alert">
                構造エディタを読み込めませんでした。上の SMILES 入力はそのまま使用できます。
              </div>
            }
            onError={() => setEditorStatus("failed")}
          >
            <Suspense fallback={<div className="editor-fallback">構造エディタを準備中…</div>}>
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
