import { useEffect, useRef } from "react";
import type { Ketcher } from "ketcher-core";
import { Editor } from "ketcher-react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import "ketcher-react/dist/index.css";

import { fromEditorStructure, toEditorStructure } from "./structureFormat";

const structureService = new StandaloneStructServiceProvider();

interface Props {
  value: string;
  linkInvent: boolean;
  syncVersion: number;
  onChange: (smiles: string) => void;
  onError: (message: string) => void;
  onReady: () => void;
}

export default function KetcherEditor({
  value,
  linkInvent,
  syncVersion,
  onChange,
  onError,
  onReady
}: Props) {
  const ketcherRef = useRef<Ketcher | null>(null);
  const changeHandlerRef = useRef<(() => void) | null>(null);
  const previousSyncVersion = useRef(syncVersion);
  const linkInventRef = useRef(linkInvent);
  const syncingFromTextRef = useRef(false);
  linkInventRef.current = linkInvent;

  useEffect(() => {
    if (previousSyncVersion.current === syncVersion) return;
    previousSyncVersion.current = syncVersion;
    const ketcher = ketcherRef.current;
    if (!ketcher || !value.trim()) return;
    syncingFromTextRef.current = true;
    void ketcher
      .setMolecule(toEditorStructure(value, linkInvent))
      .catch(() => {
        onError("構造エディタへ反映できませんでした。SMILES を確認してください。");
      })
      .finally(() => {
        syncingFromTextRef.current = false;
      });
  }, [linkInvent, onError, syncVersion, value]);

  useEffect(
    () => () => {
      if (ketcherRef.current && changeHandlerRef.current) {
        ketcherRef.current.changeEvent.remove(changeHandlerRef.current);
      }
    },
    []
  );

  const handleInit = (ketcher: Ketcher) => {
    ketcherRef.current = ketcher;
    const handleChange = () => {
      if (syncingFromTextRef.current) return;
      void ketcher
        .getSmiles()
        .then((smiles) => onChange(fromEditorStructure(smiles, linkInventRef.current)))
        .catch(() => onError("描画した構造を SMILES に変換できませんでした。"));
    };
    changeHandlerRef.current = handleChange;
    ketcher.changeEvent.add(handleChange);
    onReady();
  };

  return (
    <Editor
      staticResourcesUrl="/"
      structServiceProvider={structureService}
      onInit={handleInit}
      errorHandler={() => onError("構造エディタで処理エラーが発生しました。")}
      disableMacromoleculesEditor
    />
  );
}
