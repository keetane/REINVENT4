import { useEffect, useRef } from "react";
import { useTranslation } from "react-i18next";
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
  const { t } = useTranslation();
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
        onError(t("parent.editorSyncError"));
      })
      .finally(() => {
        syncingFromTextRef.current = false;
      });
  }, [linkInvent, onError, syncVersion, t, value]);

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
        .catch(() => onError(t("parent.editorConvertError")));
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
      errorHandler={() => onError(t("parent.editorError"))}
      disableMacromoleculesEditor
    />
  );
}
