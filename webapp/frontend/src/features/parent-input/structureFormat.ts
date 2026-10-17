export function toEditorStructure(value: string, linkInvent: boolean): string {
  return linkInvent ? value.replaceAll("|", ".") : value;
}

export function fromEditorStructure(value: string, linkInvent: boolean): string {
  const cxsmilesStart = value.indexOf(" |");
  const smiles = (cxsmilesStart === -1 ? value : value.slice(0, cxsmilesStart)).trim();
  return linkInvent ? smiles.replaceAll(".", "|") : smiles;
}
