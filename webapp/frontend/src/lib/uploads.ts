export const MAX_DATA_UPLOAD_BYTES = 50 * 1024 * 1024;

export function isDataUploadTooLarge(file: File): boolean {
  return file.size > MAX_DATA_UPLOAD_BYTES;
}

export function suggestSmilesColumn(columns: string[]): string {
  return columns.find((column) => /^smiles$/i.test(column.trim()))
    ?? columns.find((column) => /smiles/i.test(column))
    ?? columns[0]
    ?? "";
}

export function suggestQsarTargetColumn(columns: string[], smilesColumn: string): string {
  return columns.find((column) => /^pchembl\s+value$/i.test(column.trim()))
    ?? columns.find((column) => column !== smilesColumn)
    ?? "";
}
