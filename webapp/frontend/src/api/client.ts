export interface HealthResponse {
  status: "ok";
  database: "ok";
  version: string;
}

export type JobStatus =
  | "queued"
  | "running"
  | "completed"
  | "failed"
  | "cancelling"
  | "cancelled"
  | "interrupted";

export interface Job {
  id: string;
  status: JobStatus;
  created_at: string;
  queued_at: string;
  started_at: string | null;
  finished_at: string | null;
  updated_at: string;
  comment: string;
  parent_smiles: string;
  parent_smiles_canonical: string | null;
  progress: {
    current: number | null;
    total: number | null;
    percent: number | null;
    message: string | null;
  };
  exit_code: number | null;
  error_summary: string | null;
}

export interface JobListResponse {
  items: Job[];
  total: number;
  limit: number;
  offset: number;
}

export interface ParentValidationResponse {
  input_smiles: string;
  canonical_smiles: string;
  molblock: string;
  atom_count: number;
  molecular_formula: string;
}

export interface PubChemLookupResponse {
  cid: number;
  molecular_name: string;
  isomeric_smiles: string;
}

export interface JobLogResponse {
  job_id: string;
  lines: string[];
  truncated: boolean;
}

export interface DescriptorSchema {
  name: string;
  source: string;
  aliases: string[];
  minimum: number | null;
  maximum: number | null;
  missing_count: number;
}

export interface ResultSchemaResponse {
  job_id: string;
  columns: string[];
  descriptors: DescriptorSchema[];
  total_rows: number;
  valid_rows: number;
  invalid_smiles: number;
  warnings: string[];
  rdkit_version: string;
}

export interface ResultMolecule {
  row_number: number;
  smiles: string;
  canonical_smiles: string;
  structure_svg: string;
  values: Record<string, number | null>;
}

export interface MoleculeListResponse {
  job_id: string;
  items: ResultMolecule[];
  total: number;
  limit: number;
  offset: number;
  invalid_smiles: number;
}

export interface NumericSummary {
  count: number;
  missing_count: number;
  minimum: number | null;
  maximum: number | null;
  mean: number | null;
  median: number | null;
  standard_deviation: number | null;
  q1: number | null;
  q3: number | null;
}

export interface StatisticsResponse {
  job_id: string;
  total_rows: number;
  filtered_rows: number;
  invalid_smiles: number;
  summaries: Record<string, NumericSummary>;
  histograms: Record<string, { start: number; end: number; count: number }[]>;
  scatter_x: string | null;
  scatter_y: string | null;
  scatter_sampled: boolean;
  scatter: { row_number: number; x: number; y: number }[];
}

export type ResultQuery = Record<string, string | number | null | undefined>;

interface ErrorEnvelope {
  error?: { code?: string; message?: string };
}

export class ApiError extends Error {
  constructor(
    message: string,
    public readonly status: number,
    public readonly code: string
  ) {
    super(message);
  }
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(path, {
    ...init,
    headers: {
      Accept: "application/json",
      ...(init?.body ? { "Content-Type": "application/json" } : {}),
      ...init?.headers
    }
  });
  if (!response.ok) {
    let body: ErrorEnvelope = {};
    try {
      body = (await response.json()) as ErrorEnvelope;
    } catch {
      // Keep the fallback uniform instead of exposing raw server text.
    }
    throw new ApiError(
      body.error?.message ?? "サーバーで処理できませんでした。",
      response.status,
      body.error?.code ?? "request_failed"
    );
  }
  return (await response.json()) as T;
}

export function getHealth(signal?: AbortSignal): Promise<HealthResponse> {
  return request<HealthResponse>("/api/health", { signal });
}

export function validateParent(
  smiles: string,
  modelKey: string
): Promise<ParentValidationResponse> {
  return request<ParentValidationResponse>("/api/parents/validate", {
    method: "POST",
    body: JSON.stringify({ smiles, model_key: modelKey })
  });
}

export function lookupPubChem(name: string): Promise<PubChemLookupResponse> {
  const parameters = new URLSearchParams({ name });
  return request<PubChemLookupResponse>(`/api/parents/pubchem?${parameters.toString()}`);
}

export function createJob(payload: {
  parent_smiles: string;
  comment: string;
  run_conditions: Record<string, unknown>;
}): Promise<Job> {
  return request<Job>("/api/jobs", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export function listJobs(signal?: AbortSignal): Promise<JobListResponse> {
  return request<JobListResponse>("/api/jobs?limit=100&offset=0", { signal });
}

export function getJob(jobId: string, signal?: AbortSignal): Promise<Job> {
  return request<Job>(`/api/jobs/${encodeURIComponent(jobId)}`, { signal });
}

export function getJobLogs(jobId: string, signal?: AbortSignal): Promise<JobLogResponse> {
  return request<JobLogResponse>(
    `/api/jobs/${encodeURIComponent(jobId)}/logs?tail=200`,
    { signal }
  );
}

export function cancelJob(jobId: string): Promise<Job> {
  return request<Job>(`/api/jobs/${encodeURIComponent(jobId)}/cancel`, {
    method: "POST"
  });
}

export function resultQueryString(query: ResultQuery): string {
  const parameters = new URLSearchParams();
  Object.entries(query).forEach(([key, value]) => {
    if (value !== null && value !== undefined && value !== "") {
      parameters.set(key, String(value));
    }
  });
  const encoded = parameters.toString();
  return encoded ? `?${encoded}` : "";
}

export function getResultSchema(jobId: string, signal?: AbortSignal): Promise<ResultSchemaResponse> {
  return request<ResultSchemaResponse>(`/api/jobs/${encodeURIComponent(jobId)}/results/schema`, { signal });
}

export function getResultMolecules(
  jobId: string,
  query: ResultQuery,
  signal?: AbortSignal
): Promise<MoleculeListResponse> {
  return request<MoleculeListResponse>(
    `/api/jobs/${encodeURIComponent(jobId)}/results/molecules${resultQueryString(query)}`,
    { signal }
  );
}

export function getResultStatistics(
  jobId: string,
  query: ResultQuery,
  signal?: AbortSignal
): Promise<StatisticsResponse> {
  return request<StatisticsResponse>(
    `/api/jobs/${encodeURIComponent(jobId)}/results/statistics${resultQueryString(query)}`,
    { signal }
  );
}

export function resultDownloadUrl(jobId: string, format: "csv" | "sdf", query: ResultQuery): string {
  return `/api/jobs/${encodeURIComponent(jobId)}/results/download.${format}${resultQueryString(query)}`;
}
