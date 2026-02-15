# GRA-18 BFF Aggregation Boundary (First Slice)

## Scope

This document defines the first implementation slice for `GRA-18`.

Implemented in this slice:

- one BFF server-function path for job status reads
- DTO normalization from upstream payloads to a single web DTO
- auth guard enforcement at the BFF boundary

Out of scope for this slice:

- Convex runtime query integration
- AiiDA detail endpoint fan-out
- write-path aggregation

## Ownership Boundary

- `Convex` owns product-facing projection fields (job list/read model metadata).
- `AiiDA/Postgres` owns execution internals and provenance details.
- Web BFF server functions own:
  - auth guard at web boundary,
  - upstream orchestration (fan-out in later slices),
  - DTO normalization for stable web-consumer contracts.

## First Path Contract

Server function: `fetchAggregatedZpeJobStatus`

Input contract:

```ts
{
  jobId: string
  token: string | null
}
```

Rules:

- `token` is required for this protected path.
- missing/empty token is rejected before upstream calls.

Output contract (web DTO):

```ts
{
  status: 'queued' | 'started' | 'finished' | 'failed'
  detail: string | null
  updated_at: string | null
}
```

## Upstream DTOs Accepted by Normalizer

Legacy FastAPI status shape:

```ts
{
  status: 'queued' | 'started' | 'finished' | 'failed'
  detail?: string | null
  updated_at?: string | null
}
```

Convex-style projection shape (forward-compatible normalization):

```ts
{
  jobId: string
  status: 'queued' | 'started' | 'finished' | 'failed'
  detail?: string | null
  updatedAt?: string | null
  eventTime?: string | null
}
```

Normalization behavior:

- always emit `updated_at` in snake_case for current web client compatibility
- map `updatedAt` or `eventTime` to `updated_at` for Convex-style payloads
- reject payloads that do not match either accepted upstream shape

## Error Mapping

BFF path maps/propagates errors with these semantics:

- missing token -> `unauthorized` (`Authentication required`)
- upstream non-2xx -> pass through API error message from standard backend error payload
- invalid upstream DTO shape -> internal BFF error (`Unsupported status payload from upstream`)

## TODO Boundaries

- Add Convex read query for projection metadata and merge with AiiDA/adapter status read.
- Add explicit BFF error type with stable machine-readable code for all normalization failures.
- Expand contract to include source metadata once UI consumes it.
