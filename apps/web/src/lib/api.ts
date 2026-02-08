import { createApiClient } from '@chem-model/api-client'

import { getAuthToken } from './auth'

import type { ApiRequest } from '@chem-model/api-client'
import { requestApi } from '@/server/api'

const request = async <T>(params: ApiRequest): Promise<T> => {
  return (await requestApi({ data: params })) as T
}

const api = createApiClient({
  request,
  getToken: () => getAuthToken() ?? undefined,
})

// API_BASE should be a host root or already end with "/api" (no "/api/v1" style path).
const normalizeApiBase = (base: string) => {
  const trimmed = base.endsWith('/') ? base.slice(0, -1) : base
  return trimmed.endsWith('/api') ? trimmed : `${trimmed}/api`
}

const resolveApiBase = (): string => {
  if (typeof window !== 'undefined' && window.__API_BASE__) {
    return normalizeApiBase(window.__API_BASE__)
  }
  return normalizeApiBase(
    import.meta.env.VITE_API_BASE ?? 'http://localhost:8000',
  )
}

/** QE 入力テキストを構造とパラメータへ解析する。 */
export const parseQeInput = api.parseQeInput
/** QE 入力テキストから永続化済み構造を作成する。 */
export const createStructureFromQe = api.createStructureFromQe
/** ID を指定して永続化済み構造を取得する。 */
export const getStructure = api.getStructure
/** 構造 ID から QE 入力テキストをエクスポートする。 */
export const exportQeInput = api.exportQeInput
/** 構造を CIF テキストとしてエクスポートする。 */
export const exportStructureCif = api.exportStructureCif
/** 小セル/大セル入力から Δ 移植結果を生成する。 */
export const deltaTransplant = api.deltaTransplant
/** supercell オプションに基づいて新規構造を構築する。 */
export const buildSupercell = api.buildSupercell
/** ZPE 入力を解析する。 */
export const parseZpeInput = api.parseZpeInput
/** 新規ユーザーを登録する。 */
export const registerAccount = api.registerAccount
/** ユーザーをログインさせる。 */
export const loginAccount = api.loginAccount
/** 現在の認証ユーザー情報を取得する。 */
export const fetchAuthMe = api.fetchAuthMe
/** 現在のセッションをログアウトする。 */
export const logoutAccount = api.logoutAccount
/** 計算リソース登録用トークンを発行する。 */
export const createEnrollToken = api.createEnrollToken
/** 利用可能なキューターゲット一覧を取得する。 */
export const fetchQueueTargets = api.fetchQueueTargets
/** 現在ユーザーのアクティブなキューターゲットを選択する。 */
export const selectQueueTarget = api.selectQueueTarget
/** ZPE ジョブを作成する。 */
export const createZpeJob = api.createZpeJob
/** ZPE ジョブ状態を ID で取得する。 */
export const fetchZpeStatus = api.fetchZpeStatus
/** ZPE ジョブ結果を ID で取得する。 */
export const fetchZpeResult = api.fetchZpeResult
/** ZPE ジョブ出力ファイルをダウンロードする。 */
export const downloadZpeFile = api.downloadZpeFile

/** 構造 ID から公開表示用 URL を組み立てる。 */
export function structureViewUrl(
  structureId: string,
  params?: {
    format?: 'cif'
  },
): string {
  return `${resolveApiBase()}${api.structureViewPath(structureId, params)}`
}
