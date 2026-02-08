import { clsx } from 'clsx'
import { twMerge } from 'tailwind-merge'
import type { ClassValue } from 'clsx'

/** 条件付きクラス名を結合し、Tailwind の競合を解決する。 */
export function cn(...inputs: Array<ClassValue>) {
  return twMerge(clsx(inputs))
}
