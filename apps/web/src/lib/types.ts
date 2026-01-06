export type Atom = {
  symbol: string
  x: number
  y: number
  z: number
}

export type Structure = {
  atoms: Array<Atom>
  lattice?: Lattice | null
}

export type Vector3 = {
  x: number
  y: number
  z: number
}

export type Lattice = {
  a: Vector3
  b: Vector3
  c: Vector3
}

export type LatticeParams = {
  a: number
  b: number
  c: number
  alpha: number
  beta: number
  gamma: number
}

export type SupercellMeta = {
  na: number
  nb: number
  layers: number
  overlapCount?: number
}
