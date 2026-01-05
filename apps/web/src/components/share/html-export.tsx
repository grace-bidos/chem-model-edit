import type { Atom } from '../../lib/types'
import { atomsToPdb } from '../../lib/pdb'

type ShareStructure = {
  id: string
  name: string
  atoms: Atom[]
  opacity: number
  isVisible: boolean
}

type ShareModel = {
  name: string
  pdb: string
  opacity: number
}

function buildHtml(models: ShareModel[]): string {
  const serializedModels = JSON.stringify(models)
  return `<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Chem Model Share</title>
  <script src="https://unpkg.com/ngl@2/dist/ngl.js"></script>
  <style>
    html, body { margin: 0; height: 100%; background: #0b1120; }
    #viewport { width: 100vw; height: 100vh; }
    .badge { position: absolute; top: 16px; left: 16px; padding: 6px 12px; background: rgba(15,23,42,0.7); color: #e2e8f0; border-radius: 999px; font-family: system-ui, sans-serif; font-size: 12px; }
  </style>
</head>
<body>
  <div class="badge">Chem Model Share (NGL)</div>
  <div id="viewport"></div>
  <script>
    const stage = new NGL.Stage('viewport', { backgroundColor: '#0b1120' });
    const models = ${serializedModels};
    const loaders = models.map((model) => {
      const blob = new Blob([model.pdb], { type: 'text/plain' });
      return stage.loadFile(blob, { ext: 'pdb', name: model.name }).then((comp) => {
        comp.addRepresentation('ball+stick', { opacity: model.opacity });
      });
    });
    Promise.all(loaders).then(() => {
      stage.autoView();
    });
  </script>
</body>
</html>`
}

export function downloadShareHtml(
  structures: ShareStructure[],
  options: { activeId: string; overlayEnabled: boolean; filename?: string },
) {
  const targets = options.overlayEnabled
    ? structures
    : structures.filter((structure) => structure.id === options.activeId)
  const models = targets
    .filter((structure) => structure.isVisible && structure.atoms.length > 0)
    .map((structure) => ({
      name: structure.name,
      pdb: atomsToPdb(structure.atoms),
      opacity: Math.min(1, Math.max(0, structure.opacity)),
    }))
  if (models.length === 0) {
    return
  }
  const html = buildHtml(models)
  const blob = new Blob([html], { type: 'text/html' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = options.filename ?? 'chem-model-share.html'
  link.click()
  URL.revokeObjectURL(url)
}
