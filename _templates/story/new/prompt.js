const toPascalCase = (value) =>
  String(value || '')
    .split(/[^a-zA-Z0-9]+/)
    .filter(Boolean)
    .map((part) => part.charAt(0).toUpperCase() + part.slice(1))
    .join('')

const normalizePosix = (value) => String(value || '').replace(/\\/g, '/')
const toTitleSegment = (segment) => {
  const lower = String(segment || '').toLowerCase()
  if (lower === 'ui') {
    return 'UI'
  }
  if (lower === 'api') {
    return 'API'
  }
  return toPascalCase(segment)
}

module.exports = {
  prompt: ({ inquirer, args }) => {
    const defaults = {
      name: args.name || '',
      dir: normalizePosix(args.dir || 'apps/web/src/components/ui'),
      exportStyle: args.exportStyle || args.export_style || 'named',
      exportName: args.exportName || args.export_name || '',
      titlePrefix: args.titlePrefix || args.title_prefix || '',
    }

    const hasAllRequiredArgs =
      Boolean(defaults.name) &&
      Boolean(defaults.dir) &&
      Boolean(defaults.exportStyle) &&
      (defaults.exportStyle !== 'named' || Boolean(defaults.exportName))

    const normalizeAnswers = (answers) => {
      const name = String(answers.name || defaults.name).trim()
      const dir = normalizePosix(answers.dir || defaults.dir)
      const fileBase = name.replace(/\.(tsx|ts|jsx|js)$/i, '')
      const componentName = toPascalCase(fileBase)
      const exportStyle = answers.exportStyle || defaults.exportStyle
      const exportName =
        exportStyle === 'named'
          ? String(answers.exportName || defaults.exportName || componentName).trim()
          : ''
      const dirSegments = dir.split('/').filter(Boolean)
      const componentsIndex = dirSegments.lastIndexOf('components')
      const categorySegments =
        componentsIndex >= 0 ? dirSegments.slice(componentsIndex + 1) : []
      const category = categorySegments.map(toTitleSegment).join('/')
      const titlePrefix = String(answers.titlePrefix || defaults.titlePrefix).trim()
      const storyTitle = titlePrefix
        ? `${titlePrefix}/${componentName}`
        : `${category || 'Misc'}/${componentName}`
      const importPath = `./${fileBase}`

      return {
        ...answers,
        dir,
        name: fileBase,
        fileBase,
        componentName,
        exportStyle,
        exportName,
        titlePrefix,
        storyTitle,
        importPath,
      }
    }

    if (hasAllRequiredArgs) {
      return Promise.resolve(normalizeAnswers(defaults))
    }

    return inquirer
      .prompt([
        {
          type: 'input',
          name: 'name',
          message: 'Component file base name (without extension)?',
          default: defaults.name,
        },
        {
          type: 'input',
          name: 'dir',
          message: 'Directory for the component and story file?',
          default: defaults.dir,
        },
        {
          type: 'list',
          name: 'exportStyle',
          message: 'Component export style?',
          choices: ['named', 'default'],
          default: defaults.exportStyle,
        },
        {
          type: 'input',
          name: 'exportName',
          message: 'Named export identifier (named export only)?',
          default: defaults.exportName,
          when: (answers) => answers.exportStyle === 'named',
        },
      ])
      .then(normalizeAnswers)
  },
}
