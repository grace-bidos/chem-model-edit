/** @type {import('dependency-cruiser').IConfiguration} */
module.exports = {
  forbidden: [
    {
      name: 'no-ui-to-features',
      comment:
        'Base UI components must not depend on feature layer modules.',
      severity: 'error',
      from: { path: '^apps/web/src/components/ui' },
      to: { path: '^apps/web/src/features' },
    },
  ],
  options: {
    webpackConfig: {
      fileName: './dependency-cruiser.webpack.cjs',
    },
    doNotFollow: {
      path: 'node_modules',
    },
    includeOnly: '^apps/web/src',
    exclude: {
      path:
        '\\.(stories|test|a11y\\.test|fastcheck\\.test)\\.(ts|tsx)$|/storybook-static/|/dist/',
    },
    reporterOptions: {
      dot: {
        collapsePattern:
          '^(apps/web/src/(components/ui|features/editor-v2/components)|apps/web/src/[^/]+)',
      },
    },
  },
}
