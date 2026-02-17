/** @type {import('dependency-cruiser').IConfiguration} */
module.exports = {
  forbidden: [
    {
      name: 'no-circular',
      severity: 'error',
      comment: 'Avoid circular dependencies.',
      from: {},
      to: { circular: true },
    },
    {
      name: 'no-ui-to-features',
      comment:
        'Base UI components must not depend on feature layer modules.',
      severity: 'error',
      from: { path: '^apps/web/src/components/ui' },
      to: { path: '^apps/web/src/features' },
    },
    {
      name: 'no-server-to-ui-or-features',
      comment:
        'Server layer should stay independent from UI/feature presentation modules.',
      severity: 'error',
      from: { path: '^apps/web/src/server' },
      to: { path: '^apps/web/src/(components|features)' },
    },
    {
      name: 'no-feature-to-routes',
      comment: 'Feature layer should not import route definitions.',
      severity: 'error',
      from: { path: '^apps/web/src/features' },
      to: { path: '^apps/web/src/routes' },
    },
    {
      name: 'no-orphans',
      comment:
        'Highlight orphaned modules to keep dead code visible. Promoted to warn for now.',
      severity: 'warn',
      from: {
        orphan: true,
        pathNot:
          '^apps/web/src/(start\\.tsx|router\\.tsx|routeTree\\.gen\\.ts|worker\\.ts|types/.*\\.d\\.ts)$',
      },
      to: {},
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
