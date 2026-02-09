module.exports = {
  resolve: {
    extensions: ['.ts', '.tsx', '.js', '.jsx', '.mjs', '.cjs', '.json'],
    alias: {
      '@': require('path').resolve(__dirname, 'apps/web/src'),
      '@chem-model/api-client': require('path').resolve(
        __dirname,
        'packages/api-client/src',
      ),
      '@chem-model/shared': require('path').resolve(
        __dirname,
        'packages/shared/src',
      ),
    },
  },
}
