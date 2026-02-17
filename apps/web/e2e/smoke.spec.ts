import { expect, test } from '@playwright/test'

test('@smoke loads root page', async ({ page }) => {
  const response = await page.goto('/')
  expect(response?.ok()).toBeTruthy()
  await expect(page).toHaveURL(/\/(?:editor)?\/?$/)
})

test('@smoke loads editor shell', async ({ page }) => {
  await page.goto('/editor')
  await expect(page.getByRole('heading', { name: 'Structures' })).toBeVisible()
  await expect(page.getByLabel('Search tools')).toBeVisible()
})

test('@smoke shows empty workspace hint', async ({ page }) => {
  await page.goto('/editor')
  await expect(page.getByText('No files imported yet.')).toBeVisible()
  await expect(page.getByText('Import a .in file to start.')).toBeVisible()
})
