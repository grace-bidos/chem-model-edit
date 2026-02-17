import { expect, test } from '@playwright/test'

test('loads root page', async ({ page }) => {
  const response = await page.goto('/')
  expect(response?.ok()).toBeTruthy()
  await expect(page).toHaveURL(/\/(?:editor)?\/?$/)
})
