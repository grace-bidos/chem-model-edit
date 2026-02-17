import { HttpResponse, http } from 'msw'

export const defaultHandlers = [
  http.get('http://localhost:8000/api/health', () =>
    HttpResponse.json({ ok: true }),
  ),
]
