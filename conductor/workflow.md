# Workflow

## TDD Policy
**Strict** — tests must be written before implementation code. No feature or bug fix is considered complete without corresponding unit tests passing. New algorithms must have tests covering correctness, edge cases, and numerical accuracy.

## Commit Strategy
**Conventional Commits** format is required:
- `feat:` — new feature or algorithm
- `fix:` — bug fix
- `chore:` — build system, CI, tooling, deps
- `docs:` — documentation only
- `test:` — tests only
- `refactor:` — code restructuring without behavior change
- `perf:` — performance improvement

Example: `feat: add LSCM solver with sparse Cholesky backend`

## Code Review Policy
Required for non-trivial changes. Trivial changes (typo fixes, comment updates, minor doc edits) may be self-reviewed. All algorithm implementations and public API changes require review.

## Verification Checkpoints
Manual verification is required **after each phase completion**. At each checkpoint:
1. All tests pass (`ctest`)
2. Examples build and run correctly
3. Public API is reviewed for consistency and usability
4. Documentation builds without warnings (if docs are affected)

## Task Lifecycle
1. **Spec** — write spec.md describing the task, acceptance criteria, and test cases
2. **Tests** — write failing tests
3. **Implementation** — make tests pass
4. **Verify** — run full test suite, check examples, review API
5. **Done** — merge after phase checkpoint sign-off
