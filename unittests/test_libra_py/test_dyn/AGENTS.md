# Testing principles

Read the relevant architecture documentation before adding or revising tests.

- Mirror production module organization and support execution from both the
  repository root and individual test directories.
- Exercise complete workflows through the shared engine, and test reusable
  numerical kernels at their own module boundaries.
- Cover affected storage lifecycle, active-index handling, backend behavior,
  and coherent versus independent-path normalization contracts.
- Keep reference implementations independent and restricted to testing. Do not
  alter reference formulas merely to make a refactored implementation pass.
- Verify numerical outputs using independent references, analytic results,
  finite differences, or physical invariants as appropriate.
- Keep optional compiled comparisons and expensive integration checks explicit;
  ordinary interpreted execution should not acquire optional test dependencies.
- Update example references when examples move. Validate behavior and saved
  outputs rather than relying only on successful imports.
- Run checks proportionate to the change and investigate numerical regressions
  before relaxing tolerances.
