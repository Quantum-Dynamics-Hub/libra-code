# Architecture principles

Read the relevant architecture and compatibility documentation before changing
implementation boundaries. Keep concrete API details and module maps in that
documentation rather than in agent guidance.

- Route complete simulations through the shared engine. Keep method-specific
  numerical composition separate from reusable domain algorithms.
- Treat algorithm variants as options within their parent method when they share
  its lifecycle, rather than introducing separate top-level execution paths.
- Centralize orchestration, clock updates, output cadence, and history ownership.
  Do not introduce duplicate simulation loops or live numerical state containers.
- Separate storage, views, identity, and topology from physics. Higher-level
  algorithms may depend on foundational interfaces, never the reverse.
- Keep imports acyclic and package initializers focused on intentional public
  exports rather than implementations.
- Reuse existing initialization, sampling, model, numerical, and persistence
  infrastructure before introducing new abstractions or parallel implementations.
- Separate field estimation from time evolution; share compatible propagation
  kernels while preserving method-specific equations and numerical conventions.
- Preserve active-index handling, storage lifecycle, array precision, device
  placement, and differentiation support. Make conversions explicit at boundaries.
- Distinguish coherent expansions from independent paths. Preserve their respective
  normalization, representation, units, and observable definitions.
- Preserve supported interfaces; shared orchestration does not require identical
  input or result structures for physically different methods.
- Keep independent regression references separate from production dependencies.
  Verify numerical behavior and lifecycle contracts with proportionate checks.
- Keep tests structurally aligned with source modules and examples consistently
  organized. Update documentation and references when moving implementations.
