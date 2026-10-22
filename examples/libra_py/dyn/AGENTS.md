# Example principles

Read the relevant architecture and example documentation before making changes.

- Run complete simulations through the shared engine and reuse production
  initialization, models, numerical kernels, and persistence infrastructure.
- Organize examples consistently by method. Prefer a self-contained entry point
  combining simulation and plotting, with optional stage selection when useful.
- Construct inputs explicitly and explain their physical meaning with comments.
- Keep comparable command-line options, summaries, output properties, and plot
  conventions consistent across methods without obscuring physical differences.
- Label units, representations, population definitions, and normalization clearly.
- Produce individual plot artifacts without requiring interactive display.
- Do not duplicate production physics or depend on test-only reference code.
- Update documentation and referring tests when reorganizing examples.
- Verify short runs and their saved numerical records and plots. Keep generated
  validation output separate from tracked source files.
