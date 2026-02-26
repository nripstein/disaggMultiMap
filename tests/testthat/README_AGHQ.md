# AGHQ Test Context

This repository includes a gated AGHQ smoke test:

- [test-aghq-gated.R](/home/nripstein/code/disaggMultiMap/tests/testthat/test-aghq-gated.R)

It is skipped by default so AGHQ instability does not block the main test suite.

## Run AGHQ Tests

```bash
RUN_AGHQ_TESTS=true Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-aghq-gated.R', reporter='summary')"
```

Or run full suite with AGHQ enabled:

```bash
RUN_AGHQ_TESTS=true Rscript -e "devtools::load_all('.'); testthat::test_dir('tests/testthat', reporter='summary')"
```

## Why Tests Are Gated

AGHQ fit/predict behavior is currently less stable than the TMB path in this codebase (naming/order and optimizer-path issues can trigger run-specific failures).  
The gate keeps default CI/dev feedback fast and reliable while AGHQ support is being hardened.

## Next Steps For AGHQ Test Development

1. Stabilize a canonical AGHQ fit fixture (small deterministic data, known optimizer settings).
2. Add robust fit contract checks:
- class and required list members
- `model_setup` fields (`family`, `link`, `time_varying_betas`)
3. Add predict smoke checks with `predict_iid = FALSE`:
- output class
- mean and CI layer counts matching `n_times`
4. Add explicit checks for known limitations:
- field uncertainty currently mode-conditional
- IID prediction unsupported in AGHQ path
5. Ungate AGHQ tests only after repeated green runs in local and CI-like environments.

## Promotion Criteria (Ungating)

Use all of the following before moving AGHQ tests into default pass criteria:

1. No intermittent failures across repeated runs (for example, 10 consecutive runs).
2. Stable parameter-name mapping for draw extraction across supported model setups.
3. Runtime acceptable for default suite budget.
4. Clear error messages for unsupported AGHQ options.
