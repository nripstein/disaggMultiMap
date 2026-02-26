# AGHQ Test Context

AGHQ tests are split into two groups:

- [test-aghq-fast.R](/home/nripstein/code/disaggMultiMap/tests/testthat/test-aghq-fast.R): always on, fast argument-validation checks (no AGHQ fitting).
- [test-aghq-gated.R](/home/nripstein/code/disaggMultiMap/tests/testthat/test-aghq-gated.R): heavy fit/predict contract and regression checks.

## Default Behavior

Heavy AGHQ tests are enabled by default.

To opt out of heavy AGHQ tests, set `RUN_AGHQ_TESTS=false`.

## Run AGHQ Tests

Run just heavy AGHQ tests:

```bash
Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-aghq-gated.R', reporter='summary')"
```

Run full suite (includes heavy AGHQ tests by default):

```bash
Rscript -e "devtools::load_all('.'); testthat::test_dir('tests/testthat', reporter='summary')"
```

Run full suite with heavy AGHQ tests disabled:

```bash
RUN_AGHQ_TESTS=false Rscript -e "devtools::load_all('.'); testthat::test_dir('tests/testthat', reporter='summary')"
```

## Current Gated Coverage

1. AGHQ fit object contract on a deterministic one-covariate fixture.
2. AGHQ predict output contract for `predict_iid = FALSE`.
3. `new_data` covariate-name alignment validation.
4. Explicit regression assertion for the known shared-betas/two-covariate fit failure:
`Internal: some slope names not found in theta order.`

## Promotion Criteria (Ungating Heavy AGHQ Tests)

1. No intermittent failures across repeated runs (for example, 10 consecutive runs).
2. Stable fixed-effect name mapping across shared and time-varying beta configurations.
3. Runtime within default test-suite budget.
4. Clear and stable error contracts for unsupported options.
