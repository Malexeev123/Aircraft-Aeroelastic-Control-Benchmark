# Generated benchmark results

Formal and custom benchmark runs create versioned output directories here:

```text
results/<case>/<run-id>/
```

Each run can include its manifest, configuration, numerical data, metrics,
logs, checkpoints, captions, and diagnostic/publication plots. Ordinary run
contents are intentionally ignored by Git.

Compact accepted comparisons that support a regression or publication figure
may be retained under `results/validation/`. These selected files have an
explicit data contract and are included by the release inventory; timestamped
reproduction runs beneath them remain local output.
