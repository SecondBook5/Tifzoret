# Troubleshooting

- **Unset variable:** set every `${NAME}` referenced by `project.yaml`; validation names unresolved variables.
- **Offline resource error:** run once with network access or point `resources.cache` at a populated cache.
- **No gene sets remain:** verify symbol species/case, collection selection, and `min_size`/`max_size`.
- **Coefficient cannot be resolved:** confirm the contrast factor is in the design and the denominator is a real selected level.
- **Archive error:** BAM member paths are relative to `member_root`; absolute paths, links, and `..` traversal are rejected.
- **Small-sample warning:** the module ran as requested, but its result should remain exploratory; the warning is retained in the manifest.
- **Migration mismatch:** run `bulk-rna prepare` first, verify exact counts, then compare DE and displayed-data tables with `bulk-rna verify`.

