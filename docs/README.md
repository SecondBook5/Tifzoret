# Tifzoret documentation

Tifzoret is a project-agnostic engine for publication-grade bulk RNA-seq
hypothesis analysis. A study is one configuration file that points the engine at
its inputs and declares its contrasts, hypotheses, and figure recipe; the engine
supplies the reusable analysis and rendering.

Start with the study-authoring guide, then reach for the reference docs as
needed.

| Document | Read it for |
|----------|-------------|
| [authoring-a-study.md](authoring-a-study.md) | **Start here** — how to write the config that defines a study and spin up a new paper repo. |
| [cli.md](cli.md) | Complete `tifzoret` command and flag reference — the usage "API". |
| [configuration.md](configuration.md) | Every configuration key, input boundary, profile, and the direction convention. |
| [architecture.md](architecture.md) | How the CLI, Snakemake workflow, and rendering stages fit together. |
| [methods.md](methods.md) | The statistical and analytical methods, written for a methods section. |
| [figures.md](figures.md) | The hypothesis-driven figure system: panel constructors, variants, recipes. |
| [migration.md](migration.md) | Moving an existing project onto the engine and verifying it reproduces prior results. |
| [troubleshooting.md](troubleshooting.md) | Common failures and how to resolve them. |
