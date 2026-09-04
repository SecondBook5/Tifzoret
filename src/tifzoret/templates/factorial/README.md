# Factorial Design Template

This template demonstrates a **2×2 factorial design** with treatment and genotype factors, ideal for studying interactions where the effect of one factor depends on the level of another.

## Design

- **Factors:** `treatment` (control, treated) × `genotype` (WT, mutant)
- **Groups:** 4 crossed groups with 3 biological replicates each (12 samples total)
- **Question:** Does the treatment effect differ between WT and mutant genotypes?

## Estimands

The `project.yaml` defines a single family with four estimands:

1. **main_treatment** - Treatment effect averaged over both genotypes
2. **main_genotype** - Genotype effect averaged over both treatment conditions
3. **interaction** - Treatment × Genotype interaction (does treatment work differently in mutant?)
4. **simple_treated** - Simple effect of genotype in treated samples only

## Interpretation

- **Main effects** answer: "Does treatment matter (ignoring genotype)?" and "Does genotype matter (ignoring treatment)?"
- **Interaction** answers: "Does the treatment effect depend on genotype?" A significant interaction means the treatment works differently in WT vs. mutant
- **Simple effect** examines genotype differences only in the treated condition

## Data

The synthetic counts include:
- 100 genes total
- 20 genes respond to treatment in WT but not mutant (interaction genes)
- 15 genes respond to treatment in both genotypes (main treatment effect)
- 10 genes differ by genotype regardless of treatment (main genotype effect)
- 55 null genes

## Usage

```bash
# Initialize a project from this template
tifzoret init --input factorial my_factorial_study

# Validate configuration
tifzoret validate my_factorial_study/project.yaml

# Dry-run to see the DAG
tifzoret dry-run my_factorial_study/project.yaml

# Execute the analysis
tifzoret run my_factorial_study/project.yaml
```

Results will include DE tables for all four estimands plus term tests for the two main effects, the interaction, and any-effect omnibus test.
