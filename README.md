# SeedMass

This repository contains data analysis scripts for manuscript

## Q1: Spatial Patterns of Seed Mass

The script `code/Q1Q5.R` looks at general phylogenetic patterns in mean seed mass and seed mass variation

Main models:
  - Calculates Pagel’s λ for mean seed mass, seed mass variation
  - Quantify phylogenetic signal in seed mass response to climate change

## Q2: Spatial Patterns of Seed Mass

The script `code/Q2.R` looks at general spatial patterns in seed mass across the western US

Main models:
  - Fits GLM and GLMM models to test how seed mass varies by plant growth form & duration (annual/perennial)
  - Assesses geographic gradients in seed mass using latitude, longitude, and elevation
  - IDW analysis across the entire western US

## Q3: Life history and environmental predictors of seed mass across species and clades

The script `code/Q3ab.R` looks at environmental predictors of seed mass

Main models:
  - Fits SAR models for species collected at least 3 times
  - Fits GLM and GLMM for well collected clades
  - Fits GLM to test how seed mass varied over time
  - Uses Anscombe-Glynn test to assess distribution of coefficients from GLM for years

## Q4: Temporal, environmental, climate change predictors of seed mass within species

The script `code/Q4.R` looks at environmental predictors of seed mass

Main models:
  - Fits SAR models for species collected at least 20 times
