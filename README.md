# Microplastics in Terrestrial Mammals

This repository contains the data and scripts required to reproduce all analyses and figures in the paper "Human disturbance drives microplastic abundance in terrestrial wildlife" by Ferraz et al.

## Dependencies

All analyses were conducted in R version 4.3.3 (2024-02-29). The following packages are required:

```r
install.packages(c(
  "ctmm",       # v1.2.1   - movement modelling
  "mgcv",       # v1.9-1   - statistical models
  "weights",       # v1.0.4 - weighted proportions
  "ggplot2",    # v3.5.2   - figures
  "terra",      # v1.8-54  - raster processing
  "sf",         # v1.0-21  - spatial data
  "ggspatial",      # v1.1.9 - scale bars
  "tidyterra",    # v0.7.0 - maps in ggplot
  "rphylopic",    # v1.5.0 - adding animal silhouettes to figures
  "fmsb",      # v0.7.6 - spider plots
  "fields",     # v16.3 - spatial data processing
  "gridExtra",  # v2.3 - combining figures
  "lwgeom"      # v0.2-14 - spatial data processing
))
```

## Repository Structure

```
data/
    mp_data/
        blood_microplastics_polymer_profiles.csv      # Per-sample polymer-specific particle counts
        blood_microplastics_individual_summary.csv    # Per-individual summary data
        blood_microplastics_particle_sizes.csv        # Particle size measurements for each particle
        blank_controls_particle_sizes.csv             # Particle size measurements for each particle in the blank controls
        polymer_classification_key.csv                # Polymer classification reference

    movement_summaries/
        anteater_land_use_summary.csv                 # Land use within anteater home ranges
        anteater_movement_summary.csv                 # Anteater movement and home range statistics
        armadillo_land_use_summary.csv                # Land use within armadillo home ranges
        armadillo_movement_summary.csv                # Armadillo movement and home range statistics
        tapir_land_use_summary.csv                    # Land use within tapir home ranges
        tapir_movement_summary.csv                    # Tapir movement and home range statistics

results/
    anteater_movement/
        akdes/                                        # AKDE home range estimates (.Rda)
        figures/                                      # Model diagnostic figures
        movement_models/                              # Fitted ctmm movement models (.Rda)
    armadillo_movement/
        akdes/
        figures/
        movement_models/
    tapir_movement/
        akdes/
        figures/
        movement_models/

scripts/
    00_custom_project_functions.R
    01_1_anteater_movement_modelling.R
    01_2_armadillo_movement_modelling.R
    01_3_tapir_movement_modelling.R
    02_habitat_use_figures.R
    03_calculate_hfi_change.R
    04_microplastics_data_import.R
    05_interspecific_variation_models.R
    06_microplastics_habitat_use_models.R
    07_microplastics_biological_drivers.R
    08_microplastics_biome_comparisons.R
    09_polymer_composition_analyses.R
    10_sensitivity_analyses.R
    11_figure_1.R
    12_figure_2.R
    13_supplementary_figures.R
    14_supplementary_figure_4.R

figures/
    main_text/                                        # Composite manuscript figures
    main_text_panels/                                 # Individual panels for journal submission
    supplementary/                                    # Supplementary figures
    tracking_data/                                    # Individual animal tracking panels for Figure 2A
```

## Data Availability

The blood microplastics data and movement summary statistics required to reproduce all analyses are available in this repository.

GPS tracking data will be deposited on Movebank prior to publication. The Movebank study ID will be provided here upon acceptance.

The following environmental layers are required to reproduce the analyses but are too large to store on GitHub. They are publicly available at the URLs below:

**Human Footprint Index (main analyses)**
Gassert et al. (2023). An operational approach to near real time global high resolution mapping of the terrestrial Human Footprint. *Frontiers in Remote Sensing*. https://doi.org/10.3389/frsen.2023.1130896
- Data: https://source.coop/repositories/vizzuality/hfp-100/description

**Human Footprint Index (change over time analyses)**
Wildlife Conservation Society Human Footprint dataset.
- Data: https://wcshumanfootprint.org/data-access

**Land cover classification (Brazil)**
MapBiomas Brazil land use and land cover classification.
- Data: https://brasil.mapbiomas.org/en/

## Workflow

Scripts should be run in the following order. Scripts 05 onward require scripts 00–04 to have been run first but are otherwise independent of each other.

| Script | Description |
|--------|-------------|
| `00_custom_project_functions.R` | Custom functions used across all scripts |
| `01_1_anteater_movement_modelling.R` | Fits ctmm movement models and estimates AKDE home ranges for giant anteaters |
| `01_2_armadillo_movement_modelling.R` | Fits ctmm movement models and estimates AKDE home ranges for giant armadillos |
| `01_3_tapir_movement_modelling.R` | Fits ctmm movement models and estimates AKDE home ranges for lowland tapirs |
| `02_habitat_use_figures.R` | Generates individual-level habitat use figures |
| `03_calculate_hfi_change.R` | Calculates change in Human Footprint Index between 2001 and 2020 |
| `04_microplastics_data_import.R` | Imports and prepares microplastics data for analysis |
| `05_interspecific_variation_models.R` | Models interspecific variation in microplastic concentrations and polymer types |
| `06_microplastics_habitat_use_models.R` | Models relationships between HFI and blood microplastic concentrations |
| `07_microplastics_biological_drivers.R` | Models effects of age, sex, and body size on microplastic concentrations |
| `08_microplastics_biome_comparisons.R` | Compares microplastic concentrations across biomes |
| `09_polymer_composition_analyses.R` | Analyses polymer composition across species and biomes |
| `10_sensitivity_analyses.R` | Sensitivity analyses including species-specific models and influential point diagnostics |
| `11_figure_1.R` | Generates Figure 1 |
| `12_figure_2.R` | Generates Figure 2 |
| `13_supplementary_figures.R` | Generates supplementary figures 1-3, 6, and 8 |
| `14_supplementary_figure_5.R` | Generates supplementary figure 5 |

## Citation

Ferraz et al. (in revision). Human disturbance drives microplastic abundance in terrestrial wildlife.

## License

MIT
