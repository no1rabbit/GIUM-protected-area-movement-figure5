# Protected Area Coverage Along Ungulate Movement Trajectories

## Overview

This repository contains the R code used to generate **Figure 5** of the GIUM manuscript.

The script illustrates the workflow used to quantify and visualize the overlap between protected areas (PAs) and the movement trajectories of migratory and nomadic ungulate populations.

## Description

The analysis uses animal movement data and spatial boundaries of protected areas to characterize protection status along movement trajectories.

The workflow includes:

1. **Movement data processing:** Importing and processing animal movement trajectories for migratory and nomadic populations.
2. **Spatial sampling:** Sampling movement trajectories at 1-km intervals.
3. **Protected area overlap:** Classifying sampled locations according to whether they fall within or outside protected areas.
4. **Protection summaries:** Calculating the percentage of sampled movement locations within protected areas, including individual-level and population-level summaries.
5. **Visualization:** Generating proportional ribbon plots showing changes in protected area coverage along seasonal migration trajectories or annual movement cycles of nomadic populations.

The resulting figure distinguishes between seasonal migrants and nomadic populations and illustrates the proportion of movement occurring within and outside protected areas.

## Repository contents

- `figure5_protected_area_ribbon.R` — R script illustrating the analytical workflow and figure generation.

## R packages

The script uses the following R packages:

- `sf`
- `tidyverse`
- `stringr`
- `mapview`
- `ggstream`
- `cowplot`
- `grid`
- `gridExtra`
- `ggspatial`
- `patchwork`
- `ggplotify`

## Data requirements

The script requires animal movement datasets and spatial boundaries of protected areas.

The underlying movement and protected area datasets are not included in this repository.

File paths in the script are provided as placeholders (e.g., `yourpath`) and must be adapted to the user's local directory structure.

**Note:** The code is provided to document the analytical and visualization workflow used in the manuscript. It is not intended as a standalone, directly executable analysis without the corresponding input datasets and configuration.

## Citation

If using or adapting this code, please cite the associated manuscript:

*Full manuscript citation and DOI to be added upon publication.*
