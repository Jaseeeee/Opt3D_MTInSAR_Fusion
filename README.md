# Opt3D_MTInSAR_Fusion
MATLAB implementation of an Opt3D-assisted fusion framework for ascending and descending MT-InSAR data alignment and object-oriented deformation anomaly detection.

## Overview
This repository provides the core implementation of the proposed
framework, including:

- edge-structure based global co-registration;
- SAR-directional scatterer-to-building association;
- LoD1 building height correction;
- dual-track MT-InSAR fusion;
- object-oriented deformation anomaly detection using MHT.

## Framework
The workflow consists of four main modules:

1. Global co-registration
2. Bidirectional local correction
3. Dual-track fusion
4. Object-oriented deformation anomaly detection

## Data Availability
The TerraSAR-X data and GlobalBuildingAtlas (GBA) building models
used in the study are not included in this repository due to data
access restrictions.

Users should provide their own MT-InSAR products and LoD1 building
models as inputs.

## Requirements
MATLAB R2025a or later

## Citation
If you use this code, please cite

## Code Availability
This repository currently provides the core MATLAB implementation of
the single-track Opt3D-assisted correction module, including global
co-registration, SAR-directional scatterer-to-building association,
LoD1 height correction, and PS repositioning.

The dual-track fusion and object-level deformation anomaly detection
modules will be released upon publication of the associated paper.
