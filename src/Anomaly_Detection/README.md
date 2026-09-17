# Anomaly Detection

Core MATLAB implementation of the object-level deformation anomaly
detection module.

This module includes:

- PS-level deformation-model classification using multiple-hypothesis
  testing (MHT);
- aggregation of PS-level labels to building-level dominant classes;
- consolidation into four object-level deformation classes:
  non-anomalous, step-like, velocity-change, and complex;
- cross-track class consistency assessment;
- anomaly screening using matched non-linear classes and persistent
  cross-track mismatches under adequate scatterer sampling.

The main entry point is:

```matlab
run_anomaly_detection.m
