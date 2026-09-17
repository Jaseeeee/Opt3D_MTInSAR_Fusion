# Dual-track Fusion

Core MATLAB implementation of the dual-track fusion module for
ascending- and descending-track MT-InSAR data.

This module includes:

- identification of the geographic overlap between the two tracks;
- association of building objects observed by both tracks;
- fusion of single-track LoD1 building-height estimates;
- relative LOS-velocity reference alignment between tracks;
- track-specific PS repositioning using the fused building geometry.

The main entry point is:

```matlab
run_dual_track_fusion.m
