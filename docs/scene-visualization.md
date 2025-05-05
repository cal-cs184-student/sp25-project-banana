---
layout: page
title: Scene Visualization
nav_order: 4
---

# Thin Film Scene Visualization

This page visualizes the spheres scene created in our spectral renderer.

## Scene Structure

The scene contains multiple spheres with different materials and positions:

- **Sphere 1**: Glass-like thin film sphere (r=0.3) at position (-0.4, 0.3, 0.3)
  - Using thin film material with refraction index 1.33, thickness 300nm
  
- **Sphere 2**: Lambertian sphere (r=0.3) at position (0.4, -0.3, 0.3)
  - Using lambertian material with reflectance (0.5, 0.5, 0.5)

- **Sphere 3**: Smaller lambertian sphere (r=0.15) inside Sphere 1
  - Demonstrates nested geometry for complex optical effects

- **Sphere 4**: Medium lambertian sphere (r=0.25) at position (0.8, 0.7, 0.25)

- **Sphere 5**: Small lambertian sphere (r=0.2) at position (-0.8, -0.6, 0.2)

- **Sphere 6**: Large lambertian sphere (r=0.28) at position (0, 0.8, 0.28)

The camera is positioned to view all spheres without overlaps.

## Visual Representation

