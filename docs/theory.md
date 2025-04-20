---
title: Theory
layout: page
nav_order: 4
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

## Theoretical Basis

Basically treating the rays as actual light with wavelengths.


### Spectrial Power Distributions
### LMS CIE1931

### Thin Films Base Case
#### Refraction - Snell's Law 
To do refract a ray, we'll use Snell's Law.
This is easier for 2D case as we just follow

$$ n_1 \sin{\theta_1} = n_2 \sin{\theta_2} $$

Along the plane created by the incident and refracted ray.

However, let's convert this to 3D vectors. We use a combination of properties:

$$ 
\begin{align*}
\eta_1 \sin{\theta_1} &= \eta_2 \sin{\theta_2}\\
N \times E &= |N| |E| \sin{\theta_1} \eta \\
T &= \sin{(\theta_2)}(N \times \eta) - (N \cdot T) N\\
&= \text{Projection on plane of incidence} + \text{on Normal vector}\\
&= \frac{\eta_1}{\eta_2} \sin{(\theta_1)} (N \times \eta ) - N |N| |T| \cos{(\theta_2)}\\
&= \frac{\eta_1}{\eta_2} (N \times \sin{(\theta_1)} \eta ) - N \cos{(\theta_2)}\\
&= \frac{\eta_1}{\eta_2} (N \times (-N \times E)) - N \sqrt{1 - \sin^2(\theta_2)}\\
&= \frac{\eta_1}{\eta_2} (N \times (-N \times E)) - N \sqrt{1 - \left(\frac{\eta_1}{\eta_2} \sin(\theta_1) \right)^2}\\
\end{align*}
$$

(Fundamentals of Computer Graphics pg. 304)

#### Reflection - Fresnel and Schlick
Rays not only refract, but they reflect in a thin film model.

To find the direction the ray reflects in, we simply take:

$$ R = I − 2 (IN) N $$

Where R is the reflected ray direction, I is incoming ray direction, and N is the surface normal the ray reflects off of.

However, this equation doesn't account for how much of the ray is reflected and how much of it is transmitted and refracted through the material.
For this, we turn to Fresnel equations. Fresnel equations give you the coefficient for the amount of light reflected off a surface, and it's inverse, which is the amount of refracted light.

$$ R(\theta)  ,   (1 - R(\theta)) $$

For the sake of rendering power, we opt to use the Schlick approximation, which says that 

$$ R(\theta) = R_0 + (1 - R_0)(1 - \cos \theta)^5 $$

where $\theta$ is the angle between the viewpoint and the surface normal and $R_0$ is:

$$ R_0 = \left( \frac{n_1 - n_2}{n_1 + n_2} \right)^2 $$

### Soap Bubbles
#### Mesh generation and physics


### Spectral Data Gathering
### White light caustics

### Sampling Techniques

### Redshifting
