---
title: Theory
layout: page
permalink: /theory
nav_order: 4
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

## Theoretical Basis

Basically treating the rays as actual light with wavelengths.



### Spectrial Power Distributions
There we have different types of cone cells LMS.

$$ 
L = \int_0^\infty S(\lambda) \cdot \bar{L}(\lambda) d\lambda\\
M = \int_0^\infty S(\lambda) \cdot \bar{M}(\lambda) d\lambda\\
S = \int_0^\infty S(\lambda) \cdot \bar{S}(\lambda) d\lambda
$$

Each cone cell captures a different section of the visible spectrum.

Note that light is often a combination of different wavelengths.
For example, white light is a combination of all wavelengths in the visible spectrum.
Over the Spectral Power Distribution (SPD) of the light source, or the "intensity" profile of the light at different wavelengths.
We integrate with the spectral sensitivity function of the cone cells to get the response of each cone cell to the light source.
Hence we get the discrete values LMS. For each cone cell, this gives us <L,M, S> a trichromatic value that describes light natural to how the human vision system works.


Unfortunately, that isn't how computer screens work (limited by technological factors).
Instead, we have RGB values.
We can do something very similar.

$$ 
R = \int_0^\infty S(\lambda) \cdot \bar{R}(\lambda) d\lambda\\
G = \int_0^\infty S(\lambda) \cdot \bar{G}(\lambda) d\lambda\\
B = \int_0^\infty S(\lambda) \cdot \bar{B}(\lambda) d\lambda
$$

#### Spectral Profile

- Flat uniform distribution (inaccurate)
- Black body radiation model
- Can be collected by data.

#### Binning
Obviously we can't integrate. We need to discretize.
Binning resolution.

#### Sampling
Obviously, we don't want to actually integrate. We sample instead.
- Binned sampling (sample towards)
- Monte Carlo sampling
- Hero sampling



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

### BSDF function Take 1

We basically do what we do for glass.
We sample using the fresnel coefficient for probability of being reflected or transmitted.
Not only do we need to sample this, but also sample over wavelengths.

Uniform Hemisphere sampling for wo + russian roulette (multiple bounces).
For each sample, we also loop over sampling the wavelengths to compute the RGB

For now:
- uniform SPD
- uniform hemisphere sampling
- uniform wavelength sampling
- one bounce

$$ L_{out} = L_{e}(p, \lambda, \omega_0) + \int_{H^2} f(p, \lambda, \omega_i \rightarrow \omega_o) L_o(p, \lambda, \omega_o) \cos\theta d\omega $$

Discretizing this, we get:
$$ L_{out} = L_{amb?} + \frac{1}{N} \sum \frac{f_r() L_o() \cos\theta}{p(\omega_j)} $$

Where f is amount of reflectance or transmittance of this material intersection.
L is the incoming radiance achieved through that sampled direction.

To achieve f, we sample for the reflectance or transmittance.
Then we also sample for wavelengths.
Once we have that, we now sample across wavelengths:
$$
\begin{bmatrix}
R\\G\\B
\bend{bmatrix}
=
\begin{bmatrix}
\int_0^\infty S(\lambda) \cdot \bar{R}(\lambda) d\lambda\\
\int_0^\infty S(\lambda) \cdot \bar{G}(\lambda) d\lambda\\
\int_0^\infty S(\lambda) \cdot \bar{B}(\lambda) d\lambda
\bend{bmatrix}
$$
Where S is uniform. R is given. 

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

$$ S(\lambda)  = \frac{2 h c^2}{\lambda^5} \left[ exp\left( \frac{hc}{\lambda k_B T} \right) -1 \right]^{-1} $$

### Sampling Techniques

### Redshifting
