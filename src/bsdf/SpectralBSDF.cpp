#include "CGL/vector3D.h"
using CGL::Vector3D;

Vector3D SpectralBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
    // (A) sample one hero wavelength and build its small-spectrum
    const double λhero = 380.0 + random_uniform() * (780.0 - 380.0);
    const auto λs = hero_sampler(λhero);
    double Ravg = 0.0, invN = 1.0/λs.size();
    Vector3D xyzR(0), xyzT(0);

    for (double λ : λs) {
        double Rλ = airy_reflectance(wo, λ, thickness);
        Ravg += Rλ;
        xyzR += CGL::ColorSpace::wave2xyz_table(λ) * Rλ;
        xyzT += CGL::ColorSpace::wave2xyz_table(λ) * (1.0 - Rλ);
    }
    Ravg *= invN;
    Vector3D rgbR = CGL::ColorSpace::xyz2rgb(xyzR * invN);
    Vector3D rgbT = CGL::ColorSpace::xyz2rgb(xyzT * invN);

    // (B) Russian roulette on film: reflect vs transmit
    if (coin_flip(Ravg)) {
        Vector3D new_wi;
        reflect(wo, &new_wi);
        *pdf = Ravg;
        *wi  = new_wi;
        return reflectance * rgbR / (abs_cos_theta(new_wi) * (*pdf));
    }

    // (C1) refract into film
    double nbar = 0;
    for (double λ : λs) nbar += film_ior + 0.003/(λ*λ);
    nbar *= invN;
    Vector3D new_wi;
    if (!refract(wo, &new_wi, nbar)) {
        *pdf = 1.0; 
        return Vector3D(0.0);
    }
    double Tavg = 1.0 - Ravg;

    // (C2) sample base BSDF
    double base_pdf = 1.0;
    Vector3D base_f(1.0);
    if (base_bsdf) {
        Vector3D wi_base;
        base_f = base_bsdf->sample_f(new_wi, &wi_base, &base_pdf);
        new_wi = wi_base;
    }

    // (C3) combine pdf & throughput
    *pdf = Tavg * base_pdf;
    *wi  = new_wi;
    return (transmittance * rgbT * base_f) / (*pdf);
}