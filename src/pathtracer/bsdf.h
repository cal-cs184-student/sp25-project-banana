#ifndef CGL_STATICSCENE_BSDF_H
#define CGL_STATICSCENE_BSDF_H

#include "CGL/CGL.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"

#include "pathtracer/sampler.h"
#include "util/image.h"

#include <algorithm>
#include <vector>

namespace CGL {

// Helper math functions. Assume all vectors are in unit hemisphere //
#define PLANCK_CONSTANT 6.62607015e-34
#define SPEED_OF_LIGHT 2.99792458e8
#define BOLTZMANN_CONSTANT 1.38064852e-23

inline double clamp (double n, double lower, double upper) {
  return std::max(lower, std::min(n, upper));
}

inline double cos_theta(const Vector3D w) {
  return w.z;
}

inline double abs_cos_theta(const Vector3D w) {
  return fabs(w.z);
}

inline double sin_theta2(const Vector3D w) {
  return fmax(0.0, 1.0 - cos_theta(w) * cos_theta(w));
}

inline double sin_theta(const Vector3D w) {
  return sqrt(sin_theta2(w));
}

inline double cos_phi(const Vector3D w) {
  double sinTheta = sin_theta(w);
  if (sinTheta == 0.0) return 1.0;
  return clamp(w.x / sinTheta, -1.0, 1.0);
}

inline double sin_phi(const Vector3D w) {
  double sinTheta = sin_theta(w);
  if (sinTheta) return 0.0;
  return clamp(w.y / sinTheta, -1.0, 1.0);
}

void make_coord_space(Matrix3x3& o2w, const Vector3D n);

/**
 * Interface for BSDFs.
 * BSDFs (Bidirectional Scattering Distribution Functions)
 * describe the ratio of incoming light scattered from
 * incident direction to outgoing direction.
 * Scene objects are initialized with a BSDF subclass, used
 * to represent the object's material and associated properties.
 */
class BSDF {
 public:

  /**
   * Evaluate BSDF.
   * Given incident light direction wi and outgoing light direction wo. Note
   * that both wi and wo are defined in the local coordinate system at the
   * point of intersection.
   * \param wo outgoing light direction in local space of point of intersection
   * \param wi incident light direction in local space of point of intersection
   * \return reflectance in the given incident/outgoing directions
   */
  virtual Vector3D f (const Vector3D wo, const Vector3D wi) = 0;

  /**
   * Evaluate BSDF.
   * Given the outgoing light direction wo, samplea incident light
   * direction and store it in wi. Store the pdf of the sampled direction in pdf.
   * Again, note that wo and wi should both be defined in the local coordinate
   * system at the point of intersection.
   * \param wo outgoing light direction in local space of point of intersection
   * \param wi address to store incident light direction
   * \param pdf address to store the pdf of the sampled incident direction
   * \return reflectance in the output incident and given outgoing directions
   */
  virtual Vector3D sample_f (const Vector3D wo, Vector3D* wi, double* pdf) = 0;

  /**
   * Get the emission value of the surface material. For non-emitting surfaces
   * this would be a zero energy Vector3D.
   * \return emission Vector3D of the surface material
   */
  virtual Vector3D get_emission () const = 0;

  /**
   * If the BSDF is a delta distribution. Materials that are perfectly specular,
   * (e.g. water, glass, mirror) only scatter light from a single incident angle
   * to a single outgoing angle. These BSDFs are best described with alpha
   * distributions that are zero except for the single direction where light is
   * scattered.
   */
  virtual bool is_delta() const = 0;

  virtual void render_debugger_node() {};

  /**
   * Reflection helper
   */
  virtual void reflect(const Vector3D wo, Vector3D* wi);

  /**
   * Refraction helper
   */
  virtual bool refract(const Vector3D wo, Vector3D* wi, double ior);
  virtual bool refract(const Vector3D wo, Vector3D* wi, double ior_o, double ior_i);

  const HDRImageBuffer* reflectanceMap;
  const HDRImageBuffer* normalMap;

}; // class BSDF

/**
 * Diffuse BSDF.
 */
class DiffuseBSDF : public BSDF {
 public:

  /**
   * DiffuseBSDFs are constructed with a Vector3D as input,
   * which is stored into the member variable `reflectance`.
   */
  DiffuseBSDF(const Vector3D a) : reflectance(a) { }

  Vector3D f(const Vector3D wo, const Vector3D wi);
  Vector3D sample_f(const Vector3D wo, Vector3D* wi, double* pdf);
  Vector3D get_emission() const { return Vector3D(); }
  bool is_delta() const { return false; }

  void render_debugger_node();

private:
  /*
   * Reflectance is also commonly called the "albedo" of a surface,
   * which ranges from [0,1] in RGB, representing a range of
   * total absorption(0) vs. total reflection(1) per color channel.
   */
  Vector3D reflectance;
  /*
   * A sampler object that can be used to obtain
   * a random Vector3D sampled according to a 
   * cosine-weighted hemisphere distribution.
   * See pathtracer/sampler.cpp.
   */
  CosineWeightedHemisphereSampler3D sampler;

}; // class DiffuseBSDF

/**
 * Microfacet BSDF.
 */

class MicrofacetBSDF : public BSDF {
public:

  MicrofacetBSDF(const Vector3D eta, const Vector3D k, double alpha)
    : eta(eta), k(k), alpha(alpha) { }

  double getTheta(const Vector3D w) {
    return acos(clamp(w.z, -1.0 + 1e-5, 1.0 - 1e-5));
  }

  double Lambda(const Vector3D w) {
    double theta = getTheta(w);
    double a = 1.0 / (alpha * tan(theta));
    return 0.5 * (erf(a) - 1.0 + exp(-a * a) / (a * PI));
  }

  Vector3D F(const Vector3D wi);

  double G(const Vector3D wo, const Vector3D wi);

  double D(const Vector3D h);

  Vector3D f(const Vector3D wo, const Vector3D wi);
  Vector3D sample_f(const Vector3D wo, Vector3D* wi, double* pdf);
  Vector3D get_emission() const { return Vector3D(); }
  bool is_delta() const { return false; }

  void render_debugger_node();

private:
  Vector3D eta, k;
  double alpha;
  UniformGridSampler2D sampler;
  CosineWeightedHemisphereSampler3D cosineHemisphereSampler;
}; // class MicrofacetBSDF

/**
 * Mirror BSDF
 */
class MirrorBSDF : public BSDF {
 public:

  MirrorBSDF(const Vector3D reflectance) : reflectance(reflectance) { }

  Vector3D f(const Vector3D wo, const Vector3D wi);
  Vector3D sample_f(const Vector3D wo, Vector3D* wi, double* pdf);
  Vector3D get_emission() const { return Vector3D(); }
  bool is_delta() const { return true; }

  void render_debugger_node();

private:

  double roughness;
  Vector3D reflectance;

}; // class MirrorBSDF*/

/**
 * Refraction BSDF.
 */
class RefractionBSDF : public BSDF {
 public:

  RefractionBSDF(const Vector3D transmittance, double roughness, double ior)
    : transmittance(transmittance), roughness(roughness), ior(ior) { }

  Vector3D f(const Vector3D wo, const Vector3D wi);
  Vector3D sample_f(const Vector3D wo, Vector3D* wi, double* pdf);
  Vector3D get_emission() const { return Vector3D(); }
  bool is_delta() const { return true; }

  void render_debugger_node();

 private:

  double ior;
  double roughness;
  Vector3D transmittance;

}; // class RefractionBSDF

/**
 * Glass BSDF.
 */
class GlassBSDF : public BSDF {
 public:

  GlassBSDF(const Vector3D transmittance, const Vector3D reflectance,
            double roughness, double ior) :
    transmittance(transmittance), reflectance(reflectance),
    roughness(roughness), ior(ior) { }

  Vector3D f(const Vector3D wo, const Vector3D wi);
  Vector3D sample_f(const Vector3D wo, Vector3D* wi, double* pdf);
  Vector3D get_emission() const { return Vector3D(); }
  bool is_delta() const { return true; }

  void render_debugger_node();

 private:

  double ior;
  double roughness;
  Vector3D reflectance;
  Vector3D transmittance;

}; // class GlassBSDF

/**
 * Emission BSDF.
 */
class EmissionBSDF : public BSDF {
 public:

  EmissionBSDF(const Vector3D radiance) : radiance(radiance) { }

  Vector3D f(const Vector3D wo, const Vector3D wi);
  Vector3D sample_f(const Vector3D wo, Vector3D* wi, double* pdf);
  Vector3D get_emission() const { return radiance; }
  bool is_delta() const { return false; }

  void render_debugger_node();

 private:

  Vector3D radiance;
  CosineWeightedHemisphereSampler3D sampler;

}; // class EmissionBSDF

/**
 * Spectral BSDF.
 */
class SpectralBSDF : public BSDF {
 public:

  SpectralBSDF(const Vector3D transmittance, const Vector3D reflectance, double ior) :
    transmittance(transmittance), reflectance(reflectance), ior(ior) { }

  Vector3D f(const Vector3D wo, const Vector3D wi);
  Vector3D sample_f(const Vector3D wo, Vector3D* wi, double* pdf);
  Vector3D get_emission() const { return Vector3D(); };
  bool is_delta() const { return false; };

  double uniform_spd(double lambda) const { return 1.0 / (700 - 380); };
  double black_body_spd(double lambda);
  double custom_spd(double lambda);
  Vector3D sample_lambda();
  Vector3D to_xyz(double lambda);

  std::vector<double> spd_values = {
    49.98, 52.31, 54.65, 57.00, 59.34, 61.69, 64.03, 66.38, 68.73, 71.07,
    73.42, 75.77, 78.12, 80.47, 82.81, 85.16, 87.51, 89.85, 92.20, 94.55,
    96.89, 99.24, 101.59, 103.94, 106.29, 108.63, 110.98, 113.33, 115.68, 118.03,
    120.37, 122.72, 125.07, 127.42, 129.77, 132.11, 134.46, 136.81, 139.16, 141.51,
    143.85, 146.20, 148.55, 150.90, 153.25, 155.60, 157.94, 160.29, 162.64, 164.99,
    167.34, 169.68, 172.03, 174.38, 176.73, 179.08, 181.43, 183.77, 186.12, 188.47,
    190.82, 193.17, 195.51, 197.86, 200.21, 202.56, 204.91, 207.26, 209.60, 211.95,
    214.30, 216.65, 219.00, 221.34, 223.69, 226.04, 228.39, 230.74, 233.08, 235.43,
    237.78, 236.62, 235.45, 234.29, 233.13, 231.96, 230.80, 229.64, 228.47, 
    227.31, 226.15
}; // must be ordered

std::vector<double> spd_lambda = {
  380, 385, 390, 395, 400, 405, 410, 415, 420, 425,
  430, 435, 440, 445, 450, 455, 460, 465, 470, 475,
  480, 485, 490, 495, 500, 505, 510, 515, 520, 525,
  530, 535, 540, 545, 550, 555, 560, 565, 570, 575,
  580, 585, 590, 595, 600, 605, 610, 615, 620, 625,
  630, 635, 640, 645, 650, 655, 660, 665, 670, 675,
  680, 685, 690, 695, 700, 705, 710, 715, 720, 725,
  730, 735, 740, 745, 750, 755, 760, 765, 770, 775,
  780, 785, 790, 795, 800, 805, 810, 815, 820, 
  825, 830
};

  void render_debugger_node();

 private:

  double ior;
  Vector3D reflectance;
  Vector3D transmittance;

}; // class SpectralBSDF

}  // namespace CGL

#endif  // CGL_STATICSCENE_BSDF_H
