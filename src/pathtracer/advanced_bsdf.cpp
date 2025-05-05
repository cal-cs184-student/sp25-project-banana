#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>
#include <complex>
#include <fstream>
#include <mutex>

#include "application/visual_debugger.h"
#include "color_vision.cpp"

using std::max;
using std::min;
using std::swap;
using std::cout;
using std::endl;

namespace CGL {

// Create a debug logging system
static std::ofstream debug_log_file;
static std::mutex debug_log_mutex;

// Initialize debug logging
static void init_debug_log() {
    static bool initialized = false;
    if (!initialized) {
        debug_log_file.open("/Users/ashvinverma/Documents/Courses/cs184/sp25-project-banana/thin_film_debug.log", std::ios::out | std::ios::trunc);
        initialized = true;
    }
}

// Debug logging macro
#define DEBUG_LOG(message) \
    do { \
        static int debug_counter = 0; \
        if (debug_counter++ % 10000 == 0) { \
            std::lock_guard<std::mutex> lock(debug_log_mutex); \
            init_debug_log(); \
            if (debug_log_file.is_open()) { \
                debug_log_file << message << std::endl; \
                debug_log_file.flush(); \
            } \
        } \
    } while (0)

// Define this to enable forced reflectance reduction for better visual appearance
// Comment out to see raw thin film reflectance values
//#define REDUCE_THINFILM_REFLECTANCE

// Returns the RGB spectral primaries basis for a given wavelength (nm)
static Vector3D spectral_primaries_basis(double lambda) {
    // Mallet & Yuksel, Dawson: Gaussian basis for RGB
    double r = exp(-0.5 * pow((lambda - 600.0) / 40.0, 2));
    double g = exp(-0.5 * pow((lambda - 550.0) / 30.0, 2));
    double b = exp(-0.5 * pow((lambda - 450.0) / 30.0, 2));
    return Vector3D(r, g, b);
}

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1.0;
  return reflectance / abs_cos_theta(*wi);
}

void MirrorBSDF::render_debugger_node() {
  if (ImGui::TreeNode(this, "Mirror BSDF")) {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO: proj3-2, part 3
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  
  return 1.0;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO: proj3-2, part 3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.

  double cosTheta = cos_theta(wi);
  
  return Vector3D();
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO: proj3-2, part 3
  // Implement microfacet model here.

  return Vector3D();
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO: proj3-2, part 3
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  *wi = cosineHemisphereSampler.get_sample(pdf);
  return MicrofacetBSDF::f(wo, *wi);
}

void MicrofacetBSDF::render_debugger_node() {
  if (ImGui::TreeNode(this, "Micofacet BSDF")) {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // Implement RefractionBSDF
  double eta = wo.z > 0 ? 1 / ior : ior;

  // Internal reflection
  if (1 - eta * eta * sin_theta2(wo) < 0) {
    reflect(wo, wi);
    *pdf = 1.0;
    return Vector3D();
  }

  *pdf = 1.0;
  refract(wo, wi, ior);
  return transmittance / abs_cos_theta(*wi) / (eta*eta);
}

void RefractionBSDF::render_debugger_node() {
  if (ImGui::TreeNode(this, "Refraction BSDF")) {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // Compute Fresnel coefficient and either reflect or refract based on it.
  if (!refract(wo, wi, ior)) {
    reflect(wo, wi);
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
  }

  double eta = wo.z > 0 ? 1 / ior : ior;
  double R0 = powf((ior - 1) / (ior + 1), 2);
  double R = R0 + (1 - R0) * powf(1 - abs_cos_theta(wo), 5);

  if (coin_flip(R)) {
    reflect(wo, wi);
    *pdf = R;
    return reflectance / abs_cos_theta(*wi);
  }

  *pdf = 1 - R;
  return (1 - R) * transmittance / abs_cos_theta(*wi) / (eta*eta);
}

void GlassBSDF::render_debugger_node() {
  if (ImGui::TreeNode(this, "Refraction BSDF")) {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {
  double eta = wo.z > 0 ? 1 / ior : ior;

  if (1 - eta * eta * sin_theta2(wo) < 0) {
    reflect(wo, wi);
    return false;
  }

  double sign = wo.z > 0 ? -1 : 1;
  *wi = Vector3D(-eta * wo.x, -eta * wo.y, sign * sqrt(1 - eta*eta*(1 - wo.z*wo.z)));
  return true;
}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior_o, double ior_i) {
  // Use Snell's Law to refract wo surface and store result ray in wi.
  double eta = ior_i / ior_o;

  if (1 - eta * eta * sin_theta2(wo) < 0) {
    reflect(wo, wi);
    return false;
  }

  double sign = wo.z > 0 ? -1 : 1;
  *wi = Vector3D(-eta * wo.x, -eta * wo.y, sign * sqrt(1 - eta*eta*(1 - wo.z*wo.z)));
  return true;
}

// Spectral Distribution Function
void SpectralBSDF::render_debugger_node() {
  if (ImGui::TreeNode(this, "Thin Film (Spectral) BSDF")) {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("Film IOR", &film_ior, 0.005);
    DragDouble("Thickness (nm)", &thickness, 1.0);
    DragDouble("Base IOR", &base_ior, 0.005);
    ImGui::TreePop();
  }
}

// Spectral primaries basis (Mallet & Yuksel, Dawson Eq. 3)
double SpectralBSDF::spd(double lambda) {
  Vector3D basis = spectral_primaries_basis(lambda);
  return reflectance.x * basis.x + reflectance.y * basis.y + reflectance.z * basis.z;
}

// Airy reflectance for thin film interference
double SpectralBSDF::airy_reflectance(const Vector3D& wo, double lambda, double thickness) const {
  double n1 = 1.0; // air
  double n2 = film_ior + 0.003/(lambda*lambda); // Cauchy dispersion, B=0.003
  double n3 = base_ior;
  double cos_theta1 = abs_cos_theta(wo);
  double sin_theta1 = sqrt(1 - cos_theta1*cos_theta1);
  double sin_theta2 = n1/n2 * sin_theta1;
  if (sin_theta2 > 1.0) sin_theta2 = 1.0; // TIR guard
  double cos_theta2 = sqrt(1 - sin_theta2*sin_theta2);
  double sin_theta3 = (n2/n3) * sin_theta2; // FIXED: correct Snell's law for n2->n3
  if (sin_theta3 > 1.0) sin_theta3 = 1.0;
  double cos_theta3 = sqrt(1 - sin_theta3*sin_theta3);
  
  // Fresnel coefficients (perpendicular) - using amplitude Fresnel equations
  double r12 = (n1*cos_theta1 - n2*cos_theta2) / (n1*cos_theta1 + n2*cos_theta2);
  double r23 = (n2*cos_theta2 - n3*cos_theta3) / (n2*cos_theta2 + n3*cos_theta3);
  
  double phi = 4.0 * M_PI * thickness * cos_theta2 / lambda; // no fudge factor, thickness in nm

  double r12_2 = r12*r12, r23_2 = r23*r23;
  double numerator   = r12_2 + r23_2 + 2.0*r12*r23*cos(phi);
  double denominator = 1.0 + r12_2*r23_2 - 2.0*r12*r23*cos(phi); // sign fixed
  double R = numerator / denominator;
  R = std::max(0.0, std::min(0.999, R)); // only clamp to avoid FP blow-ups
  
#ifdef REDUCE_THINFILM_REFLECTANCE
  // Apply scaling factor to lower overall reflectance
  // This can be adjusted between 0.2-0.6 for different visual effects
  const double REFLECTANCE_SCALE = 0.4;
  DEBUG_LOG("  Using reflectance reduction factor: " << REFLECTANCE_SCALE);
  R *= REFLECTANCE_SCALE;
  
  // Ensure reasonable bounds
  R = std::max(0.05, std::min(0.7, R));
#else
  // Without reduction, we still need to avoid extreme values
  // that could cause rendering artifacts
  R = std::max(0.01, std::min(0.99, R));
  DEBUG_LOG("  Using raw reflectance value: " << R);
#endif
  
  return R;
}

std::vector<double> SpectralBSDF::hero_sampler(double lambda_hero) {
  std::vector<double> samples;
  
  // Ensure hero wavelength is in our primary sampling range
  lambda_hero = std::max(400.0, std::min(700.0, lambda_hero));
  
  // Add the hero wavelength first
  samples.push_back(lambda_hero);
  
  // Add 9 more random samples distributed around the hero wavelength
  // within a 300nm range (roughly 400-700nm visible spectrum)
  const double range_width = 300.0;
  const int num_samples = 19;  // 9 more for a total of 10
  
  for (int i = 0; i < num_samples; i++) {
    // Create samples that vary by up to +/- 150nm from the hero
    // using random distribution to better capture spectral detail
    double offset = (random_uniform() * range_width) - (range_width/2.0);
    double sample = lambda_hero + offset;
    
    // Ensure samples stay in visible range
    sample = std::max(380.0, std::min(780.0, sample));
    samples.push_back(sample);
  }
  
  return samples;
}

Vector3D SpectralBSDF::f(const Vector3D wo, const Vector3D wi) {
  // Not used for delta BSDFs, return black
  return Vector3D();
}

Vector3D SpectralBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // Hero wavelength sampling - ensure good distribution across spectrum
  double lambda = 380.0 + random_uniform() * (780.0 - 380.0);
  std::vector<double> wavelengths = hero_sampler(lambda);
  
  // Debug: Print hero wavelength and samples - limit debug output
  static int debug_count = 0;
  bool should_debug = (debug_count < 5) && (random_uniform() < 0.001); // Only debug 0.1% of calls
  
  if (should_debug) {
    debug_count++;
    DEBUG_LOG("Hero λ: " << lambda << "nm, Samples: " 
              << wavelengths[0] << "nm " 
              << wavelengths[1] << "nm " 
              << wavelengths[2] << "nm " 
              << wavelengths[3] << "nm");
  }
  
  // Average Airy reflectance across hero samples
  double R = 0;
  for (double l : wavelengths) {
    if (should_debug) DEBUG_LOG("Computing reflectance for λ=" << l << "nm");
    double r_l = airy_reflectance(wo, l, thickness);
    R += r_l;
    
    // Debug: Print reflectance for each wavelength
    if (should_debug) DEBUG_LOG("  λ=" << l << "nm: R=" << r_l);
  }
  R /= wavelengths.size();
  
#ifdef REDUCE_THINFILM_REFLECTANCE  
  // Force reasonable reflectance values - lower the upper bound only if reduction is enabled
  R = std::min(0.6, std::max(0.1, R));
#endif
  
  // Debug: Print average reflectance
  if (should_debug) {
    DEBUG_LOG("  Avg R=" << R << ", film_ior=" << film_ior << ", thickness=" << thickness << "nm");
  }
  
  // Use Russian roulette for termination
  double russian_roulette_prob = 0.9; // 90% chance to continue, 10% to terminate
  
  // Reflection or transmission based on probability
  if (coin_flip(R)) {
    reflect(wo, wi);
    *pdf = R;
    if (should_debug) DEBUG_LOG("  REFLECT: pdf=" << *pdf);
    Vector3D F = reflectance * R; // scale the colour mask by Fresnel R
    return F / abs_cos_theta(*wi) / (*pdf); // correct energy balance
  } else if (base_bsdf && coin_flip(russian_roulette_prob)) {
    // Apply Russian roulette for recursive base material
    // Transmission to base material, use base BSDF
    double base_pdf;
    Vector3D f = base_bsdf->sample_f(wo, wi, &base_pdf);
    *pdf = (1 - R) * base_pdf * russian_roulette_prob;
    if (should_debug) DEBUG_LOG("  TRANSMIT TO BASE: pdf=" << *pdf);
    Vector3D F = transmittance * (1 - R); // scale by transmission probability
    return f * F / russian_roulette_prob / (*pdf); // correct energy balance
  } else {
    // Simple refraction through the film
    double n = 0;
    for (double l : wavelengths) {
      n += film_ior + 0.003/(l*l);
    }
    n /= wavelengths.size();
    refract(wo, wi, n);
    *pdf = 1 - R;
    double eta = wo.z > 0 ? 1/n : n;
    if (should_debug) DEBUG_LOG("  REFRACT: n=" << n << ", pdf=" << *pdf);
    Vector3D F = transmittance * (1 - R);
    return F / abs_cos_theta(*wi) / (eta*eta) / (*pdf); // correct energy balance
  }
}

Vector3D SpectralBSDF::sample_lambda() {
  // Simple caching to prevent excessive recalculation
  static Vector3D cached_result;
  static int call_count = 0;
  
  // Only recalculate every 100 calls to improve performance
  if (call_count++ % 100 == 0) {
    Vector3D f;
    double lambda = 380.0 + random_uniform() * (400.0 - 380.0); // Narrower range for better efficiency
    std::vector<double> sample = hero_sampler(lambda);
    for (double wavelength : sample) {
      f += 1.0 * ColorSpace::wave2xyz_table(wavelength);
    }
    f /= sample.size();
    cached_result = f;
  }
  
  return cached_result;
}

} // namespace CGL
