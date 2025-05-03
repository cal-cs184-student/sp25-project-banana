#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"
#include "color_vision.cpp"

using std::max;
using std::min;
using std::swap;
using std::cout;
using std::endl;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return reflectance / abs_cos_theta(wi);
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO:
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1.0;
  return f(wo, *wi);
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
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
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.



  *wi = cosineHemisphereSampler.get_sample(pdf);

  return MicrofacetBSDF::f(wo, *wi);
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return transmittance;
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO:
  // Implement RefractionBSDF

    *pdf = 1.0;
	refract(wo, wi, ior);
  return f(wo, *wi);
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
	double R0 = powf((ior - 1) / (ior + 1), 2);
	double R = R0 + (1 - R0) * powf(1 - abs_cos_theta(wo), 5);

    if (coin_flip(R)) {
		return reflectance;
    }
    return transmittance;
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO:
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305

	double R0 = powf((ior - 1) / (ior + 1), 2);
	double R = R0 + (1 - R0) * powf(1 - abs_cos_theta(wo), 5);

    if (coin_flip(R)) {
        std::cout << "r";
        *pdf = R;
		reflect(wo, wi);
		return reflectance;
    }

    std::cout << "t";
    *pdf = 1 - R;
	refract(wo, wi, ior);
  return transmittance;
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {
    return refract(wo, wi, ior, 1);
}
bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior_o, double ior_i) {

  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

    if (ior_o*ior_o*sin_theta2(wo) > ior_i*ior_i) {
        return false;
    }

	Vector3D N = Vector3D(0, 0, 1);
	*wi = ior_o/ior_i * cross(N, cross(-N, wo))
       - N * sqrt(1 - (ior_o*ior_o/(ior_i*ior_i)) * sin_theta2(wo));

  return true;

}

// Spectral Distribution Function
void SpectralBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

double SpectralBSDF::black_body_spd(double lambda) {
	  lambda = lambda * 1e-9; // convert nm to m
	  double K = 2 * PLANCK_CONSTANT * SPEED_OF_LIGHT * SPEED_OF_LIGHT / pow(lambda, 5);
      double T = 500;
	  double exp_term = exp(PLANCK_CONSTANT * SPEED_OF_LIGHT / (lambda * BOLTZMANN_CONSTANT * T));
	  return K / (exp_term - 1);
}

 double SpectralBSDF::custom_spd(double lambda) { 
   // assume spd is ordered
 	for (int i = 0; i < spd.size() - 5; i++) {
 		if (lambda < spd[i]) {
 			return spd[i];
 		}
 	}
 	return spd[spd.size() - 5];
 };

Vector3D SpectralBSDF::f(const Vector3D wo, const Vector3D wi) {
	std::cout << "f called" << std::endl;

	double R0 = powf((ior - 1) / (ior + 1), 2);
	double R = R0 + (1 - R0) * powf(1 - abs_cos_theta(wo), 5);
    Vector3D spectral_response = sample_lambda();
	if (coin_flip(R)) {
    // takes the reflectance relative to the wavelength, converts to RGB before entering the pipeline
		return reflectance * spectral_response;
	}
	return transmittance * spectral_response;
}

Vector3D SpectralBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

	// Fresnel coefficient 
	// Fundamentals of Computer Graphics page 305
	double R0 = powf((ior - 1) / (ior + 1), 2);
	double R = R0 + (1 - R0) * powf(1 - abs_cos_theta(wo), 5);
    Vector3D spectral_response = sample_lambda(); 

    if (coin_flip(R)) {
        *pdf = R;
		reflect(wo, wi);
        //std::cout << "reflected" << std::endl;
        // reflectance is relative to the wavelength and we convert back to RGB before entering the graphics pipeline
        return reflectance;
    }

    *pdf = 1 - R;
	refract(wo, wi, ior);
    //std::cout << "transmitted" << std::endl;

    // same as reflectance
    return transmittance;
}

std::vector<double> SpectralBSDF::hero_sampler(double lambda) {

  // not a very good hero sampler
  // just gets the closest 5 wavelengths
  // to the given wavelength
  double lambda2 = lambda;
  if (lambda - 20 < 380) {
    lambda2 += 20;
  } else if (lambda + 20 > 830) {
    lambda2 -= 20;
  } 
  std::vector<double> result;
  result = {lambda2 - 20, lambda2 - 10, lambda, lambda2 + 10, lambda2 + 20};
  return result;
}

Vector3D SpectralBSDF::sample_lambda() {
    // hero sample! >> straightforward imo
    Vector3D f;
    // gets a random lambda between 380 and 830
    double lambda = 380.0 + random_uniform() * (830.0 - 380.0);
    std::vector<double> sample = hero_sampler(lambda);
    for (double wavelength : sample) {
      // this should not be 1.0, this is a placeholder value because i noticed that f was getting way to small with the uniform sampler
      // not sure how to improve
      f += 1.0 * ColorSpace::wave2xyz_table(wavelength);
    }
    f /= sample.size();
	return f;
}




} // namespace CGL
