#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

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
  return Vector3D(0);
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO:
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305

	double R0 = powf((ior - 1) / (ior + 1), 2);
	double R = R0 + (1 - R0) * powf(1 - abs_cos_theta(wo), 5);

    if (coin_flip(R)) {
        *pdf = R;
		reflect(wo, wi);
		return reflectance;
    }

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

  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

    if (sin_theta2(wo) > ior*ior) {
        return false;
    }

	Vector3D N = Vector3D(0, 0, 1);
	*wi = ior * cross(N, cross(-N, wo))
       - N * sqrt(1 - ior * ior * sin_theta2(wo));

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
	for (int i = 0; i < spd.size(); i++) {
		if (lambda < spd[i]) {
			return spd[i];
		}
	}
	return spd[spd.size() - 1];
};

Vector3D SpectralBSDF::f(const Vector3D wo, const Vector3D wi) {
	return ;
}

Vector3D SpectralBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

	// Fresnel coefficient 
	// Fundamentals of Computer Graphics page 305

	double R0 = powf((ior - 1) / (ior + 1), 2);
	double R = R0 + (1 - R0) * powf(1 - abs_cos_theta(wo), 5);

    if (coin_flip(R)) {
        *pdf = R;
		reflect(wo, wi);
		return reflectance * sample_lambda();
    }

    *pdf = 1 - R;
	refract(wo, wi, ior);
  return transmittance * sample_lambda();
}

Vector3D SpectralBSDF::sample_lambda() {
    int N = 10;
    Vector3D f;
    for (int i = 0; i < N; i++) {
    	double lambda = random_uniform();
		f += uniform_spd(lambda) * to_xyz(lambda);
    }
    f /= N;
	return f;
}

Vector3D SpectralBSDF::to_xyz(double lambda) {
    return Vector3D();
}


} // namespace CGL
