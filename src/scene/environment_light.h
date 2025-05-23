#ifndef CGL_STATICSCENE_ENVIRONMENTLIGHT_H
#define CGL_STATICSCENE_ENVIRONMENTLIGHT_H

#include "pathtracer/sampler.h"
#include "util/image.h"
#include "scene.h"

namespace CGL { namespace SceneObjects {

// An environment light can be thought of as an infinitely big sphere centered
// around your scene, radiating light on the scene in a pattern matching some
// image. This is commonly used for low-cost renderings of complex backrounds or
// environments (fancy church, overgrown forest, etc) that would be difficult to
// model in the scene.
class EnvironmentLight : public SceneLight {
public:
  EnvironmentLight(const HDRImageBuffer* envMap);
  ~EnvironmentLight();
  /**
    * In addition to the work done by sample_dir, this function also has to
    * choose a ray direction to sample along. You should initially use a uniform
    * randomly distributed direction, but later you'll be required to implement
    * this with importance sampling. The requirements are:
    * - Sample within a pixel with probability proportional to the total power
    *   going through that pixel (so proportional to the pixel's area when
    *   projected onto the sphere, multiplied by the sum of the components of its
    *   spectrum). This is a discrete probability distribution, so there's no
    *   formula.
    * - Don't take linear time to generate a single sample! You'll be calling
    *   this a LOT; it should be fast.
    */
  Vector3D sample_L(const Vector3D p, Vector3D* wi, double* distToLight,
    double* pdf) const;
  bool is_delta_light() const { return false; }
  /**
    * Returns the color found on the environment map by travelling in a specific
    * direction. This entails:
    * - Converting the 3-d ray direction into a 2-d image coordinate.
    * - Performing bilinear interpolation to determine the color at the given
    *   coordinate.
    * - Handling the edge cases correctly (what happens if you wrap around the
    *   environment map horizontally? What about vertically?).
    */
  Vector3D sample_dir(const Ray& r) const;

  // New methods for controlling environment lighting
  void set_constant_color(bool use_constant, const Vector3D& color = Vector3D(1.0, 1.0, 1.0));
  void set_brightness(double scale);

private:
  const HDRImageBuffer* envMap;
  UniformGridSampler2D sampler_uniform2d;
  UniformSphereSampler3D sampler_uniform_sphere;

  void init();
  double* pdf_envmap, * marginal_y, * conds_y;
  double brightness_scale;    // scaling factor for controlling HDR brightness
  
  // New fields for constant environment mode
  bool use_constant_color;
  Vector3D constant_color;

  Vector2D dir_to_theta_phi(const Vector3D dir) const;
  Vector3D theta_phi_to_dir(const Vector2D& theta_phi) const;

  Vector2D theta_phi_to_xy(const Vector2D& theta_phi) const;
  Vector2D xy_to_theta_phi(const Vector2D& xy) const;

  Vector3D bilerp(const Vector2D& xy) const;
  void save_probability_debug();

}; // class EnvironmentLight

} // namespace SceneObjects
} // namespace CGL

#endif //CGL_STATICSCENE_ENVIRONMENTLIGHT_H
