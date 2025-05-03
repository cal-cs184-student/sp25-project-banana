#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 

  // Just copying the slides

  L_out = Vector3D(0);


  double pdf = 1/(2*PI); // p(w_j)
  for (int j = 0; j < num_samples; j++) {
      // f_r(p, w_j -> w_r)
	  Vector3D w_j = hemisphereSampler->get_sample();
	  Vector3D f_r = isect.bsdf->f(w_out, w_j);

      // Light from w_j direction
	  Intersection* light_j = new Intersection();
	  Ray* ray_j = new Ray(hit_p, w_j.unit());
      ray_j->min_t = EPS_F;
      if (!bvh->intersect(*ray_j, light_j)) continue;
      Vector3D L_i = light_j->bsdf->get_emission();

      L_out += f_r * L_i * dot(w_j, isect.n);
  }

  L_out /= pdf;
  L_out /= num_samples;

  return L_out;

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  // TODO (Part 3.4): Write your sampling loop here

  L_out = Vector3D(0);
  int N = 0;

  for (SceneLight* light : scene->lights) {
      int nsamples = light->is_delta_light() ? 1 : ns_area_light;
      N += nsamples;

      for (int i = 0; i < nsamples; i++) {
          Vector3D wi;
          double distToLight, pdf;
          Vector3D L_i = light->sample_L(hit_p, &wi, &distToLight, &pdf);
          wi.normalize();

          // Light behind hit surface
          if ((w2o * wi).z <= 0) {
              continue;
          }

          // Check if exists object closer than light source
          Intersection* isect_i = new Intersection();
          Ray r(hit_p, wi);
          r.min_t = EPS_F; // avoid self-intersection
          r.max_t = distToLight - EPS_F;
          if (bvh->intersect(r, isect_i)) {
              continue;
          }

		  // f = f_r(w_i -> w_r) * L_i(w_i) * cos(theta_i)
          // p = pdf (uniform)
          // ~ \sum f/p (found on slides)
          L_out += isect.bsdf->f(w_out, wi) * L_i * dot(wi, isect.n) / pdf;
          //L_out += isect.bsdf->f(w_out, wi) * L_i * dot(wi, isect.n) / pdf;
      }
  }

  return L_out / N;
}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r, const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light
  return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

    if (direct_hemisphere_sample) {
		return estimate_direct_lighting_hemisphere(r, isect);
    }
  return estimate_direct_lighting_importance(r, isect);
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out(0, 0, 0);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.


    // Always compute direct lighting first
    Vector3D direct = one_bounce_radiance(r, isect);
    
    if (isAccumBounces) {
        // ACCUMULATION MODE: Add direct + continue bouncing
        L_out += direct;
        
        // Russian Roulette only for indirect paths
        if (r.depth > 1 && coin_flip(0.7)) {
            Vector3D w_in;
            double pdf;
            Vector3D f = isect.bsdf->sample_f(w_out, &w_in, &pdf);
            
            // Convert to world space without explicit normalize()
            Vector3D wi_world = o2w * w_in;
            
            // Create ray using your existing pattern
            Ray new_ray;
            new_ray.o = hit_p + wi_world * EPS_F;
            new_ray.d = wi_world;
            new_ray.min_t = EPS_F;
            new_ray.max_t = INF_D;
            new_ray.depth = r.depth - 1;

            Intersection next_isect;
            if (bvh->intersect(new_ray, &next_isect)) {
                double cos_theta = fabs(dot(isect.n, wi_world));
                Vector3D indirect = at_least_one_bounce_radiance(new_ray, next_isect);
                L_out += (f * cos_theta * indirect) / (pdf * 0.7);
            }
        }
    } else {
        // VISUALIZATION MODE: Only count if at target depth
        if (r.depth == 1) {  // Assuming depth=1 means "first bounce"
            L_out += direct;
        }
    }


  return L_out;

}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.

  if (!bvh->intersect(r, &isect))
    return envLight ? envLight->sample_dir(r) : L_out;


  // TODO (Part 1): If the ray intersects a light, return the light's emission
  //L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.
  L_out = zero_bounce_radiance(r, isect);
  //L_out += one_bounce_radiance(r, isect);

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct
  L_out += at_least_one_bounce_radiance(r, isect);

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
    // TODO (Part 1.2):
    // Make a loop that generates num_samples camera rays and traces them
    // through the scene. Return the average Vector3D.
    // You should call est_radiance_global_illumination in this function.

    // Set true for part 5
    // Just keeping the previous solution
    // Although enable adaptive sampling when running the script
    bool adaptiveSampling = false;

    if (!adaptiveSampling) {

        Vector3D total_radiance = Vector3D(0);
		for (int i = 0; i < ns_aa; i++) {
			Vector2D sample = gridSampler->get_sample();
			double normalized_x = (x + sample.x) / sampleBuffer.w;
			double normalized_y = (y + sample.y) / sampleBuffer.h;
			Ray ray = camera->generate_ray(normalized_x, normalized_y);
      
			Vector3D sample_radiance = est_radiance_global_illumination(ray);
			total_radiance += sample_radiance;
		}
		Vector3D avg_radiance = total_radiance / ns_aa;
		sampleBuffer.update_pixel(avg_radiance, x, y);
		sampleCountBuffer[x + y * sampleBuffer.w] = ns_aa;
        return;
    }


  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"
  int num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel

  Vector3D total_radiance = Vector3D(0);

  int i = 0;
  double mu = 0, var = 0;
  while (i < num_samples) {
      double s1 = 0, s2 = 0;
      int n = samplesPerBatch;
      for (int k = 0; k < n; k++) {
          Vector2D sample = gridSampler->get_sample();
          double normalized_x = (x + sample.x) / sampleBuffer.w;
          double normalized_y = (y + sample.y) / sampleBuffer.h;
          Ray ray = camera->generate_ray(normalized_x, normalized_y);
          ray.depth = max_ray_depth;
          Vector3D radiance = est_radiance_global_illumination(ray);
          total_radiance += radiance;

          float x_k = radiance.illum();
          s1 += x_k;
          s2 += x_k * x_k;
      }
      i += n;

      mu = (mu + s1 / n) / 2;
      var = (var + (s2 - (s1 * s1) / n) / (n - 1)) / 2;

      // Check if converged (in specs)
      double I = 1.96 * sqrt(var / i);

      if (I <= maxTolerance * mu) break;
  }

  total_radiance /= i;
  sampleBuffer.update_pixel(total_radiance, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = i; // num samples = i

}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
