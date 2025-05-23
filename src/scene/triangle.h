#ifndef CGL_STATICSCENE_TRIANGLE_H
#define CGL_STATICSCENE_TRIANGLE_H

#include "object.h"
#include "primitive.h"

namespace CGL { namespace SceneObjects {

/**
 * A single triangle from a mesh.
 * To save space, it holds a pointer back to the data in the original mesh
 * rather than holding the data itself. This means that its lifetime is tied
 * to that of the original mesh. The primitive may refer back to the mesh
 * object for other information such as normal, texcoord, material.
 */
class Triangle : public Primitive {
public:

  /**
   * Constructor.
   * Construct a mesh triangle with the given indicies into the triangle mesh.
   * \param mesh pointer to the mesh the triangle is in
   * \param v1 index of triangle vertex in the mesh's attribute arrays
   * \param v2 index of triangle vertex in the mesh's attribute arrays
   * \param v3 index of triangle vertex in the mesh's attribute arrays
   */
  Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3);

  Triangle() {};

  /**
   * Barycentric coordinates.
   * Compute the barycentric coordinates of point p with respect to the triangle
   * defined by vertices a, b, and c. The barycentric coordinates are returned
   * in u, v, and w.
   */
  void barycentric(const Vector3D& p, const Vector3D& a, const Vector3D& b, const Vector3D& c, float& u, float& v, float& w) const;

  /**
   * Get the world space bounding box of the triangle.
   * \return world space bounding box of the triangle
   */
  BBox get_bbox() const;

  /**
   * Ray - Triangle intersection.
   * Check if the given ray intersects with the triangle, no intersection
   * information is stored.
   * \param r ray to test intersection with
   * \return true if the given ray intersects with the triangle,
             false otherwise
   */
  bool has_intersection(const Ray& r) const;

  /**
   * Ray - Triangle intersection 2.
   * Check if the given ray intersects with the triangle, if so, the input
   * intersection data is updated to contain intersection information for the
   * point of intersection.
   * \param r ray to test intersection with
   * \param i address to store intersection info
   * \return true if the given ray intersects with the triangle,
             false otherwise
   */
  bool intersect(const Ray& r, Intersection* i) const;

  /**
   * Get BSDF.
   * In the case of a triangle, the surface material BSDF is stored in 
   * the mesh it belongs to. 
   */
  BSDF* get_bsdf() const { return bsdf; }

  /**
   * Draw with OpenGL (for visualizer)
   */
  void draw(const Color& c, float alpha) const;

  /**
   * Draw outline with OpenGL (for visualizer)
   */
  void drawOutline(const Color& c, float alpha) const;

  Vector3D p1, p2, p3;
  Vector3D n1, n2, n3;
  
  BSDF* bsdf;

  BBox bbox;
}; // class Triangle

} // namespace SceneObjects
} // namespace CGL

#endif //CGL_STATICSCENE_TRIANGLE_H
