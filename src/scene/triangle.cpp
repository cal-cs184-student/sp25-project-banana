#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

void Triangle::barycentric(const Vector3D& p, const Vector3D& a, const Vector3D& b, const Vector3D& c, float& u, float& v, float& w) const
{
	Vector3D v0 = b - a, v1 = c - a, v2 = p - a;
	float d00 = dot(v0, v0);
	float d01 = dot(v0, v1);
	float d11 = dot(v1, v1);
	float d20 = dot(v2, v0);
	float d21 = dot(v2, v1);
	float denom = d00 * d11 - d01 * d01;
	v = (d11 * d20 - d01 * d21) / denom;
	w = (d00 * d21 - d01 * d20) / denom;
	u = 1.0f - v - w;
}

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.

	Vector3D e1 = p2 - p1, e2 = p3 - p1, s = r.o - p1;
	Vector3D s1 = cross(r.d, e2), s2 = cross(s, e1);

	double denom = dot(s1, e1);
	if (denom == 0)
		return false;
	double t = dot(s2, e2) / denom;
	if (t < r.min_t || t > r.max_t)
		return false;

	float alpha = dot(s1, s) / denom;
	float beta = dot(s2, r.d) / denom;
	if (alpha < 0 || beta < 0 || alpha + beta > 1)
		return false;

	r.max_t = t;
	return true;
}

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly 

	// t = dot(p - r.o, n)/dot(r.d, n)
	//	 = dot(r.o - p1, -n) / dot(r.d, n)
	//	 = dot(r.o - p1, cross(e1, e2)) / dot(r.d, cross(e2, e1))
	//	 = dot(cross(r.o - p1, e1), e2) / dot(cross(r.d, e2), e1)

	// alpha = dot(cross(r.o - p1, e1), r.d) / dot(cross(r.d, e2), e1)
	// beta = dot(cross(r.o - p1, e2), r.d) / dot(cross(r.d, e1), e2)
	Vector3D e1 = p2 - p1, e2 = p3 - p1, s = r.o - p1;
	Vector3D s1 = cross(r.d, e2), s2 = cross(s, e1);

	double denom = dot(s1, e1);
	if (denom == 0)
		return false;
	double t = dot(s2, e2) / denom;
	if (t < r.min_t || t > r.max_t)
		return false;

	float alpha = dot(s1, s) / denom;
	float beta = dot(s2, r.d) / denom;
	if (alpha < 0 || beta < 0 || alpha + beta > 1)
		return false;

	r.max_t = t;
	isect->bsdf = get_bsdf();
	isect->t = t;
	isect->primitive = this;
	isect->n = n1 * alpha + n2 * beta + n3 * (1 - alpha - beta);
	return true;
}

void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
