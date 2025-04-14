#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bounding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

	if (t0 > t1) std::swap(t0, t1); // make sure t0 is less than t1

	double tmin = (min.x - r.o.x) / r.d.x;
	double tmax = (max.x - r.o.x) / r.d.x;

	if (tmin > tmax) std::swap(tmin, tmax); // might not be necessary but for safety

	double tmin_y = (min.y - r.o.y) / r.d.y;
	double tmax_y = (max.y - r.o.y) / r.d.y;
	if (tmin_y > tmax_y) std::swap(tmin_y, tmax_y);
	tmin = std::max(tmin, tmin_y);
	tmax = std::min(tmax, tmax_y);

	double tmin_z = (min.z - r.o.z) / r.d.z;
	double tmax_z = (max.z - r.o.z) / r.d.z;
	if (tmin_z > tmax_z) std::swap(tmin_z, tmax_z);
	tmin = std::max(tmin, tmin_z);
	tmax = std::min(tmax, tmax_z);

	if (tmin > tmax || tmax < t0 || tmin > t1) return false;

	// boundary of the box only
	// don't change r.max_t, there is no guarantee of intersection with an object
	t0 = std::max(tmin, t0);
	t1 = std::min(tmax, t1);

  return true;
}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
