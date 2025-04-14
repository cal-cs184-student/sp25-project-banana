#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
    namespace SceneObjects {

        BVHAccel::BVHAccel(const std::vector<Primitive*>& _primitives,
            size_t max_leaf_size) {

            primitives = std::vector<Primitive*>(_primitives);
            root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
        }

        BVHAccel::~BVHAccel() {
            if (root)
                delete root;
            primitives.clear();
        }

        BBox BVHAccel::get_bbox() const { return root->bb; }

        void BVHAccel::draw(BVHNode* node, const Color& c, float alpha) const {
            if (node->isLeaf()) {
                for (auto p = node->start; p != node->end; p++) {
                    (*p)->draw(c, alpha);
                }
            }
            else {
                draw(node->l, c, alpha);
                draw(node->r, c, alpha);
            }
        }

        void BVHAccel::drawOutline(BVHNode* node, const Color& c, float alpha) const {
            if (node->isLeaf()) {
                for (auto p = node->start; p != node->end; p++) {
                    (*p)->drawOutline(c, alpha);
                }
            }
            else {
                drawOutline(node->l, c, alpha);
                drawOutline(node->r, c, alpha);
            }
        }

        BVHNode* BVHAccel::construct_bvh(std::vector<Primitive*>::iterator start,
            std::vector<Primitive*>::iterator end,
            size_t max_leaf_size) {

            // TODO (Part 2.1):
            // Construct a BVH from the given vector of primitives and maximum leaf
            // size configuration. The starter code build a BVH aggregate with a
            // single leaf node (which is also the root) that encloses all the
            // primitives.

            BBox bbox;
            for (auto p = start; p != end; p++) {
                bbox.expand((*p)->get_bbox());
            }

            BVHNode* node = new BVHNode(bbox);

            size_t size = end - start;
            if (size <= max_leaf_size) {
                node->start = start;
                node->end = end;
                return node;
            }

            // Find optimal split using median partitioning
            int axis = bbox.extent.x > bbox.extent.y ?
                (bbox.extent.x > bbox.extent.z ? 0 : 2) :
                (bbox.extent.y > bbox.extent.z ? 1 : 2);

            auto mid = start + size / 2;
            std::nth_element(start, mid, end,
                [axis](Primitive* a, Primitive* b) {
                    return a->get_bbox().centroid()[axis] < b->get_bbox().centroid()[axis];
                });

            node->l = construct_bvh(start, mid, max_leaf_size);
            node->r = construct_bvh(mid, end, max_leaf_size);

            return node;
        }

        bool BVHAccel::has_intersection(const Ray& ray, BVHNode* node) const {
            // TODO (Part 2.3):
            // Fill in the intersect function.
            // Take note that this function has a short-circuit that the
            // Intersection version cannot, since it returns as soon as it finds
            // a hit, it doesn't actually have to find the closest hit.

            double t_min = ray.min_t, t_max = ray.max_t;
            if (!node->bb.intersect(ray, t_min, t_max)) return false;

            if (node->isLeaf()) {
                for (auto p = node->start; p != node->end; ++p) {
                    total_isects++;
                    if ((*p)->has_intersection(ray)) return true;
                }
                return false;
            }

            // Traverse closer child first
            BVHNode* first = node->l, * second = node->r;
            double t1, t2;
            bool hit1 = first->bb.intersect(ray, t_min, t1);
            bool hit2 = second->bb.intersect(ray, t_min, t2);

            if (hit1 && (t1 < t2 || !hit2)) {
                if (has_intersection(ray, first)) return true;
                return hit2 && has_intersection(ray, second);
            }
            else if (hit2) {
                if (has_intersection(ray, second)) return true;
                return hit1 && has_intersection(ray, first);
            }
            return false;
            return has_intersection(ray, node->l) || has_intersection(ray, node->r);
        }

        bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode* node) const {
            // TODO (Part 2.3):
            // Fill in the intersect function.

            std::stack<BVHNode*> stack;
            stack.push(node);
            bool hit = false;

            while (!stack.empty()) {
                BVHNode* curr = stack.top();
                stack.pop();

                double t_min = ray.min_t, t_max = ray.max_t;
                if (!curr->bb.intersect(ray, t_min, t_max)) continue;

                if (curr->isLeaf()) {
                    for (auto p = curr->start; p != curr->end; ++p) {
                        total_isects++;
                        hit = (*p)->intersect(ray, i) || hit;
                    }
                }
                else {
                    // Push children in reverse order for LIFO stack
                    stack.push(curr->r);
                    stack.push(curr->l);
                }
            }
            return hit;

        } // namespace SceneObjects
    } // namespace CGL
}