#ifndef CGAL_H
#define CGAL_H
// CGAL
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef CGAL::Aff_transformation_3<K> Transformation;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;


Transformation rotationMatrixZ(double angle);

Polyhedron rotation(Point& a, Point& b, Point& c, Point& d, Transformation& rotate);

#endif // CGAL_H
