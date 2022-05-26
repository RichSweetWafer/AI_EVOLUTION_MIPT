#include "cgal.h"

Transformation rotationMatrixZ(double angle)
{
    const double cosa = cos(angle);
    const double sina = sin(angle);
    return Transformation(
            cosa, -sina, 0.0,
            sina,  cosa, 0.0,
            0.0,   0.0,  1.0);
}

Polyhedron rotation(Point& a, Point& b, Point& c, Point& d, Transformation& rotate)
{
    Polyhedron polyhedron;
    double centroid_x = double((a.x() + b.x() + c.x() + d.x()) / 4);
    double centroid_y = double((a.y() + b.y() + c.y() + d.y()) / 4);
    double centroid_z = double((a.z() + b.z() + c.z() + d.z()) / 4);
    Point A(a.x() - centroid_x, a.y() - centroid_y, a.z() - centroid_z);
    Point B(b.x() - centroid_x, b.y() - centroid_y, b.z() - centroid_z);
    Point C(c.x() - centroid_x, c.y() - centroid_y, c.z() - centroid_z);
    Point D(d.x() - centroid_x, d.y() - centroid_y, d.z() - centroid_z);
    A = rotate(A);
    B = rotate(B);
    C = rotate(C);
    D = rotate(D);
    a = Point(A.x() + centroid_x, A.y() + centroid_y, A.z() + centroid_z);
    b = Point(B.x() + centroid_x, B.y() + centroid_y, B.z() + centroid_z);
    c = Point(C.x() + centroid_x, C.y() + centroid_y, C.z() + centroid_z);
    d = Point(D.x() + centroid_x, D.y() + centroid_y, D.z() + centroid_z);
    polyhedron.make_tetrahedron(a, b, c, d);
    return polyhedron;
}
