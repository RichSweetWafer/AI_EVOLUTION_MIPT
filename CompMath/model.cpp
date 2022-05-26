#include "model.h"

Process::Process():
    a(Point(3.0, -0.25, -0.5)),
    b(Point(2.0, 0.75, -0.5)),
    c(Point(2.0, -0.25, 0.5)),
    d(Point(2.0, -0.25, -0.5))
{
    m_data.resize(2500 * 6);
    polyhedron.make_tetrahedron(a, b, c, d);

    // Object
    polygon(a, b, c);
    polygon(a, b, d);
    polygon(a, c, d);
    polygon(b, c, d);
};

double Process::f(double rotate_angle, double ray_Yangle, double ray_Zangle, bool render)
{
    // Temp points
    Point A = a;
    Point B = b;
    Point C = c;
    Point D = d;
    Point result;

    // Transform matrix
    double rotate_rad = rotate_angle * boost::math::constants::pi<double>() / 180;
    Transformation rotate{ rotationMatrixZ(rotate_rad) };

    // Ray
    double rayY_rad = ray_Yangle * boost::math::constants::pi<double>() / 180;
    double rayZ_rad = ray_Zangle * boost::math::constants::pi<double>() / 180;
    ray = Ray(Point(0.0, 0.0, 0.0),
        Point(0.5, tan(rayY_rad) * 0.5, tan(rayZ_rad) * 0.5));

    Polyhedron current = rotation(A, B, C, D, rotate);
    // For OpenGL
    m_count = 0;
    if (render)
    {
        polygon(A, B, C);
        polygon(A, B, D);
        polygon(A, C, D);
        polygon(B, C, D);
    }
//    qDebug("A: %f %f %f", A.x(), A.y(), A.z());
//    qDebug("B: %f %f %f", B.x(), B.y(), B.z());
//    qDebug("C: %f %f %f", C.x(), C.y(), C.z());
//    qDebug("D: %f %f %f", D.x(), D.y(), D.z());

    Tree tree(CGAL::faces(current).first, CGAL::faces(current).second, current);

    boost::optional<Tree::Intersection_and_primitive_id<Segment>::Type> intersection = tree.first_intersection(ray);
    if (intersection)
    {
        Point * ptr = boost::get<Point>(&(intersection->first));
        if (ptr) result = *ptr;
        else result = Point(0, 0, 0);//ray.point(10000.0);
    }
    else result = Point(0, 0, 0);//ray.point(10000.0);
    Point v(result.x(), result.y(), result.z() + 0.01);
    if (render)
    {
        polygon(Point(0.0, 0.0, 0.0), ray.point(10000.0), v);
        polygon(Point(0.0, 0.0, 0.0), v, Point(0.0, 0.0, 0.01));
    }
//    qDebug("\nIntersection = (%f, %f, %f)\nRotation angle = %f\n(Rad = %f)\nRay angle (y-axis) = %f  (Rad = %f)\n",
//           result.x(), result.y(), result.z(), rotate_angle, rotate_rad, ray_angle, ray_rad);
    double y = CGAL::squared_distance(ray.source(), result);
    return y;
}



