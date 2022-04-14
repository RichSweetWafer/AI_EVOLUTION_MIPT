// Tree of Polyhedron Triangle Facets for Distance Queries
#include <iostream>
#include <random>
#include <cmath>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Aff_transformation_3.h>
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

Transformation rotationMatrixZ(double angle)
{
    const double cosa = cos(angle);
    const double sina = sin(angle);
    return Transformation(
            cosa, -sina, 0.0,
            sina,  cosa, 0.0,
            0.0,   0.0,  1.0);
}

void Rotation(Point& a, Point& b, Point& c, Point& d, Transformation& rotate)
{
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

}

int main()
{
    std::random_device r;
    std::default_random_engine eng(r());
    std::uniform_real_distribution<double> disr(-0.08, 0.08); // approx 5 degrees
    Transformation rotate{rotationMatrixZ(0.5)};
    Point a(3.0, -0.5, -0.5);
    Point b(2.0, 0.5, -0.5);
    Point c(2.0, -0.5, 0.5);
    Point d(2.0, -0.5, -0.5);
    Point query(disr(eng),disr(eng),disr(eng));
    //Point query(1,1,1);
    Point l(0.0, 0.0, 0.0);
    Point h(0.5, tan(disr(eng)) * 0.5, 0.0);
    Ray ray(l, h);
    std::cout << "Query point: " << query
            << "\nRay: " << ray.to_vector() << std::endl;

    // std::cout << a << ':' << b << ':' << c<< ':' << d << std::endl;
    // std::cout << "rotate\n";
    // Rotation(a, b, c, d, rotate);
    // std::cout << a << ':' << b << ':' << c<< ':' << d << std::endl;
            
    Polyhedron polyhedron;
    polyhedron.make_tetrahedron(a, b, c, d);
    Tree tree(CGAL::faces(polyhedron).first, CGAL::faces(polyhedron).second, polyhedron);
    //polyhedron.create_center_vertex();
    
    // Euclidean Distance
    double sqd = tree.squared_distance(query);
    std::cout << "\nQuery.\nsquared distance: " << sqd << '\n' << std::endl;
    
    // // closest point
    //  Point closest = tree.closest_point(query);
    //  std::cout << "closest point: " << closest << std::endl;
    
    // // closest point and primitive id
    // Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
    // Point closest_point = pp.first;
    // Polyhedron::Face_handle f = pp.second; // closest primitive id
    // std::cout << "closest point: " << closest_point << std::endl;
    // std::cout << "closest triangle: ( "
    //           << f->halfedge()->vertex()->point() << " , "
    //           << f->halfedge()->next()->vertex()->point() << " , "
    //           << f->halfedge()->next()->next()->vertex()->point()
    //           << " )" << std::endl;

    // Ray intersection
    for (size_t i = 0; i < 5; i++)
    {
        std::cout << "\nA = (" << a.x() << ',' << a.y() << ',' << a.z() << ")\n"
                     << "B = (" << b.x() << ',' << b.y() << ',' << b.z() << ")\n"
                     << "C = (" << c.x() << ',' << c.y() << ',' << c.z() << ")\n"
                     << "D = (" << d.x() << ',' << d.y() << ',' << d.z() << ')'
                     << std::endl;
        std::cout << "Ray.\n" << tree.number_of_intersected_primitives(ray)
            << " instersections with the ray" << std::endl;
        boost::optional<Tree::Intersection_and_primitive_id<Segment>::Type> intersection = tree.first_intersection(ray);
        if (intersection)
        {
            Point * ptr = boost::get<Point>(&(intersection->first));
            if (ptr)
                std::cout << "The first intersection point " << *ptr << std::endl;
        }
        Rotation(a, b, c, d, rotate);
        polyhedron.clear();
        polyhedron.make_tetrahedron(a, b, c, d);
        tree.rebuild(CGAL::faces(polyhedron).first, CGAL::faces(polyhedron).second, polyhedron);
    }
    
    
    return 0;
}

// попробовать задать луч углами // DONE
// луч в нулевой коорде // DONE 
// тетраэдр "вращается" // harder than i thought
// Qt?