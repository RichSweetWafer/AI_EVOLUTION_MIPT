// Tree of Polyhedron Triangle Facets for Distance Queries
#include <iostream>
#include <random>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
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
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
int main()
{
    std::random_device r;
    std::default_random_engine eng(r());
    std::uniform_real_distribution<double> disr(-10, 10);
    Point a(1.0, 0.0, 0.0);
    Point b(0.0, 1.0, 0.0);
    Point c(0.0, 0.0, 1.0);
    Point d(0.0, 0.0, 0.0);
    Point query(disr(eng),disr(eng),disr(eng));
    //Point query(1,1,1);
    Point l(0.2, 0.2, 0.2);
    Point h(disr(eng), disr(eng), disr(eng));
    Ray ray(h, l);
    std::cout << "Query point: " << query
            << "\nRay: " << ray << std::endl;
    Polyhedron polyhedron;
    polyhedron.make_tetrahedron(a, b, c, d);
    Tree tree(CGAL::faces(polyhedron).first, CGAL::faces(polyhedron).second, polyhedron);
    
    // Euclidean Distance
    double sqd = tree.squared_distance(query);
    std::cout << "Query.\nsquared distance: " << sqd << std::endl;
    
    // closest point
    // Point closest = tree.closest_point(query);
    // std::cout << "closest point: " << closest << std::endl;
    
    // closest point and primitive id
    Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
    Point closest_point = pp.first;
    Polyhedron::Face_handle f = pp.second; // closest primitive id
    std::cout << "closest point: " << closest_point << std::endl;
    std::cout << "closest triangle: ( "
              << f->halfedge()->vertex()->point() << " , "
              << f->halfedge()->next()->vertex()->point() << " , "
              << f->halfedge()->next()->next()->vertex()->point()
              << " )" << std::endl;

    // Ray intersection
    std::cout << "Ray.\n" << tree.number_of_intersected_primitives(ray)
            << " instersections with the ray" << std::endl;
    boost::optional<Tree::Intersection_and_primitive_id<Segment>::Type> intersection = tree.first_intersection(ray);
    if (intersection)
    {
        Point * ptr = boost::get<Point>(&(intersection->first));
        if (ptr)
        std::cout << "intersection point " << *ptr << std::endl;
    }  
    return 0;
}

// попробовать задать луч углами
// луч в нулевой коорде
// тетраэдр "вращается"
// Qt?