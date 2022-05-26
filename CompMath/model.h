#ifndef MODEL_H
#define MODEL_H

#include "axis.h"

class Process : public RenderedObject
{
    // Figure and AABB tree
    Point a;
    Point b;
    Point c;
    Point d;
    Polyhedron polyhedron;
    Ray ray;

public:
    Process();
    double f(double rotate_angle, double ray_Yangle, double ray_Zanglem, bool render = true);
};

#endif // MODEL_H
