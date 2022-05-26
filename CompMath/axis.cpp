#include "axis.h"

Axis::Axis()
{
    m_data.resize(576);
    /*polygon = 36;
     * quad = 72;
     *resize = 360 + 216
    */

    // Axis
    Point z(0.0, 0.0, 0.0);
    Point Y1(0.0, 6.0, 0.0);
    Point Y2(0.01, 6.0, 0.0);
    Point X1(6.0, 0.0, 0.0);
    Point X2(6.0, 0.01, 0.0);
    Point Z1(0.0, 0.0, 6.0);
    Point Z2(0.0, 0.01, 6.0);
    quad(z, Y1, Y2);
    quad(z, X1, X2);
    quad(z, Z1, Z2);
    // Letter X
    Point x1(0.1, 0.0, 0.1);
    Point x2(0.2, 0.0, 0.1);
    Point x3(0.3, 0.0, 0.1);
    Point x4(0.4, 0.0, 0.1);
    Point x5(0.1, 0.0, 0.4);
    Point x6(0.2, 0.0, 0.4);
    Point x7(0.3, 0.0, 0.4);
    Point x8(0.4, 0.0, 0.4);
    polygon(x1, x2, x8);
    polygon(x1, x7, x8);
    polygon(x3, x4, x5);
    polygon(x4, x5, x6);
    // Letter Y
    Point y1(0.0, 0.2, 0.1);
    Point y2(0.0, 0.3, 0.1);
    Point y3(0.0, 0.2, 0.3);
    Point y4(0.0, 0.3, 0.3);
    polygon(y1, y2, y3);
    polygon(y3, y4, y2);
    Point y5(0.0, 0.1, 0.4);
    Point y6(0.0, 0.15, 0.4);
    Point y7(0.0, 0.25, 0.3);
    polygon(y5, y6, y7);
    polygon(y3, y7, y5);
    Point y8(0.0, 0.4, 0.4);
    Point y9(0.0, 0.35, 0.4);
    Point y10(0.0, 0.25, 0.3);
    polygon(y8, y9, y10);
    polygon(y4, y10, y8);

}
