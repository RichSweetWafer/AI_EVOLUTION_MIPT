#include "renderedobject.h"

RenderedObject::RenderedObject(){}

void RenderedObject::add(const Point &v, const QVector3D &n)
{
    GLfloat *p = m_data.data() + m_count;
    *p++ = v.x() / 6;
    *p++ = v.y() / 6;
    *p++ = v.z() / 6;
    *p++ = n.x() / 6;
    *p++ = n.y() / 6;
    *p++ = n.z() / 6;
    m_count += 6;
}

void RenderedObject::polygon(Point v1, Point v2, Point v3)
{
    QVector3D n = QVector3D::normal(QVector3D(v3.x() - v1.x(), v3.y() - v1.y(), v3.z() - v1.z()),
                                    QVector3D(v2.x() - v1.x(), v2.y() - v1.y(), v2.z() - v1.z()));
    add(v1, n);
    add(v3, n);
    add(v2, n);
    n = QVector3D::normal(QVector3D(v2.x() - v1.x(), v2.y() - v1.y(), v2.z() - v1.z()),
                          QVector3D(v3.x() - v1.x(), v3.y() - v1.y(), v3.z() - v1.z()));
    add(v3, n);
    add(v1, n);
    add(v2, n);
}

void RenderedObject::quad(Point v1, Point v2, Point v3)
{
    polygon(v1, v2, v3);
    Point v(v3.x() - v2.x(), v3.y() - v2.y(), v3.z() - v2.z());
    polygon(v1, v, v3);
}
