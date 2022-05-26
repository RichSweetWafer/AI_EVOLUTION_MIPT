#ifndef RENDEREDOBJECT_H
#define RENDEREDOBJECT_H
#include "headers.h"

class RenderedObject
{
protected:
    QList<GLfloat> m_data;
    int m_count = 0;
    void add(const Point &v, const QVector3D &n);
    void polygon(Point x, Point y, Point z);
    void quad(Point v1, Point v2, Point v3);

public:
    RenderedObject();

    const GLfloat *constData() const { return m_data.constData(); }
    int count() const { return m_count; }
    int vertexCount() const { return m_count / 6; }
};

#endif // RENDEREDOBJECT_H
