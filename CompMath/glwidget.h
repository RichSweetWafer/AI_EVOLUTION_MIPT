#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMatrix4x4>
#include "AL.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = nullptr);
    ~GLWidget();

    QSize minimumSizeHint() const override;
    QSize sizeHint() const override;
    vector step();
    QPair<double, double> predict(int rotateAngle,
                                  int rayYAngle,
                                  int rayZAngle);
public slots:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void cleanup();

signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);

protected:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int width, int height) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;

private:
    void setupVertexAttribs();
    void updateBuffer();
    // Mouse movement
    QPoint m_lastPos;
    int m_xRot = 0;
    int m_yRot = 0;
    int m_zRot = 0;

    // Objects & Active Learning

    MVN prior;
    BIRegressor* AL;
    Process* model;
    Axis axis;

    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_logoVbo;

    // Shader
    QOpenGLShaderProgram *m_program = nullptr;
    // Shader parameters' locations
    int m_projMatrixLoc = 0;
    int m_mvMatrixLoc = 0;
    int m_normalMatrixLoc = 0;
    int m_lightPosLoc = 0;
    int m_colorLoc = 0;

    QMatrix4x4 m_proj;
    QMatrix4x4 m_camera;
    QMatrix4x4 m_world;
};

#endif
