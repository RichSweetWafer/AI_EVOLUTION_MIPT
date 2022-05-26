#include "window.h"

Window::Window(MainWindow *mw)
    : mainWindow(mw)
{
    glWidget = new GLWidget;

    QVBoxLayout* mainLayout = new QVBoxLayout;

    QHBoxLayout* buttonLayout = new QHBoxLayout;

    startBtn = new QPushButton(tr("Start"), this);
    connect(startBtn, &QPushButton::clicked, this, &Window::onStart);
    buttonLayout->addWidget(startBtn);

    stepBtn = new QPushButton(tr("Step"), this);
    connect(stepBtn, &QPushButton::clicked, this, &Window::onStep);
    buttonLayout->addWidget(stepBtn);

    predictBtn = new QPushButton(tr("Predict"), this);
    connect(predictBtn, &QPushButton::clicked, this, &Window::onPredict);
    buttonLayout->addWidget(predictBtn);

    dataLabel = new QLabel(tr(""), this);
    buttonLayout->addWidget(dataLabel);

    QWidget *btns = new QWidget;
    btns->setLayout(buttonLayout);
    btns->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);


    QVBoxLayout* container = new QVBoxLayout;
    container->addWidget(btns);
    container->addWidget(glWidget);

    QWidget *w = new QWidget;
    w->setLayout(container);
    mainLayout->addWidget(w);

    setLayout(mainLayout);
}

void Window::onStart()
{
    QDateTime dateTime1 = QDateTime::currentDateTime();
    for (int i = 0; i < 3000; i++)
    {
        onStep();
//        QPair<double, double> pair = glWidget->predict(19,
//                                                       -1,
//                                                       4);
        //if (i % 100 == 0)
            qDebug("i: %d", i);
//        qDebug("Predict: %d, %d, %d, act: %f  pred %f loss: %f", 19, -1, 4, pair.first, pair.second, 100.0 - pair.second / pair.first * 100);
    }
    QDateTime dateTime2 = QDateTime::currentDateTime();
    qint64 millisecondsDiff = dateTime1.msecsTo(dateTime2);
    qDebug("diff: %lld", millisecondsDiff);
}

void Window::onStep()
{
    vector ans = glWidget->step();
    stateText(ans[0], ans[1], ans[2], ans[3]);
}

void Window::onPredict()
{
    bool ok;
    int rotateAngle = QInputDialog::getInt(this,
                                           tr("Object rotatrion angle"),
                                           tr("Value (from 0 to 359):"),
                                           0, 0, 359, 1, &ok);
    if (ok)
    {
        int rayYAngle = QInputDialog::getInt(this,
                                             tr("Ray angle axis Y"),
                                             tr("Value (from -10 to 10):"),
                                             0, -10, 10, 1, &ok);
        if (ok)
        {
            int rayZAngle = QInputDialog::getInt(this,
                                                 tr("Ray angle axis Z"),
                                                 tr("Value (from -10 to 10):"),
                                                 0, -10, 10, 1, &ok);
            if (ok)
            {
                QPair<double, double> pair = glWidget->predict(rotateAngle,
                                                               rayYAngle,
                                                               rayZAngle);
                stateText(rotateAngle, rayYAngle, rayZAngle, pair);
            }
        }
    }
}

void Window::stateText(int rotateAngle, int rayYAngle,
                          int rayZAngle, QPair<double, double> pair)
{
    QString str = "Rotatrion angle: ";
    str.append(QString::number(rotateAngle));
    str.append("  Ray angle Y: ");
    str.append(QString::number(rayYAngle));
    str.append("  Ray angle Z: ");
    str.append(QString::number(rayZAngle));
    str.append("  Predicted distance: ");
    str.append(QString::number(pair.second));
    str.append("  Actual distance: ");
    str.append(QString::number(pair.first));
    str.append("\n");
    dataLabel->setText(str);
}

void Window::stateText(int rotateAngle, int rayYAngle,
                  int rayZAngle, double distance)
{
    QString str = "Rotatrion angle: ";
    str.append(QString::number(rotateAngle));
    str.append("  Ray angle Y: ");
    str.append(QString::number(rayYAngle));
    str.append("  Ray angle Z: ");
    str.append(QString::number(rayZAngle));
    str.append("  Distance: ");
    str.append(QString::number(distance));
    str.append("\n");
    dataLabel->setText(str);
}

