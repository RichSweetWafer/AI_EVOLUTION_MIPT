#include <QApplication>
#include <QSurfaceFormat>
#include <QScreen>
#include <QCommandLineParser>
#include <QCommandLineOption>

#include "glwidget.h"
#include "mainwindow.h"

class Model;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QCoreApplication::setApplicationName("Active Learning");
    /*
    QSurfaceFormat fmt;
    fmt.setDepthBufferSize(24);
    if (parser.isSet(multipleSampleOption))
        fmt.setSamples(4);
    if (parser.isSet(coreProfileOption)) {
        fmt.setVersion(3, 2);
        fmt.setProfile(QSurfaceFormat::CoreProfile);
    }
    QSurfaceFormat::setDefaultFormat(fmt);
    */
    MainWindow mainWindow;
    mainWindow.resize(mainWindow.sizeHint());
//    int desktopArea = QGuiApplication::primaryScreen()->size().width() *
//                      QGuiApplication::primaryScreen()->size().height();
//    int widgetArea = mainWindow.width() * mainWindow.height();
//    if (((float)widgetArea / (float)desktopArea) < 0.75f)
//        mainWindow.show();
//    else
//        mainWindow.showMaximized();
    mainWindow.showMaximized();
    return app.exec();
}



/* Кнопка Start - выполнить n шагов (хардкод) AL (визуализация каждого шага ~0.5 с)
 * Кнопка Step - выполнить шаг AL (загрузить и рендерить соответствующий стейт)
 * Кнопка Predict - выполнить предположение в указанной точке, сравнить результаты
 * Кнопка Reset - сбросить AL до начального состояния
 * Добавить логирование?
*/
