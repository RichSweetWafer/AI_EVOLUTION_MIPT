#include "mainwindow.h"
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>

MainWindow::MainWindow()
{
    QMenuBar *menuBar = new QMenuBar;
    QMenu *menuWindow = menuBar->addMenu(tr("&Menu"));
    QAction *resetAct = new QAction(menuWindow);
    resetAct->setText(tr("Reset"));
    resetAct->setShortcut(tr("Ctrl+R"));
    menuWindow->addAction(resetAct);
    connect(resetAct, &QAction::triggered, this, &MainWindow::onReset);

    setMenuBar(menuBar);
    onReset();
}

void MainWindow::onReset()
{
    if (centralWidget())
    {
        QWidget* ptr = centralWidget();
        delete ptr;
    }
    setCentralWidget(new Window(this));
}
