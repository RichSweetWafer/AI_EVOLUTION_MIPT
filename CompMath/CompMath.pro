HEADERS       = glwidget.h \
                AL.h \
                axis.h \
                boost.h \
                cgal.h \
                headers.h \
                model.h \
                renderedobject.h \
                window.h \
                mainwindow.h
SOURCES       = glwidget.cpp \
                AL.cpp \
                axis.cpp \
                boost.cpp \
                cgal.cpp \
                main.cpp \
                model.cpp \
                renderedobject.cpp \
                window.cpp \
                mainwindow.cpp

QT += widgets opengl openglwidgets

# install
target.path = $$[QT_INSTALL_EXAMPLES]/opengl/hellogl2
INSTALLS += target

INCLUDEPATH += /usr/include
