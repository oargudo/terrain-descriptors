QT += core gui opengl widgets
greaterThan(QT_MAJOR_VERSION, 5): QT += openglwidgets

CONFIG += c++11

INCLUDEPATH += akirmse-mountains

SOURCES += \
    core.cpp \
    heightfield-breach.cpp \
    heightfield-fill.cpp \
    heightfield-metrics-hydrology.cpp \
    heightfield-metrics-landforms.cpp \
    heightfield-metrics-local.cpp \
    heightfield-metrics-orometry.cpp \
    heightfield-metrics-visibility.cpp \
    heightfield.cpp \
    integerfield.cpp \
    isolation.cpp \
    main.cpp \
    mainwindow-draw-metrics.cpp \
    mainwindow.cpp \
    ppa.cpp \
    riversnet.cpp \
    scalarfield.cpp \
    terrainanalysis.cpp \
    terrainwidget.cpp \
    akirmse-mountains/divide_tree.cpp \
    akirmse-mountains/domain_map.cpp \
    akirmse-mountains/island_tree.cpp \
    akirmse-mountains/latlng.cpp \
    akirmse-mountains/line_tree.cpp \
    akirmse-mountains/tile.cpp \
    akirmse-mountains/tree_builder.cpp \
    akirmse-mountains/utm_coordinate_system.cpp

HEADERS += \
    core.h \
    heightfield.h \
    integerfield.h \
    isolation.h \
    mainwindow.h \
    ppa.h \
    riversnet.h \
    scalarfield.h \
    terrainanalysis.h \
    terrainwidget.h \
    akirmse-mountains/coordinate_system.h \
    akirmse-mountains/divide_tree.h \
    akirmse-mountains/domain_map.h \
    akirmse-mountains/island_tree.h \
    akirmse-mountains/latlng.h \
    akirmse-mountains/line_tree.h \
    akirmse-mountains/pixel_array.h \
    akirmse-mountains/primitives.h \
    akirmse-mountains/tile.h \
    akirmse-mountains/tree_builder.h \
    akirmse-mountains/util.h \
    akirmse-mountains/utm.h \
    akirmse-mountains/utm_coordinate_system.h

FORMS += \
    mainwindow.ui

RESOURCES += \
    resources.qrc

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
