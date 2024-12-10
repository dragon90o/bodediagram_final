QT       += core gui serialport widgets printsupport charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17 bigobj
QMAKE_CXXFLAGS += -bigobj

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    magnitudeandphase.cpp \
    main.cpp \
    widget.cpp \


HEADERS += \
    magnitudeandphase.h \
    widget.h


FORMS += \
    widget.ui

INCLUDEPATH += \
    $$PWD/exprtk \
    $$PWD/eigen-3.4.0/eigen-3.4.0

MOC_DIR = $$PWD/build/Desktop_Qt_6_7_3_MSVC2019_64bit-Debug/moc

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target


