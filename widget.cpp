#include "widget.h"
#include "ui_widget.h"
#include "magnitudeandphase.h"
#include <QDebug>
#include <iomanip>

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)

{
    ui->setupUi(this);
    connect(ui->calculatePushButton, &QPushButton::clicked, this, &Widget::onCalculatePushButtonclicked);
}

Widget::~Widget()
{
    delete ui;
}

void Widget::onCalculatePushButtonclicked()
{
    std::string _userFunction = ui->transferFunctionlineEdit->text().toStdString();
    double _freqMin = ui->frequencyMinlineEdit->text().toDouble();
    double _freqMax = ui->frequencyMaxlineEdit->text().toDouble();
    double _s_real = ui->sigmaLineEdit->text().toDouble();

    MagnitudeAndPhase mapObject(_freqMin, _freqMax, _userFunction, _s_real,
                                ui->transferFunctionlineEdit, ui->sigmaLineEdit,
                                ui->frequencyMaxlineEdit, ui->frequencyMinlineEdit, this);


    auto [angularFreqVector, freqVector] = mapObject.frequencies();

    qDebug() << "frequencias angulares Totales (w): ";
    for (const auto& w :angularFreqVector) {
        qDebug()<< QString::number(w, 'f', 4);
    }

    qDebug() << "frequencias Totales (Hz): ";
        for (const auto& freq : freqVector) {
        qDebug()<<QString::number(freq, 'f', 4);
    }

    for(const auto& w : angularFreqVector){
        double magnitude = mapObject.calculateMagnitude(w);
        double phase = mapObject.calculatePhase(w);

        qDebug() <<"--> Magnitud: "<<QString::number(magnitude, 'f', 4);
        qDebug() <<"--> Fase Total en (grados): "<< QString::number(phase, 'f', 4);



        qDebug() <<"formula (despues de traducir): ";
        std::string translated = mapObject.translatefunction(magnitude);
        qDebug() << QString::fromStdString(translated);

        double magnitudedB = mapObject.evaluatetranslatedfunction(translated);
        qDebug() << "Magnitud total en (dB): "<< QString::number(magnitudedB, 'f', 4);

    }


}
