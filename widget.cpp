#include "widget.h"
#include "ui_widget.h"
#include "magnitudeandphase.h"
#include <QDebug>
#include <iomanip>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>


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

    std::string _numerator = ui->numeratorLineEdit->text().toStdString();
    std::string _denominator = ui->denominatorLineEdit->text().toStdString();
    double _freqMin = ui->frequencyMinlineEdit->text().toDouble();
    double _freqMax = ui->frequencyMaxlineEdit->text().toDouble();
    QString sigmaText = ui->sigmaLineEdit->text();
    double _s_real = sigmaText.isEmpty() ? 0.0 : sigmaText.toDouble();

    MagnitudeAndPhase mapObject(_freqMin, _freqMax, _numerator,_denominator, _s_real,
                                ui->numeratorLineEdit,ui->denominatorLineEdit,ui->sigmaLineEdit,
                                ui->frequencyMaxlineEdit, ui->frequencyMinlineEdit, this);


    mapObject.processTransferFunction();

    auto [angularFreqVector, freqVector] = mapObject.frequencies();

    qDebug() << "\nTotal angular frequencies (w) in radians:";
    for (const auto& w : angularFreqVector) {
        qDebug() << QString::number(w, 'f', 4) << "rad";
    }

    qDebug() << "\nTotal frequencies (Hz):";
    for (const auto& freq : freqVector) {
        qDebug() << QString::number(freq, 'f', 4) << "Hz";
    }
    // Crear las series gráficas
    QLineSeries *magnitudeSeries = new QLineSeries();
    QLineSeries *phaseSeries = new QLineSeries();
    qDebug() << "\nMagnitude and Phase calculations";

    for (const auto& s : angularFreqVector) {
        // Calcular la magnitud y la función de transferencia
        std::pair<double, std::complex<double>> magnitudeResult = mapObject.calculateMagnitude(s);

        if (magnitudeResult.second == std::complex<double>(0, 0)) continue; // Prevenir división por cero

        double magnitudeDB = magnitudeResult.first; // Magnitud en dB
        double phaseDegrees = mapObject.calculatePhase(magnitudeResult.second); // Calcular fase

        // Agregar datos a las series gráficas
        magnitudeSeries->append(s, magnitudeDB);
        phaseSeries->append(s, phaseDegrees);
    }

    // Crear y configurar el gráfico de magnitud
    QChart *magnitudeChart = new QChart();
    magnitudeChart->setTitle("Diagrama de Magnitud");
    magnitudeChart->addSeries(magnitudeSeries);

    QValueAxis *axisXMag = new QValueAxis();
    axisXMag->setTitleText("Frecuencia (rad/s)");
    axisXMag->setLabelFormat("%.1f");
    magnitudeChart->addAxis(axisXMag, Qt::AlignBottom);
    magnitudeSeries->attachAxis(axisXMag);

    QValueAxis *axisYMag = new QValueAxis();
    axisYMag->setTitleText("Magnitud (dB)");
    magnitudeChart->addAxis(axisYMag, Qt::AlignLeft);
    magnitudeSeries->attachAxis(axisYMag);

    ui->magnitudeGraphicView->setChart(magnitudeChart);

    // Crear y configurar el gráfico de fase
    QChart *phaseChart = new QChart();
    phaseChart->setTitle("Diagrama de Fase");
    phaseChart->addSeries(phaseSeries);

    QValueAxis *axisXPhase = new QValueAxis();
    axisXPhase->setTitleText("Frecuencia (rad/s)");
    axisXPhase->setLabelFormat("%.1f");
    phaseChart->addAxis(axisXPhase, Qt::AlignBottom);
    phaseSeries->attachAxis(axisXPhase);

    QValueAxis *axisYPhase = new QValueAxis();
    axisYPhase->setTitleText("Fase (grados)");
    phaseChart->addAxis(axisYPhase, Qt::AlignLeft);
    phaseSeries->attachAxis(axisYPhase);

    ui->phaseGraphicView->setChart(phaseChart);


}
