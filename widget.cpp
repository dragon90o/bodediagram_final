#include "widget.h"
#include "ui_widget.h"
#include "magnitudeandphase.h"
#include <QDebug>
#include <iomanip>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QLogValueAxis>
#include <QPrinter>
#include <QPrintDialog>
#include <QFileDialog>
#include <QPainter>

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)


{
    ui->setupUi(this);
    this->setFixedWidth(950);
    this->setFixedHeight(750);
    this->setStyleSheet("background-color: #2B3040; "
                        //"border-radius: 5px;"

                        );
    ui->frequencyMaxlineEdit->setStyleSheet("background-color: white; color: black");
    ui->frequencyMinlineEdit->setStyleSheet("background-color: white; color: black");
    ui->sigmaLineEdit->setStyleSheet("background-color: white; color: black");
    ui->denominatorLineEdit->setStyleSheet("background-color: white; color: black");
    ui->numeratorLineEdit->setStyleSheet("background-color: white; color: black");
    ui->calculatePushButton->setStyleSheet(
        "background: qlineargradient(x1:0, y1:0, x2:0, y2:1, "
        "stop:0 #635EF2, stop:1 #A260FF); "
        "color: #CBC5D9; "
        "font-size: 18px; "
        "border-radius: 5px;"
        );
    ui->calculatePushButton->setFixedHeight(40);
    ui->exportButton->setStyleSheet(
        "background: qlineargradient(x1:0, y1:0, x2:0, y2:1, "
        "stop:0 #635EF2, stop:1 #A260FF); "
        "color: #CBC5D9; "
        "font-size: 18px; "
        "border-radius: 5px;"
        );
    ui->exportButton->setFixedHeight(40);
    ui->printText->setStyleSheet("background-color: white; color: black");
    ui->printTextIsStable->setStyleSheet("background-color: white; color: black");
    ui->frame_2->setStyleSheet("background-color: #252526");
    ui->groupBox_4->setStyleSheet("color: white; font-size: 14px");
    ui->groupBox_3->setStyleSheet("color: white; font-size: 14px");
    ui->groupBox_2->setStyleSheet("color: white; font-size: 14px");
    ui->groupBox->setStyleSheet("color: white; font-size: 14px");
    connect(ui->calculatePushButton, &QPushButton::clicked, this, &Widget::onCalculatePushButtonclicked);
    connect(ui->exportButton, &QPushButton::clicked, this, &Widget::onExportButtonClicked);

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
                                ui->frequencyMaxlineEdit, ui->frequencyMinlineEdit,ui->printText,
                                ui->printTextIsStable, this);


    std::string result = mapObject.processTransferFunction();
    ui->printText->setPlainText(QString::fromStdString(result));
    bool result1 = mapObject.isStable();
    std::string stabilityResult = result1 ? "Stable: All poles have negative real parts."
                                          : "Unstable: At least one pole has a non-negative real part.";
    ui->printTextIsStable->setPlainText(QString::fromStdString(stabilityResult));

    auto [angularFreqVector, freqVector] = mapObject.frequencies();

    // Crear las series gráficas
    QLineSeries *magnitudeSeries = new QLineSeries();
    QLineSeries *phaseSeries = new QLineSeries();
    qDebug() << "\nMagnitude and Phase calculations";

    for (const auto& s : angularFreqVector) {
        // Calcular la magnitud y la función de transferencia
        std::pair<double, std::complex<double>> magnitudeResult = mapObject.calculateMagnitude(s);

        if (magnitudeResult.second == std::complex<double>(0, 0)) continue;

        double magnitudeDB = magnitudeResult.first; // Magnitud en dB
        double phaseDegrees = mapObject.calculatePhase(magnitudeResult.second); // Calcular fase

        // Agregar datos a las series gráficas
        magnitudeSeries->append(s, magnitudeDB);
        phaseSeries->append(s, phaseDegrees);
    }
   // qDebug() << "\nTotal angular frequencies (w) in radians:";
    //for (const auto& w : angularFreqVector) {
     //   qDebug() << QString::number(w, 'f', 4) << "rad";
    //}

    //qDebug() << "\nTotal frequencies (Hz):";
    //for (const auto& freq : freqVector) {
       // qDebug() << QString::number(freq, 'f', 4) << "Hz";
    //}


    // Crear y configurar el gráfico de magnitud
    if (magnitudeChart) {
        delete magnitudeChart; // Limpia la instancia previa para evitar fugas
        magnitudeChart = nullptr;
    }
    magnitudeChart = new QChart();
    magnitudeChart->setTitle("Diagrama de Magnitud");
    magnitudeChart->addSeries(magnitudeSeries);

    QLogValueAxis *axisXMag = new QLogValueAxis();
    axisXMag->setTitleText("Frecuencia (rad/s)");
    axisXMag->setRange(_freqMin, _freqMax); // Ajustar el rango de frecuencia
    axisXMag->setLabelFormat("%.2f");
    magnitudeChart->addAxis(axisXMag, Qt::AlignBottom);
    magnitudeSeries->attachAxis(axisXMag);


    QValueAxis *axisYMag = new QValueAxis();
    axisYMag->setTitleText("Magnitud (dB)");
    axisYMag->setRange(-100, 100);
    magnitudeChart->addAxis(axisYMag, Qt::AlignLeft);
    magnitudeSeries->attachAxis(axisYMag);

    ui->magnitudeGraphicView->setChart(magnitudeChart);

    // Crear y configurar el gráfico de fase
    if (phaseChart) {
        delete phaseChart;
        phaseChart = nullptr;
    }
    phaseChart = new QChart();
    phaseChart->setTitle("Diagrama de Fase");
    phaseChart->addSeries(phaseSeries);

    QLogValueAxis *axisXPhase = new QLogValueAxis();
    axisXPhase->setTitleText("Frecuencia (rad/s)");
    axisXPhase->setRange(_freqMin, _freqMax); // Ajustar el rango de frecuencia
    axisXPhase->setLabelFormat("%.2f");
    phaseChart->addAxis(axisXPhase, Qt::AlignBottom);
    phaseSeries->attachAxis(axisXPhase);

    QValueAxis *axisYPhase = new QValueAxis();
    axisYPhase->setTitleText("Fase (grados)");
    axisYPhase->setRange(-180, 180);
    phaseChart->addAxis(axisYPhase, Qt::AlignLeft);
    phaseSeries->attachAxis(axisYPhase);

    ui->phaseGraphicView->setChart(phaseChart);
    qDebug() << "Cantidad de series en Magnitude Chart:" << magnitudeChart->series().size();

    qDebug() << "Cantidad de series en Phase Chart:" << phaseChart->series().size();

}
void Widget::onExportButtonClicked()
{
    if (!magnitudeChart || !phaseChart) {
        qDebug() << "Las gráficas no se han generado aún.";
        return;
    }

    QString filePath = QFileDialog::getSaveFileName(
        this,
        "Guardar Gráficas como PDF",
        QDir::homePath() + "/bodeDiagram.pdf",
        "PDF Files (*.pdf)"
        );

    if (filePath.isEmpty()) {
        qDebug() << "El usuario canceló la operación de guardar.";
        return;
    }

    if (!filePath.endsWith(".pdf")) {
        filePath += ".pdf";
    }

    QPrinter printer(QPrinter::HighResolution);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName(filePath);

    QPainter painter;

    if (!painter.begin(&printer)) {
        qDebug() << "No se pudo iniciar el pintor para las gráficas.";
        return;
    }

    // Ajustar el espacio para renderizar las gráficas
    int margin = 30; // Espacio adicional para márgenes
    int halfPageHeight = (printer.pageRect(QPrinter::Point).height() - 3 * margin) / 2; // Dividir espacio entre dos gráficas

    QRect magnitudeRect = QRect(
        margin,
        margin,
        printer.pageRect(QPrinter::Point).width() - 2 * margin,
        halfPageHeight
        );

    QRect phaseRect = QRect(
        margin,
        margin + halfPageHeight + margin,
        printer.pageRect(QPrinter::Point).width() - 2 * margin,
        halfPageHeight
        );

    // Renderizar la gráfica de magnitud
    magnitudeChart->scene()->render(&painter, magnitudeRect);
    qDebug() << "Gráfica de magnitud renderizada.";

    // Crear una nueva página para la fase
    if (!printer.newPage()) {
        qDebug() << "No se pudo crear una nueva página para la fase.";
    }

    // Renderizar la gráfica de fase
    phaseChart->scene()->render(&painter, phaseRect);
    qDebug() << "Gráfica de fase renderizada.";

    painter.end();

    qDebug() << "Archivo guardado exitosamente en:" << filePath;
}
