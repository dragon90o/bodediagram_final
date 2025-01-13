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
#include <QRect>
#include <QFontMetrics>
#include <QPdfDocument>

Widget::Widget(QWidget *parent)
    : QWidget(parent),
    magnitudeChart(nullptr),
    phaseChart(nullptr),
    ui(new Ui::Widget)


{
    ui->setupUi(this);
    this->setFixedWidth(950);
    this->setFixedHeight(750);
    this->setStyleSheet("background-color: #2B3040; "  );
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
    ui->clearButton->setStyleSheet(
        "background: qlineargradient(x1:0, y1:0, x2:0, y2:1, "
        "stop:0 #635EF2, stop:1 #A260FF); "
        "color: white; "
        "font-size: 16px; "
        "border-radius: 5px;"
        "padding: 10px 20px;"
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
    connect(ui->clearButton, &QPushButton::clicked, this, &Widget::onClearButtonClicked);

}

Widget::~Widget()
{
    delete ui;

}




void Widget::onCalculatePushButtonclicked()
{
    //get input data from the UI
    std::string _numerator = ui->numeratorLineEdit->text().toStdString();
    std::string _denominator = ui->denominatorLineEdit->text().toStdString();
    double _freqMin = ui->frequencyMinlineEdit->text().toDouble();
    double _freqMax = ui->frequencyMaxlineEdit->text().toDouble();
    QString sigmaText = ui->sigmaLineEdit->text();
    double _s_real = sigmaText.isEmpty() ? 0.0 : sigmaText.toDouble();

    // Create an object for Magnitude an Phase calculations
    MagnitudeAndPhase mapObject(_freqMin, _freqMax, _numerator,_denominator, _s_real,
                                ui->numeratorLineEdit,ui->denominatorLineEdit,ui->sigmaLineEdit,
                                ui->frequencyMaxlineEdit, ui->frequencyMinlineEdit,ui->printText,
                                ui->printTextIsStable, this);

    // Process the transfer function and display the result
    std::string result = mapObject.processTransferFunction();
    ui->printText->setPlainText(QString::fromStdString(result));
    bool result1 = mapObject.isStable();
    std::string stabilityResult = result1 ? "Stable: All poles have negative real parts."
                                          : "Unstable: At least one pole has a non-negative real part.";
    ui->printTextIsStable->setPlainText(QString::fromStdString(stabilityResult));

    // Get frequency data
    auto [angularFreqVector, freqVector] = mapObject.frequencies();

     // Create line series for Magnitude and Phase charts
    QLineSeries *magnitudeSeries = new QLineSeries();
    QLineSeries *phaseSeries = new QLineSeries();
    qDebug() << "\nMagnitude and Phase calculations";

    // Loop through angular frequencies and calculate magnitude and phase
    for (const auto& s : angularFreqVector) {
        std::pair<double, std::complex<double>> magnitudeResult = mapObject.calculateMagnitude(s);

        // Skip if magnitude is zero
        if (magnitudeResult.second == std::complex<double>(0, 0)) continue;

        double magnitudeDB = magnitudeResult.first; // Magnitud en dB
        double phaseDegrees = mapObject.calculatePhase(magnitudeResult.second); // Phase in degrees

        //Add data to the series
        magnitudeSeries->append(s, magnitudeDB);
        phaseSeries->append(s, phaseDegrees);
    }


    // Create and set up Magnitude chart
    if (magnitudeChart) {
        delete magnitudeChart;
        magnitudeChart = nullptr;
    }
    magnitudeChart = new QChart();
    magnitudeChart->setTitle("Magnitude Diagram");
    magnitudeChart->addSeries(magnitudeSeries);

    // Configure the X and Y axes for the Magnitude chart
    QLogValueAxis *axisXMag = new QLogValueAxis();
    axisXMag->setTitleText("Frecuency (rad/s)");
    axisXMag->setRange(_freqMin, _freqMax);
    axisXMag->setLabelFormat("%.2f");
    magnitudeChart->addAxis(axisXMag, Qt::AlignBottom);
    magnitudeSeries->attachAxis(axisXMag);


    QValueAxis *axisYMag = new QValueAxis();
    axisYMag->setTitleText("Magnitude (dB)");
    axisYMag->setRange(-30, 20);
    magnitudeChart->addAxis(axisYMag, Qt::AlignLeft);
    magnitudeSeries->attachAxis(axisYMag);

    // Set the chart for the Magnitude graphic view
    ui->magnitudeGraphicView->setChart(magnitudeChart);

    // Create and set up Phase chart
    if (phaseChart) {
        delete phaseChart;
        phaseChart = nullptr;
    }
    phaseChart = new QChart();
    phaseChart->setTitle("Phase Diagram");
    phaseChart->addSeries(phaseSeries); // Apply theme Blue Cerulean

     // Configure the X and Y axes for the Phase chart
    QLogValueAxis *axisXPhase = new QLogValueAxis();
    axisXPhase->setTitleText("Frecuency (rad/s)");
    axisXPhase->setRange(_freqMin, _freqMax); // Ajustar el rango de frecuencia
    axisXPhase->setLabelFormat("%.2f");
    phaseChart->addAxis(axisXPhase, Qt::AlignBottom);
    phaseSeries->attachAxis(axisXPhase);

    QValueAxis *axisYPhase = new QValueAxis();
    axisYPhase->setTitleText("Phase (degrees)");
    axisYPhase->setRange(-360, 360);
    phaseChart->addAxis(axisYPhase, Qt::AlignLeft);
    phaseSeries->attachAxis(axisYPhase);

    ui->phaseGraphicView->setChart(phaseChart);
    qDebug() << "Number of series in Magnitude Chart:" << magnitudeChart->series().size();

    qDebug() << "Number of series in Phase Chart:" << phaseChart->series().size();

}
void Widget::onExportButtonClicked()
{
    // Check if the charts have been generated
    if (!magnitudeChart || !phaseChart) {
        qDebug() << "Las gráficas no se han generado aún.";
        return;
    }

    // Open a dialog for the user to choose the file location
    QString filePath = QFileDialog::getSaveFileName(
        this,
        "Save graphics as .PDF",
        QDir::homePath() + "/bodeDiagram.pdf",
        "PDF Files (*.pdf)"
        );

     // If the user cancels the operation, exit the function
    if (filePath.isEmpty()) {
        qDebug() << "User canceled the save operation.";
        return;
    }

    // Ensure the file path ends with .pdf
    if (!filePath.endsWith(".pdf")) {
        filePath += ".pdf";
    }

    // Set up the QPrinter object for PDF output
    QPrinter printer(QPrinter::HighResolution);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName(filePath);

    // Set the page size and orientation
    printer.setPageSize(QPageSize::A4);
    printer.setPageOrientation(QPageLayout::Portrait);

    QPainter painter;
     // Start painting to the printer
    if (!painter.begin(&printer)) {
        qDebug() << "No se pudo iniciar el pintor para las gráficas.";
        return;
    }

    QRect pageRect = printer.pageRect(QPrinter::DevicePixel).toRect();

    // Draw the magnitude chart manually
    QRect magnitudeRect = QRect(pageRect.left(), pageRect.top(), pageRect.width(), pageRect.height() / 2);
    QPixmap magnitudePixmap(magnitudeChart->size().toSize());
    magnitudePixmap.fill(Qt::transparent);

    QPainter chartPainter(&magnitudePixmap);
    magnitudeChart->scene()->render(&chartPainter, QRectF(magnitudePixmap.rect()), magnitudeChart->sceneBoundingRect());
    chartPainter.end();
    painter.drawPixmap(magnitudeRect, magnitudePixmap);

    // Draw the phase chart manually
    QRect phaseRect = QRect(pageRect.left(), pageRect.top() + pageRect.height() / 2, pageRect.width(), pageRect.height() / 2);
    QPixmap phasePixmap(phaseChart->size().toSize());
    phasePixmap.fill(Qt::transparent);

    QPainter phaseChartPainter(&phasePixmap);
    phaseChart->scene()->render(&phaseChartPainter, QRectF(phasePixmap.rect()), phaseChart->sceneBoundingRect());
    phaseChartPainter.end();
    painter.drawPixmap(phaseRect, phasePixmap);

    painter.end();
    qDebug() << "File successfully saved at:" << filePath;
}
void Widget::onClearButtonClicked()
{
    // Limpiar los campos de entrada
    ui->numeratorLineEdit->clear();
    ui->denominatorLineEdit->clear();
    ui->frequencyMinlineEdit->clear();
    ui->frequencyMaxlineEdit->clear();
    ui->sigmaLineEdit->clear();

    // Limpiar los textos de resultados
    ui->printText->clear();
    ui->printTextIsStable->clear();

    // Crear un nuevo gráfico en lugar de limpiar
    magnitudeChart = new QChart();
    phaseChart = new QChart();
    }





