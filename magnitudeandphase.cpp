#include "magnitudeandphase.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <utility>
#include <sstream>
#include <string>
#include <iomanip>
#include "exprtk/exprtk.hpp"

MagnitudeAndPhase::MagnitudeAndPhase(QWidget* parent)
    : QObject(parent), transferFunctionLineEdit(nullptr), sigmaLineEdit(nullptr),
    frequencyMaxlineEdit(nullptr), frequencyMinlineEdit(nullptr),
    _s_real(0.0), _freqMin(0), _freqMax(0), _userFunction("")
{
    // Inicialización por defecto (si es necesario)
}

//contructor out of class
MagnitudeAndPhase::MagnitudeAndPhase(int freqMin, int freqMax, std::string userFunction ,double s_real, QLineEdit *transferFunctionLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent)
    : QObject(parent),
    transferFunctionLineEdit(transferFunctionLineEdit),
    sigmaLineEdit(sigmaLineEdit),
    frequencyMaxlineEdit(frequencyMaxlineEdit),
    frequencyMinlineEdit(frequencyMinlineEdit),
    _s_real(s_real),
    _freqMin(freqMin),
    _freqMax(freqMax),
    _userFunction(userFunction)

{
}
//frequencies function
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;

    double _freqMin = frequencyMinlineEdit->text().toDouble();
    double _freqMax = frequencyMaxlineEdit->text().toDouble();
    double step = 1.0;



    for (int i = 0; _freqMin + i * step <= _freqMax; ++i) {
        double frequency = _freqMin + i * step;
        freqVector.push_back(frequency);
        angularFreqVector.push_back(2 * M_PI * frequency );

        std::cout << "Frequency: " << frequency << " Hz, Angular Frequency: "
                  << angularFreqVector.back() << " rad/s" << std::endl;
    }
    return std::make_pair(angularFreqVector, freqVector);
}
//translatefunction function
std::string MagnitudeAndPhase::translatefunction(double magnitude){
    std::string _userFunction = transferFunctionLineEdit->text().toStdString();
    std::string translatedfunction = _userFunction;
    size_t pos = 0;


    std::string magnitudeStr = std::to_string(magnitude);
    while ((pos = translatedfunction.find("s", pos)) != std::string::npos) {
        translatedfunction.replace(pos, 1, magnitudeStr);
        pos += magnitudeStr.length();
    }
    return translatedfunction;

}

// calculateMagnitude function ->calcula la magnitud /s/ = |j*w|
double MagnitudeAndPhase::calculateMagnitude(double w){
    std::complex<double> s(_s_real, w);
    return std::abs(s);

}
//evaluatetranslatedfunction function
double MagnitudeAndPhase::evaluatetranslatedfunction(std::string translated){
    exprtk::expression<double> expression;
    exprtk::parser<double> parser;

    if (!parser.compile(translated, expression)) {
        qDebug() << "Error al compilar la funcion: "<<parser.error();
        return 0.0;
    }
    double result = expression.value();
    double magnitudedB = 20 *log10(result);
    return magnitudedB;
}
double MagnitudeAndPhase::calculatePhase(double w){
    std::complex<double> s(_s_real, w);
    double magnitude = std::abs(s);
    double phase = std::arg(s);

    qDebug() << "s:" << QString::number(s.real(), 'f', 4) << "+" << QString::number(s.imag(), 'f', 4) << "i -> Magnitud:"
             << QString::number(magnitude, 'f', 4) << ", Fase en (radianes):" << QString::number(phase, 'f', 4);


    return std::arg(s)*(180.0 / M_PI);

}
