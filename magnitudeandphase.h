#ifndef MAGNITUDEANDPHASE_H
#define MAGNITUDEANDPHASE_H

#include <QLineEdit>
#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <utility>
#include <iomanip>
#include <regex>
#include <QObject>
#include <Eigen/Dense>
#include <QTextEdit>


class MagnitudeAndPhase : public QObject
{
    Q_OBJECT

public:
    MagnitudeAndPhase(QWidget* parent = nullptr);

    explicit MagnitudeAndPhase(double freqMin, double freqMax, std::string numerator, std::string denominator,
                               double s_real, QLineEdit *numeratorLineEdit,
                               QLineEdit *denominatorLineEdit, QLineEdit *sigmaLineEdit,
                               QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit,
                               QTextEdit *printText, QTextEdit *printTextIsStable, QObject *parent = nullptr);

    std::string processTransferFunction();
    std::vector<double> parseCoefficients(const std :: string& polynomial);
    std::vector<std::complex<double>> findRoots(const std::vector<double>& coefficients);
    std::pair<std::vector<double>, std::vector<double>> frequencies();  // Método público para la lógica de cálculo
    std::complex<double> translateFunction(double angularFrequency, bool isNumerator);
    std::pair<double, std::complex<double>> calculateMagnitude(double angularFrequency);
    double calculatePhase(const std::complex<double>& transferFunction);
    bool isStable();




private:
    QLineEdit *numeratorLineEdit;
    QLineEdit *denominatorLineEdit;
    QLineEdit *sigmaLineEdit;
    QLineEdit *frequencyMaxlineEdit;
    QLineEdit *frequencyMinlineEdit;
    double _s_real;
    double _freqMin;
    double _freqMax;
    std::string _numerator;
    std::string _denominator;
    QTextEdit *printText;
    QTextEdit *printTextIsStable;



};

#endif // MAGNITUDEANDPHASE_H
