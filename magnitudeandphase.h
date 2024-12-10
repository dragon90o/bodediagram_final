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
#include "exprtk/exprtk.hpp"
#include <QObject>
#include <Eigen/Dense>

    class MagnitudeAndPhase : public QObject
{
    Q_OBJECT

public:
    MagnitudeAndPhase(QWidget* parent = nullptr);

    explicit MagnitudeAndPhase(int freqMin, int freqMax, std::string numerator, std::string denominator ,double s_real, QLineEdit *numeratorLineEdit, QLineEdit *denominatorLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent = nullptr);

    std::vector<double> parseCoefficients(const std :: string& polynomial);
    std::vector<std::complex<double>> findRoots(const std::vector<double>& coefficients);
    std::pair<std::vector<double>, std::vector<double>> frequencies();  // Método público para la lógica de cálculo
    std::complex<double> translateFunction(double angularFrequency, bool isNumerator);
    std::pair<double, std::complex<double>> calculateMagnitude(double angularFrequency);
    double calculatePhase(const std::complex<double>& transferFunction);
    void processTransferFunction();


private:
    QLineEdit *numeratorLineEdit;
    QLineEdit *denominatorLineEdit;
    QLineEdit *sigmaLineEdit;
    QLineEdit *frequencyMaxlineEdit;
    QLineEdit *frequencyMinlineEdit;
    double _s_real;
    int _freqMin;
    int _freqMax;
    std::string _numerator;
    std::string _denominator;
};

#endif // MAGNITUDEANDPHASE_H

