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
#include "exprtk/exprtk.hpp"
#include <QObject>

    class MagnitudeAndPhase : public QObject
{
    Q_OBJECT

public:
    MagnitudeAndPhase(QWidget* parent = nullptr);

    explicit MagnitudeAndPhase(int freqMin, int freqMax, std::string userFunction ,double s_real, QLineEdit *transferFunctionLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent = nullptr);

    std::pair<std::vector<double>, std::vector<double>> frequencies();  // Método público para la lógica de cálculo
    std::string translatefunction(double magnitude);
    double evaluatetranslatedfunction(std::string translated);
    double calculateMagnitude(double w);
    double calculatePhase(double w);

private:
    QLineEdit *transferFunctionLineEdit;
    QLineEdit *sigmaLineEdit;
    QLineEdit *frequencyMaxlineEdit;
    QLineEdit *frequencyMinlineEdit;
    double _s_real;
    int _freqMin;
    int _freqMax;
    std::string _userFunction;
};

#endif // MAGNITUDEANDPHASE_H

