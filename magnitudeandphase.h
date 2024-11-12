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
    explicit MagnitudeAndPhase(QLineEdit *transferFunctionLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent = nullptr);

    std::pair<std::vector<double>, std::vector<double>> frequencies();  // Método público para la lógica de cálculo

private:
    QLineEdit *transferFunctionLineEdit;
    QLineEdit *sigmaLineEdit;
    QLineEdit *frequencyMaxlineEdit;
    QLineEdit *frequencyMinlineEdit;
};

#endif // MAGNITUDEANDPHASE_H

