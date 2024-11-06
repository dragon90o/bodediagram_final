
// magnitudeandphase.h
#ifndef MAGNITUDEANDPHASE_H
#define MAGNITUDEANDPHASE_H

#include <QLineEdit>
#include <QObject>
#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <utility>
#include <iomanip>
#include "exprtk/exprtk.hpp"

    class MagnitudeAndPhase : public QObject
{
    Q_OBJECT

public:
    explicit MagnitudeAndPhase(QLineEdit *transferFunctionLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent = nullptr);

    void calculate();  // Método público para la lógica de cálculo

private:
    QLineEdit *transferFunctionLineEdit;
    QLineEdit *sigmaLineEdit;
    QLineEdit *frequencyMaxlineEdit;
    QLineEdit *frequencyMinlineEdit;
};

#endif // MAGNITUDEANDPHASE_H

