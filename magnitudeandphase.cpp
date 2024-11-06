
// magnitudeandphase.cpp
#include "magnitudeandphase.h"
#include <iostream>

MagnitudeAndPhase::MagnitudeAndPhase(QLineEdit *transferFunctionLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent)
    : QObject(parent),
    transferFunctionLineEdit(transferFunctionLineEdit),
    sigmaLineEdit(sigmaLineEdit),
    frequencyMaxlineEdit(frequencyMaxlineEdit),
    frequencyMinlineEdit(frequencyMinlineEdit)
{
}

void MagnitudeAndPhase::calculate() {
    // Implementación del método
    std::cout << "Calculando Magnitud y Fase" << std::endl;
}

#include "magnitudeandphase.moc"
