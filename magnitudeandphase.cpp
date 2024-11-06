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

std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    // Implementación del método
    std::cout << "Calculando Magnitud y Fase" << std::endl;
}

#include "magnitudeandphase.moc"
