#include "magnitudeandphase.h"
#include <iostream>
#include <cmath>

MagnitudeAndPhase::MagnitudeAndPhase(QLineEdit *transferFunctionLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent)
    : QObject(parent),
    transferFunctionLineEdit(transferFunctionLineEdit),
    sigmaLineEdit(sigmaLineEdit),
    frequencyMaxlineEdit(frequencyMaxlineEdit),
    frequencyMinlineEdit(frequencyMinlineEdit)
{
}

std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;

    double _freqMin = frequencyMinlineEdit->text().toDouble();
    double _freqMax = frequencyMaxlineEdit->text().toDouble();
    double step = 1.0;



    for (int i = 0; _freqMin + i * step <= _freqMax; ++i) {
        double frequency = _freqMin + i * step;
        freqVector.push_back(i);
        angularFreqVector.push_back(2 * M_PI * i );

        std::cout << "Frequency: " << frequency << " Hz, Angular Frequency: "
                  << angularFreqVector.back() << " rad/s" << std::endl;
    }
    return std::make_pair(angularFreqVector, freqVector);
}


