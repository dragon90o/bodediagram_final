#ifndef MAGNITUDEANDPHASE_H
#define MAGNITUDEANDPHASE_H

#include <iostream>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <utility>
#include <iomanip>
#include "exprtk/exprtk.hpp"

class MagnitudeAndPhase
{
    public:
   explicit MagnitudeAndPhase(QLineEdit *transferFunctionLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent = nullptr);


    private:


};

#endif // MAGNITUDEANDPHASE_H
