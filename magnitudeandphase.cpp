#include "magnitudeandphase.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <utility>
#include <sstream>
#include <string>
#include <iomanip>
#include <regex>
#include "exprtk/exprtk.hpp"
#include <Eigen/Dense>


MagnitudeAndPhase::MagnitudeAndPhase(QWidget* parent)
    : QObject(parent), numeratorLineEdit(nullptr), denominatorLineEdit(nullptr), sigmaLineEdit(nullptr),
    frequencyMaxlineEdit(nullptr), frequencyMinlineEdit(nullptr),
    _s_real(0.0), _freqMin(0), _freqMax(0), _numerator(""), _denominator("")
{
    // Inicialización por defecto (si es necesario)
}

//contructor out of class
MagnitudeAndPhase::MagnitudeAndPhase(int freqMin, int freqMax, std::string numerator, std::string denominator ,double s_real, QLineEdit *numeratorLineEdit, QLineEdit *denominatorLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit, QLineEdit *frequencyMinlineEdit, QObject *parent)
    : QObject(parent),
    numeratorLineEdit(numeratorLineEdit),
    denominatorLineEdit(denominatorLineEdit),
    sigmaLineEdit(sigmaLineEdit),
    frequencyMaxlineEdit(frequencyMaxlineEdit),
    frequencyMinlineEdit(frequencyMinlineEdit),
    _s_real(s_real),
    _freqMin(freqMin),
    _freqMax(freqMax),
    _numerator(numerator),
    _denominator(denominator)

{
}

//Method: parseCoefficients
std::vector<double> MagnitudeAndPhase::parseCoefficients(const std::string& polynomial) {
    std::vector<double> coefficients;
    std::string modifiedPoly = polynomial;


    // Use a regular expression to separate the expression and properly handle the operators.
    // Regex to handle coefficients with s^n or s
    std::regex termRegex("([+-]?\\d*\\.?\\d+)?(s(\\^\\d+)?)?");
    auto termsBegin = std::sregex_iterator(modifiedPoly.begin(), modifiedPoly.end(), termRegex);
    auto termsEnd = std::sregex_iterator();
    //DEBUGING -->
    std::cout << "Debugging: Found " << std::distance(termsBegin, termsEnd) << " matches." << std::endl;
    //END DEBBUGING -->
    //  Iterate over terms and extract coefficients
    int maxPower = 0;
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        if (match.str().empty()) continue;

        std::string sTerm = match[2].str();

        if (!sTerm.empty()) { // Check if it's a term with 's'
            if (sTerm.find('^') != std::string::npos) {
                int power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
                maxPower = std::max(maxPower, power);
            } else {
                maxPower = std::max(maxPower, 1); // 's' without '^n' means power 1
            }
        }
    }

    // Initialize the coefficients vector with zeros
    coefficients.resize(maxPower + 1, 0.0);

    // Fill coefficients in the proper positions
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        if (match.str().empty() || (match[1].str().empty() && match[2].str().empty())) continue;

        std::string coefStr = match[1].str();  // Coefficient part
        std::string sTerm = match[2].str();    // 's' term

        // Parse coefficient
        double coefficient = 1.0; // Default is 1
        if (!coefStr.empty()) {
            coefficient = std::stod(coefStr);
        }

        // Determine power of s
        int power = 0;
        if (!sTerm.empty()) {
            if (sTerm.find('^') != std::string::npos) {
                power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
            } else {
                power = 1; // 's' without '^n' means power 1
            }
        }

        // Place coefficient in the correct position
        coefficients[maxPower - power] = coefficient;
        //DEBUGING
        // Debugging: Show each term
        std::cout << "Debugging: term: " << match.str()
                  << ", coefficient: " << coefficient
                  << ", power: " << power
                  << ", maxPower: " << maxPower << std::endl;
    }

    // Debugging: Show the final coefficients vector
    std::cout << "Debugging: Coefficients vector: ";
    for (const auto& coef : coefficients) {
        std::cout << coef << " ";
    }
    std::cout << std::endl;
    //END DEBUGING
    return coefficients;
}

//Method: findRoots
std::vector<std::complex<double>> MagnitudeAndPhase::findRoots(const std::vector<double>& coefficients) {
    size_t degree = coefficients.size() - 1;
    //DEBUGING -->
    std::cout << "DEBUGING METODO FINDROOTS"<< std::endl;
    std::cout << "Debugging: Degree of polynomial: " << degree << std::endl;
    // Mostrar los coeficientes
    std::cout << "Debugging: Coefficients: ";

    for (const auto& coef : coefficients) {
        std::cout << coef << " ";
    }
    std::cout << std::endl;

    // END DEBUGING -->

    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();

    for (size_t i = 0; i < degree; ++i) {
        // Asignar coeficientes en la última columna
        std::cout << "Assigning coefficient[" << degree - i << "] = "
                  << coefficients[degree - i] << " to companion matrix.\n";

        companion(i, degree - 1) = -coefficients[degree - i] / coefficients.front();

        // Asignar las subdiagonales
        if (i < degree - 1) {
            companion(i + 1, i) = 1.0;
        }
    }
    //DEBUGING -->
    // Mostrar la matriz compañera
    std::cout << "Debugging: Companion matrix:\n" << companion << std::endl;
    //END DEBUGING-->
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::VectorXcd roots = solver.eigenvalues();

    // Ajustar las raíces con sigma si es necesario
    std::vector<std::complex<double>> result(roots.size());
    for (size_t i = 0; i < roots.size(); ++i) {
        if (std::abs(_s_real) > 1e-9) { // Si sigma se especificó
            result[i] = roots[i] + _s_real;
        } else { // Sin desplazamiento
            result[i] = roots[i];
        }

        // Debugging: Mostrar cada raíz ajustada
        std::cout << "Root before adjustment: " << roots[i]
                  << ", after adjustment with sigma: " << result[i] << std::endl;
    }
    // Mostrar las raíces ajustadas dependiendo de si sigma fue proporcionado
    std::cout << "Debugging: Roots (final, adjusted or original):" << std::endl;
    for (size_t i = 0; i < roots.size(); ++i) {
        if (std::abs(_s_real) > 1e-9) { // Si sigma fue especificado, mostrar el ajuste
            std::cout << result[i] << std::endl;
        } else { // Si no fue especificado, mostrar la raíz original
            std::cout << roots[i] << std::endl;
        }
    }

    return result;
}
void MagnitudeAndPhase::processTransferFunction() {
    // Convertir cadenas a coeficientes
    std::vector<std::complex<double>> zeros;
    std::vector<std::complex<double>> poles;
    auto numCoefficients = parseCoefficients(_numerator);
    auto denCoefficients = parseCoefficients(_denominator);

    // Calculate zeros and poles
    zeros = findRoots(numCoefficients);
    poles = findRoots(denCoefficients);

    //Print results
    std::cout << "Zeros (Roots of the Numerator): ";
    for (const auto& zero : zeros) {
        std::cout << zero << " ";
    }
    std::cout << std::endl;

    std::cout << "Poles (Roots of the Denominator): ";
    for (const auto& pole : poles) {
        std::cout << pole << " ";
    }
    std::cout << std::endl;
}
//Method: frequencies
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;
    //DEBUGING
    std::cout << "DEBUGING FREQUENCIES"<< std::endl;
        // Debug: Verificar límites
    std::cout << "Debugging: Frequencies method" << std::endl;
    std::cout << "Debugging: Frequency range - Min: " << _freqMin << ", Max: " << _freqMax << std::endl;
    //END DEBUGING
    for (double i = _freqMin; i <= _freqMax; i += 1.0) {
        freqVector.push_back(i);
        // Angular frequency calculation
        angularFreqVector.push_back(2 * M_PI * i);
        //DEBUGING
        // Debug: Mostrar valores en cada iteración
        std::cout << "Frequency: " << i << ", Angular Frequency: " << angularFreqVector.back() << "rad/s."<< std::endl;
        //END DEBUGING
    }
    //DEBUGING
    // Debug: Mostrar tamaños de los vectores
    std::cout << "Debugging: Frequency vector size: " << freqVector.size() << std::endl;
    std::cout << "Debugging: Angular Frequency vector size: " << angularFreqVector.size() << std::endl;
    //END DEBUGING
    // Returns the pair of vectors
    return std::make_pair(angularFreqVector, freqVector);
}

// Method: translateFunction
std::complex<double> MagnitudeAndPhase::translateFunction(double angularFrequency, bool isNumerator) {
    std::string functionStr = isNumerator ? _numerator : _denominator;
    std::complex<double> jw(0.0, angularFrequency); // j * omega
    std::complex<double> s = (std::abs(_s_real) > 1e-9)? std::complex<double>(_s_real, angularFrequency): jw;

    // Evaluar cada término en la expresión
    std::regex termRegex("([+-]?\\d*\\.?\\d+)?(s(\\^\\d+)?)?");
    auto termsBegin = std::sregex_iterator(functionStr.begin(), functionStr.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    std::complex<double> result(0.0, 0.0); // Resultado acumulado

    // Mostrar la función a evaluar
    std::cout << "Debugging: Evaluating function: " << functionStr << std::endl;
    std::cout << "Debugging: Angular frequency (w): " << angularFrequency << std::endl;

    for (auto it = termsBegin; it != termsEnd; ++it) {
        std::smatch match = *it;
        if (match.str().empty()) continue;

        // Obtener coeficiente
        double coefficient = 1.0; // Por defecto
        if (!match[1].str().empty()) {
            coefficient = std::stod(match[1].str());
        }

        // Obtener potencia de s
        int power = 0;
        if (!match[2].str().empty()) {
            if (match[2].str().find('^') != std::string::npos) {
                power = std::stoi(match[2].str().substr(match[2].str().find('^') + 1));
            } else {
                power = 1;
            }
        }

        // Evaluar el término
        std::complex<double> termValue = coefficient * std::pow(s, power);
        result += termValue;

        // Mostrar información de depuración por cada término
        std::cout << "Debugging: Term: " << match.str()
                  << ", Coefficient: " << coefficient
                  << ", Power: " << power
                  << ", Term Value: " << termValue
                  << ", Accumulated Result: " << result
                  << std::endl;
    }

    // Mostrar resultado final
    std::cout << "Debugging: Final translated function result: " << result << std::endl;

    return result;
}

// Method:calculateMagnitude
std::pair<double, std::complex<double>> MagnitudeAndPhase::calculateMagnitude(double angularFrequency)  {
    // Calcular el numerador y denominador traducidos
    // Ajuste de la frecuencia compleja s = sigma + jw
    std::complex<double> s = _s_real + std::complex<double>(0, angularFrequency);

    // Mostrar s en consola para debugging
    std::cout << "Debugging: s = sigma + jw = " << s << std::endl;
    std::cout << "Debugging: Magnitude of s = |s| = " << std::abs(s) << std::endl;

    std::complex<double> resultNumerator = translateFunction(angularFrequency, true);
    std::complex<double> resultDenominator = translateFunction(angularFrequency, false);

    // Evitar división por un denominador cercano a cero
    if (std::abs(resultDenominator) < 1e-12) {
        std::cerr << "Error: Denominator is too close to zero!" << std::endl;
        return {};
    }


    // Calcular la función de transferencia (numerador / denominador)
    std::complex<double> transferFunction = resultNumerator / resultDenominator;

    // Calcular magnitud en decibeles (dB)
    double magnitudeDB = 20 * std::log10(std::abs(transferFunction));

    // Mostrar resultados
    std::cout << "Numerator: " << resultNumerator << std::endl;
    std::cout << "Denominator: " << resultDenominator << std::endl;
    std::cout << "Transfer Function: " << transferFunction << std::endl;
    std::cout << "Magnitude (dB): " << magnitudeDB << std::endl;
    // Llamar a la función calculatePhase con la transferencia calculada

    return {magnitudeDB, transferFunction};
}

// Method: calculatePhase
double MagnitudeAndPhase::calculatePhase(const std::complex<double>& transferFunction){
    double phaseRadians = std::atan2(transferFunction.imag(), transferFunction.real());

    // Convertir la fase a grados
    double phaseDegrees = phaseRadians * (180.0 / M_PI);

    std::cout << "DEBUGING PHASE"<<std::endl;
    // Debugging para asegurarse de obtener los resultados esperados
    std::cout << "Debugging: Transfer Function = " << transferFunction << std::endl;
    std::cout << "Debugging: Real = " << transferFunction.real() << std::endl;
    std::cout << "Debugging: Imaginary = " << transferFunction.imag() << std::endl;
    std::cout << "Debugging: Phase (radians) = " << phaseRadians << std::endl;
    std::cout << "Debugging: Phase (degrees) = " << phaseDegrees << std::endl;

    return phaseDegrees;
}





