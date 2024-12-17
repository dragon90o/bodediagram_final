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
#include <QDebug>
#include <QTextEdit>



MagnitudeAndPhase::MagnitudeAndPhase(QWidget* parent)
    : QObject(parent), numeratorLineEdit(nullptr), denominatorLineEdit(nullptr), sigmaLineEdit(nullptr),
    frequencyMaxlineEdit(nullptr), frequencyMinlineEdit(nullptr),
    _s_real(0.0), _freqMin(0), _freqMax(0), _numerator(""), _denominator("")
{
    // Inicialización por defecto (si es necesario)
}

//contructor out of class
MagnitudeAndPhase::MagnitudeAndPhase(int freqMin, int freqMax, std::string numerator, std::string denominator ,double s_real, QLineEdit *numeratorLineEdit, QLineEdit *denominatorLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit,
                                     QLineEdit *frequencyMinlineEdit, QTextEdit *printText, QTextEdit *printTextIsStable, QObject *parent)
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
    _denominator(denominator),
    printText(printText),
    printTextIsStable(printTextIsStable)


{
}
//Method: processTransferFunction
std::string MagnitudeAndPhase::processTransferFunction() {
    // Convertir cadenas a coeficientes
    qDebug()<< "Debugging Method processTransferFunction: ";

    std::vector<std::complex<double>> poles;
    std::vector<std::complex<double>> zeros;

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

    // Constuct result
    std::ostringstream result;
    result << "Zeros (Roots of the Numerator): ";
    for (const auto& zero : zeros) {
        result << zero << " ";
    }
    result << "\nPoles (Roots of the Denominator): ";
    for (const auto& pole : poles) {
        result << pole << " ";

    }
    return result.str();
}


//Method: parseCoefficients
std::vector<double> MagnitudeAndPhase::parseCoefficients(const std::string& polynomial) {
    std::vector<double> coefficients;
    std::string modifiedPoly = polynomial;

    qDebug() << "Debugging Method parseCoefficients: ";

    std::regex termRegex(R"(([+-]?\s*\d*\.?\d*)\s*(s(\^\d+)?)?)");
    auto termsBegin = std::sregex_iterator(modifiedPoly.begin(), modifiedPoly.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    int maxPower = 0;
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        if (match.str().empty()) continue;

        std::string sTerm = match[2].str();
        if (!sTerm.empty()) {
            if (sTerm.find('^') != std::string::npos) {
                int power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
                maxPower = std::max(maxPower, power);
            } else {
                maxPower = std::max(maxPower, 1);
            }
        }
    }

    coefficients.resize(maxPower + 1, 0.0);

    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        if (match.str().empty() || (match[1].str().empty() && match[2].str().empty())) continue;

        std::string coefStr = match[1].str();
        std::string sTerm = match[2].str();

        // Eliminar espacios en blanco
        coefStr.erase(remove_if(coefStr.begin(), coefStr.end(), ::isspace), coefStr.end());

        // Parsear coeficiente
        double coefficient = 1.0;
        if (!coefStr.empty() && coefStr != "+" && coefStr != "-") {
            coefficient = std::stod(coefStr);
        } else if (coefStr == "-") {
            coefficient = -1.0;
        }

        // Determinar potencia de s
        int power = 0;
        if (!sTerm.empty()) {
            if (sTerm.find('^') != std::string::npos) {
                power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
            } else {
                power = 1;
            }
        }

        // Colocar coeficiente en la posición correcta
        coefficients[maxPower - power] = coefficient;

        // debugMessage
        std::ostringstream debugMessage;
        debugMessage << "Debugging: term: " << match.str()
                     << ", coefficient: " << coefficient
                     << ", power: " << power
                     << ", maxPower: " << maxPower;
        debugMessage << "\n";
        qDebug() << QString::fromStdString(debugMessage.str());
    }

    //debug Message1
    std::ostringstream debugMessage1;
    debugMessage1 << "Debugging: Coefficients vector: ";
    for (const auto& coef : coefficients) {
        debugMessage1 << coef << " ";
    }
    debugMessage1 << "\n";
    qDebug() << QString::fromStdString(debugMessage1.str());

    return coefficients;
}

//Method: findRoots
std::vector<std::complex<double>> MagnitudeAndPhase::findRoots(const std::vector<double>& coefficients) {
    size_t degree = coefficients.size() - 1;
    qDebug() << "Debugging Method findRoots: ";

    // Validar el grado del polinomio y los coeficientes
    if (degree < 1) throw std::invalid_argument("The polynomial degree must be at least 1.");
    if (coefficients[0] == 0) throw std::invalid_argument("The leading coefficient cannot be zero.");
    if (coefficients.size() != degree + 1) throw std::invalid_argument("The size of the coefficients vector must be degree + 1.");

    //debug Message
    std::ostringstream debugMessage;
    debugMessage << "Debugging: Degree of polynomial: " << degree << std::endl
                 << "Debugging: Coefficients: ";
    for (const auto& coef : coefficients) {
        debugMessage << coef << " ";
    }
    debugMessage << "\n";
    qDebug() << QString::fromStdString(debugMessage.str());


    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();


    // Asignar 1 a las subdiagonales
    for (std::size_t i = 1; i < degree; ++i) {
        companion( i - 1 , i ) = 1.0;
        std::cout << "Subdiagonal: companion(" << (i - 1) << ", " << i << ") = 1.0" << std::endl;
    }
    // Asignar los coeficientes correctamente a la última fila de la matriz
    for (std::size_t i = 0; i < degree; ++i) {
        companion(degree - 1, i) = -coefficients[degree - i] / coefficients[0];
        std::cout << "Last row: companion(" << (degree - 1) << ", "
                  << i << ") = "
                  << -coefficients[degree - i] / coefficients[0] << std::endl;
    }

    debugMessage.str("");
    debugMessage << "Debugging: Companion matrix:\n" << companion << std::endl;
    qDebug() << QString::fromStdString(debugMessage.str());

    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::VectorXcd roots = solver.eigenvalues();

    std::vector<std::complex<double>> result(roots.size());
    for (std::size_t i = 0; i < roots.size(); ++i) {
        result[i] = (std::abs(_s_real) > 1e-9) ? roots[i] + _s_real : roots[i];
        debugMessage.str("");
        debugMessage << "Root before adjustment: " << roots[i] << ", after adjustment with sigma: " << result[i] << std::endl;
        qDebug() << QString::fromStdString(debugMessage.str());
    }

    qDebug() << "Debugging: Roots (final, adjusted or original): " ;
    for (const auto& root : result) {
        std::cout << root << std::endl;
    }

    return result;
}
//Method: frequencies
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;
    qDebug() << "Debugging Method Frequencies: " ;

    //debug Message
    std::ostringstream debugMessage;
    debugMessage << "Debugging: Frequencies method" << std::endl
                 << "Debugging: Frequency range - Min: " << _freqMin << ", Max: " << _freqMax << std::endl;
    qDebug() << QString::fromStdString(debugMessage.str());


    double freq = _freqMin;
    while (freq <= _freqMax) {
        freqVector.push_back(freq);
        // Angular frequency calculation
        angularFreqVector.push_back(2 * M_PI * freq);


        //debugging Message1
        std::ostringstream debugMessage1;
        debugMessage1 << "Frequency: " << freq << ", Angular Frequency: " << angularFreqVector.back() << "rad/s."<< std::endl;
        qDebug() << QString::fromStdString(debugMessage1.str());

        freq *=2.0;

    }
    //debugging Message2
    std::ostringstream debugMessage2;
    debugMessage2 << "Debugging: Frequency vector size: " << freqVector.size() << std::endl
                  << "Debugging: Angular Frequency vector size: " << angularFreqVector.size() << std::endl;
    qDebug() << QString::fromStdString(debugMessage2.str());


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
    qDebug() << "\033[32mDebugging Method translateFunction: \033[0m";

    //debug Message
    std::ostringstream debugMessage;
    debugMessage << "Debugging: Evaluating function: " << functionStr << std::endl                 << "Debugging: Angular frequency (w): " << angularFrequency << std::endl;
    qDebug() << QString::fromStdString(debugMessage.str());

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
        //debug Message1
        std::ostringstream debugMessage1;
        debugMessage1 << "Debugging: Term: " << match.str()
                      << ", Coefficient: " << coefficient
                      << ", Power: " << power
                      << ", Term Value: " << termValue
                      << ", Accumulated Result: " << result << std::endl;

        qDebug() << QString::fromStdString(debugMessage1.str());

    }
    //debug Message2
    std::ostringstream debugMessage2;
    debugMessage2 << "Debugging: Final translated function result: " << result << std::endl;
    qDebug() << QString::fromStdString(debugMessage2.str());

    return result;
}
// Method:calculateMagnitude
std::pair<double, std::complex<double>> MagnitudeAndPhase::calculateMagnitude(double angularFrequency) {
    // Calcular el numerador y denominador traducidos
    // Ajuste de la frecuencia compleja s = sigma + jw
    std::complex<double> s = _s_real + std::complex<double>(0, angularFrequency);

    qDebug() << "Debugging Method CalculateMagnitude: ";

    //debug Message
    std::ostringstream debugMessage1;
    debugMessage1 << "Debugging: s = sigma + jw = " << s << std::endl
                  << "Debugging: Magnitude of s = |s| = " << std::abs(s)
                  << std::endl;

    qDebug() << QString::fromStdString(debugMessage1.str());


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
    //debug Message2
    std::ostringstream debugMessage2;
    debugMessage2 << "Numerator: " << resultNumerator << std::endl
                  << "Denominator: " << resultDenominator << std::endl
                  << "Transfer Function |H(jw)| : " << transferFunction << std::endl;
    qDebug() << QString::fromStdString(debugMessage2.str());

    std::cout << "Transfer Function |H(jw)| : " << transferFunction 			       << std::endl
              << "Magnitude (dB): " << magnitudeDB << std::endl;
    return {magnitudeDB, transferFunction};
}

// Method: calculatePhase
double MagnitudeAndPhase::calculatePhase(const std::complex<double>& transferFunction){
    double phaseRadians = std::atan2(transferFunction.imag(), transferFunction.real());

    // Convertir la fase a grados
    double phaseDegrees = phaseRadians * (180.0 / M_PI);

    qDebug()<<"Debugging Method CalculatePhase: ";

    //debug Message
    std::ostringstream debugMessage;
    debugMessage << "Debugging: Transfer Function = " << transferFunction
                 << std::endl
                 << "Debugging: Real = " << transferFunction.real() << std::endl
                 << "Debugging: Imaginary = " << transferFunction.imag()
                 << std::endl;
    qDebug() << QString::fromStdString(debugMessage.str());

    std::cout << "Debugging: Phase (radians) = " << phaseRadians << std::endl
              << "Debugging: Phase (degrees) = " << phaseDegrees << std::endl;

    return phaseDegrees;
}
//Method: isStable
bool MagnitudeAndPhase::isStable(){

    std::vector<std::complex<double>> poles;

    if (poles.empty()) {
        // Si no hay polos, asumimos que el sistema no puede ser estable
        return false;
    }

    for (const auto& pole : poles) {
        if (pole.real() >=0) {
            return false; // unstable system if a real pole part is not negative
        }
    }
    return true; // stable system if all poles have a real negative part
}
