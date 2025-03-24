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
#include <Eigen/Dense>
#include <QDebug>
#include <QTextEdit>
#include <algorithm>
#include <QString>



// Constructor for default initialization
MagnitudeAndPhase::MagnitudeAndPhase(QWidget* parent)
    : QObject(parent), numeratorLineEdit(nullptr), denominatorLineEdit(nullptr), sigmaLineEdit(nullptr),
    frequencyMaxlineEdit(nullptr), frequencyMinlineEdit(nullptr),
    _s_real(0.0), _freqMin(0), _freqMax(0), _numerator(""), _denominator("")
{
    // Default initialization if needed
}

// Constructor with additional parameters for the transfer function
MagnitudeAndPhase::MagnitudeAndPhase(double freqMin, double freqMax, std::string numerator, std::string denominator ,double s_real, QLineEdit *numeratorLineEdit, QLineEdit *denominatorLineEdit, QLineEdit *sigmaLineEdit, QLineEdit *frequencyMaxlineEdit,
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
    qDebug() << "\033[32mDebugging Method processTransferFunction: \033[0m";

    std::vector<std::complex<double>> poles;
    std::vector<std::complex<double>> zeros;

    auto numCoefficients = parseCoefficients(_numerator); // calls parseCoefficients
    auto denCoefficients = parseCoefficients(_denominator);

    // Calculate zeros and poles
    zeros = findRoots(numCoefficients);  // Parse numerator coefficients
    poles = findRoots(denCoefficients);  // Parse denominator coefficients

    // Debug output for zeros (roots of the numerator) using qDebug
    qDebug() << "Zeros (Roots of the Numerator):";
    for (const auto& zero : zeros) {
        qDebug() << QString("(%1, %2)").arg(zero.real()).arg(zero.imag());
    }

    // Debug output for poles (roots of the denominator) using qDebug
    qDebug() << "Poles (Roots of the Denominator):";
    for (const auto& pole : poles) {
        qDebug() << QString("(%1, %2)").arg(pole.real()).arg(pole.imag());
    }

    // Create result string for further use (returning the result as a string)
    std::ostringstream result;
    result << "Zeros (Roots of the Numerator):";
    for (const auto& zero : zeros) {
        result << zero << " ";
    }
    result << "Poles (Roots of the Denominator):";
    for (const auto& pole : poles) {
        result << pole << " ";
    }

    return result.str();
}



//Method: parseCoefficients (para polos y zeros)
std::vector<double> MagnitudeAndPhase::parseCoefficients(const std::string& polynomial) {
    std::vector<double> coefficients;
    std::string modifiedPoly = polynomial;

    qDebug() << "Debugging Method parseCoefficients: ";
    // Regular expression to match terms of the polynomial
    // viejo std::regex termRegex(R"(([+-]?\s*\d*\.?\d*)\s*(s(\^\d+)?)?)");
     std::regex termRegex ("([+-]?\\d*\\.?\\d+)?([+-]?s(?:\\^\\d+)?(?:\\*\\d+)?)?");
    auto termsBegin = std::sregex_iterator(modifiedPoly.begin(), modifiedPoly.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    int maxPower = 0;
    // Determine the maximum power of 's' from the polynomial
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

    coefficients.resize(maxPower + 1, 0.0); // Resize the coefficients vector based on max power
    // Parse each term and assign coefficients to the appropriate power of 's'
    for (std::sregex_iterator i = termsBegin; i != termsEnd; ++i) {
        std::smatch match = *i;
        if (match.str().empty() || (match[1].str().empty() && match[2].str().empty())) continue;

        std::string coefStr = match[1].str();
        std::string sTerm = match[2].str();


        //  Parse the coefficient (considering signs and empty strings)
        coefStr.erase(remove_if(coefStr.begin(), coefStr.end(), ::isspace), coefStr.end());
        double coefficient = 1.0;

        if (!coefStr.empty() && coefStr != "+" && coefStr != "-") {
            coefficient = std::stod(coefStr);
        } else if (coefStr == "-") {
            coefficient = -1.0;
        }
        else if (sTerm[0] == '-') {
            // Handle cases where the term starts with '-s' but coefStr is empty
            coefficient = -1.0;
        }

        // Extract and apply multiplier if present
        size_t multPos = sTerm.find('*');
        if (multPos != std::string::npos) {
            std::string multiplierStr = sTerm.substr(multPos + 1);
            double multiplier = std::stod(multiplierStr);
            coefficient *= multiplier;
        }
        // Determine the power of 's' based on the term
        int power = 0;
        if (!sTerm.empty()) {
            if (sTerm.find('^') != std::string::npos) {
                power = std::stoi(sTerm.substr(sTerm.find('^') + 1));
            } else {
                power = 1;
            }
        }

        // Store the coefficient at the correct position in the vector
        coefficients[maxPower - power] = coefficient;

        // Mensaje de depuración
        QString debugMessage = QString("Debugging: term: %1, coefficient: %2, power: %3, maxPower: %4")
                                   .arg(QString::fromStdString(match.str()))
                                   .arg(coefficient)
                                   .arg(power)
                                   .arg(maxPower);
        qDebug() << debugMessage;
    }

    // Mensaje de depuración de los coeficientes
    QString debugMessage1 = "Debugging: Coefficients vector: ";
    for (const auto& coef : coefficients) {
        debugMessage1.append(QString::number(coef) + " ");
    }
    qDebug() << debugMessage1;

    return coefficients;
}

//Method: findRoots
std::vector<std::complex<double>> MagnitudeAndPhase::findRoots(const std::vector<double>& coefficients) {
    size_t degree = coefficients.size() - 1;
    qDebug() << "Debugging Method findRoots: ";

    // Validate the degree and coefficients
    if (degree < 1) throw std::invalid_argument("The polynomial degree must be at least 1.");
    if (coefficients[0] == 0) throw std::invalid_argument("The leading coefficient cannot be zero.");
    if (coefficients.size() != degree + 1) throw std::invalid_argument("The size of the coefficients vector must be degree + 1.");

    //debug Message
    QString debugMessage = QString("Debugging: Degree of polynomial: %1\nDebugging: Coefficients: ")
                               .arg(degree);
    for (const auto& coef : coefficients) {
        debugMessage.append(QString::number(coef) + " ");
    }
    qDebug() << debugMessage;


   // Create a companion matrix for solving the polynomial roots
    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();


    //Create a companion matrix for solving the polynomial roots
    for (std::size_t i = 1; i < degree; ++i) {
        companion( i - 1 , i ) = 1.0;
        qDebug() << QString("Subdiagonal: companion(%1, %2) = 1.0").arg(i - 1).arg(i);
    }
    //  Assign the last row of the companion matrix using the coefficients
    for (std::size_t i = 0; i < degree; ++i) {
        companion(degree - 1, i) = -coefficients[degree - i] / coefficients[0];
        qDebug() << QString("Last row: companion(%1, %2) = %3")
                        .arg(degree - 1)
                        .arg(i)
                        .arg(-coefficients[degree - i] / coefficients[0]);
    }

    // Mensaje de depuración para la matriz compañera
    std::stringstream ss;
    ss << companion.format(Eigen::IOFormat(Eigen::StreamPrecision, 0, " ", " "));  // Formatea la matriz y la convierte a string

    debugMessage = QString("Debugging: Companion matrix: ") + QString::fromStdString(ss.str());  // Convierte std::string a QString
    qDebug() << debugMessage;


    // Solve the companion matrix to find the roots (eigenvalues)
    Eigen::EigenSolver<Eigen::MatrixXd> solver(companion);
    Eigen::VectorXcd roots = solver.eigenvalues();

    // Adjust roots based on sigma (real part of 's')
    std::vector<std::complex<double>> result(roots.size());
    for (Eigen::VectorXcd::Index i = 0; i < roots.size(); ++i) {
        result[i] = (std::abs(_s_real) > 1e-9) ? roots[i] + _s_real : roots[i];
        debugMessage = QString("Root before adjustment: %1 + %2i, after adjustment with sigma: %3 + %4i")
                           .arg(roots[i].real())   // Parte real de la raíz original
                           .arg(roots[i].imag())   // Parte imaginaria de la raíz original
                           .arg(result[i].real())  // Parte real de la raíz ajustada
                           .arg(result[i].imag()); // Parte imaginaria de la raíz ajustada

        qDebug() << debugMessage;
    }
    // Mensaje de depuración de las raíces finales (ajustadas o originales)
    qDebug() << "Debugging: Roots (final, adjusted or original): ";
    for (const auto& root : result) {
        // Imprimir la parte real y la parte imaginaria por separado
        qDebug() << QString("Root: %1 + %2i")
                        .arg(root.real())  // Parte real
                        .arg(root.imag()); // Parte imaginaria
    }

    return result;
}
//Method: frequencies
std::pair<std::vector<double>, std::vector<double>> MagnitudeAndPhase::frequencies() {
    std::vector<double> freqVector;
    std::vector<double> angularFreqVector;
    qDebug() << "Debugging Method Frequencies: " ;

    //debug Message
    qDebug() << "Debugging: Frequency range - Min:" << _freqMin << ", Max:" << _freqMax;


    // Calculate log-spaced frequencies
    int numPoints = 100;   // Number of logarithmic frequency points
    for (int i = 0; i < numPoints; i++) {
        double logFreq = std::pow(10, std::log10(_freqMin) + i * (std::log10(_freqMax) - std::log10(_freqMin)) / (numPoints - 1));
        freqVector.push_back(logFreq);
        angularFreqVector.push_back(2 * M_PI * logFreq);

        //debugging Message1
        qDebug() << "Frequency:" << logFreq << ", Angular Frequency:" << angularFreqVector.back() << "rad/s.";
    }

    //debugging Message2
    qDebug() << "Debugging: Frequency vector size:" << freqVector.size();
    qDebug() << "Debugging: Angular Frequency vector size:" << angularFreqVector.size();


    // Returns the pair of vectors
    return std::make_pair(angularFreqVector, freqVector);
}
// Method: translateFunction
std::complex<double> MagnitudeAndPhase::translateFunction(double angularFrequency, bool isNumerator) {
    std::string functionStr = isNumerator ? _numerator : _denominator;
    std::complex<double> jw(0.0, angularFrequency); // j * omega
    std::complex<double> s = (std::abs(_s_real) > 1e-9) ? std::complex<double>(_s_real, angularFrequency) : jw;

    // Expresión regular para evaluar términos de la función
    std::regex termRegex("([+-]?\\d*\\.?\\d+)?([+-]?s(?:\\^\\d+)?(?:\\*\\d+)?)?");
    auto termsBegin = std::sregex_iterator(functionStr.begin(), functionStr.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    std::complex<double> result(0.0, 0.0); // Acumulador del resultado
    qDebug() << "\033[32mDebugging Method translateFunction: \033[0m";

    // Mensaje de depuración inicial
    QString functionStrQt = QString::fromStdString(functionStr);
    qDebug() << "Debugging: Evaluating function:" << functionStrQt;
    qDebug() << "Debugging: Angular frequency (w):" << angularFrequency;

    // Evaluar cada término en la función
    for (auto it = termsBegin; it != termsEnd; ++it) {
        std::smatch match = *it;
        if (match.str().empty()) continue;

        // Depuración del término actual
        QString matchStr = QString::fromStdString(match[0].str());
        qDebug() << "Match[0]:" << matchStr;

        // Extraer coeficientes
        double coefficient = 1.0; // Por defecto
        std::string coefStr = match[1].str();
        if (!coefStr.empty()) {
            coefficient = std::stod(coefStr); // Convertir a número
        } else if (coefStr == "-" || coefStr == "+") {
            coefficient = (coefStr == "-") ? -1.0 : 1.0;
        }

        if (match[0].str()[0] == '-') {
            coefficient = -std::abs(coefficient);
        }

        // Extraer la potencia de "s"
        int power = 0;
        if (!match[2].str().empty()) {
            if (match[2].str().find('^') != std::string::npos) {
                power = std::stoi(match[2].str().substr(match[2].str().find('^') + 1));
            } else {
                power = 1;
            }
        }

        // Evaluar el término y acumular
        std::complex<double> termValue = coefficient * std::pow(s, power);
        result += termValue;

        // Depuración del valor del término y el resultado acumulado
        QString termValueStr = QString("(%1, %2)").arg(termValue.real()).arg(termValue.imag());
        QString resultStr = QString("(%1, %2)").arg(result.real()).arg(result.imag());

        qDebug() << "M Debugging: Term:" << matchStr
                 << ", M Coefficient:" << coefficient
                 << ", M Power:" << power
                 << ", M Term Value:" << termValueStr
                 << ", M Accumulated Result:" << resultStr;
    }

    // Mensaje de depuración final con el resultado acumulado
    QString finalResultStr = QString("Debugging: Final translated function result: (%1, %2)").arg(result.real()).arg(result.imag());
    qDebug() << finalResultStr;

    return result;
}
// Method:calculateMagnitude
// Function to convert a std::complex<double> to QString
QString complexToQString(const std::complex<double>& c) {
    return QString("(%1, %2)").arg(c.real()).arg(c.imag());
}
// Method: calculateMagnitude
std::pair<double, std::complex<double>> MagnitudeAndPhase::calculateMagnitude(double angularFrequency) {
    // Calculate the complex frequency s = sigma + jw
    std::complex<double> s = _s_real + std::complex<double>(0, angularFrequency);

    qDebug() << "[Debug] Method: CalculateMagnitude";

    // Debug Message 1
    qDebug() << "[Debug] s = sigma + jw =" << complexToQString(s);
    qDebug() << "[Debug] Magnitude of s = |s| =" << std::abs(s);

    // Get the results for numerator and denominator functions
    std::complex<double> resultNumerator = translateFunction(angularFrequency, true);
    std::complex<double> resultDenominator = translateFunction(angularFrequency, false);

    // Prevent division by a very small denominator
    if (std::abs(resultDenominator) < 1e-12) {
        std::cerr << "Error: Denominator is too close to zero!" << std::endl;
        return {};
    }

    // Compute transfer function as numerator/denominator
    std::complex<double> transferFunction = resultNumerator / resultDenominator;
    double magnitude = std::abs(transferFunction);
    qDebug() << "[Debug] Magnitude (absolute value):" << magnitude;

    // Compute magnitude in dB
    double magnitudeDB = (magnitude > 0) ? 20 * std::log10(magnitude) : -std::numeric_limits<double>::infinity();
    qDebug() << "[Debug] Magnitude (dB):" << magnitudeDB;

    // Debug Message 2
    qDebug() << "[Debug] Numerator:" << complexToQString(resultNumerator);
    qDebug() << "[Debug] Denominator:" << complexToQString(resultDenominator);
    qDebug() << "[Debug] Transfer Function |H(jw)| :" << complexToQString(transferFunction);
    qDebug() << "[Debug] Magnitude (dB):" << magnitudeDB;

    return {magnitudeDB, transferFunction};
}

// Method: calculatePhase
double MagnitudeAndPhase::calculatePhase(const std::complex<double>& transferFunction){
    // Calculate phase in radians using atan2
    double phaseRadians = std::atan2(transferFunction.imag(), transferFunction.real());

     // Convert phase to degrees
    double phaseDegrees = phaseRadians * (180.0 / M_PI);
    //not sure if i need to adjust  phaseDegrees += 360;   to positive range
    if (phaseDegrees < 0) {

    }

    // Debug messages
    qDebug() << "[Debug] Method: CalculatePhase";
    qDebug() << "[Debug] Phase in degrees =" << phaseDegrees;
    qDebug() << "[Debug] Transfer Function: Real =" << transferFunction.real()
             << ", Imaginary =" << transferFunction.imag();

    qDebug() << "[Debug] Phase (radians) =" << phaseRadians;
    qDebug() << "[Debug] Phase (degrees) =" << phaseDegrees;
    return phaseDegrees;
}
//Method: isStable


bool MagnitudeAndPhase::isStable() {
    qDebug() << "\033[32mDebugging isStable: Calling processTransferFunction()\033[0m";

    // Call processTransferFunction to get the result as a string
    std::string result = processTransferFunction();  // This function returns a std::string
    qDebug() << "\033[34mResult from processTransferFunction():\033[0m" << QString::fromStdString(result);

    // Process the result to separate zeros and poles
    size_t pos = result.find("Poles (Roots of the Denominator):");
    if (pos == std::string::npos) {
        qDebug() << "\033[31mNo poles found. The system is considered unstable.\033[0m";
        return false;  // If the "Poles" string is not found, the system is unstable
    }

    // Extract the poles part
    std::string polesPart = result.substr(pos + std::string("Poles (Roots of the Denominator):").length());

    // Now process the poles with regex
    std::vector<std::complex<double>> poles;
    std::regex polePattern(R"(\((-?\d*\.\d+|\d+),\s*(-?\d*\.\d+|\d+)\))");

    std::sregex_iterator iter(polesPart.begin(), polesPart.end(), polePattern);
    std::sregex_iterator end;

    while (iter != end) {
        double realPart = std::stod((*iter)[1].str());
        double imagPart = (*iter)[2].str().empty() ? 0.0 : std::stod((*iter)[2].str()); // If the imaginary part is empty, set it to 0
        poles.push_back(std::complex<double>(realPart, imagPart));
        qDebug() << "Extracted pole: (" << realPart << ", " << imagPart << ")";  // Debug for each extracted pole
        ++iter;
    }

    // Verify if poles were extracted correctly
    if (poles.empty()) {
        qDebug() << "\033[31mNo poles found. The system is considered unstable.\033[0m";
        return false;  // If no poles are found, the system is unstable
    }

    // Verify stability based on poles
    qDebug() << "\033[34mChecking poles for stability:\033[0m";
    bool stable = true; // Assume stable initially

    for (const auto& pole : poles) {
        qDebug() << QString("Pole: (%1, %2)").arg(pole.real()).arg(pole.imag());
        if (pole.real() > 0) {
            qDebug() << "\033[31mPole has positive real part - contributes to instability.\033[0m";
            stable = false; // Mark as unstable but continue checking other poles
        }
    }

    if (stable) {
        qDebug() << "\033[32mAll poles have negative real part. The system is stable.\033[0m";
    } else {
        qDebug() << "\033[31mSystem is unstable due to one or more poles with positive real part.\033[0m";
    }

    return stable;
}
