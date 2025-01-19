# *** Bode Diagram Viewer Project ***

# Requirements for compiling the project

## This project is compiled using qmake. It is recommended to have Qt version 6.7.3 installed along with the Eigen library.
## This project was developed on Windows. If compiling on a different operating system,
## please verify the directory structure on bodeDiagram_Final.pro and ensure the Eigen include paths are correctly configured.

INCLUDEPATH += "path/to/eigen"  # Replace with the correct path to Eigen

---------------------------------------------------------------------------
this project is graphically representing the mathematical calculations generated by my functions, showing magnitude and phase graphs, as well as the zeros and poles, and determining if the system is stable or unstable. The project is divided into two main parts: a visual part, responsible for the user interface, and a mathematical part, which handles the calculations. The files `widget.h` and `widget.cpp` correspond to the visual part, while the mathematical calculations are contained in the files `magnitudeandphase.cpp` and `magnitudeandphase.h`. 

The `main.cpp` file is solely responsible for managing the widget display.

## Visual Part: `widget.h` and `widget.cpp`

In `widget.h` and `widget.cpp`, interactive elements such as input fields (QLine) are managed, where the user can input the numerator, denominator, and the minimum and maximum frequencies, as well as a sigma value. If no sigma is entered, a default value is used. The widget offers the user the flexibility to define the frequency range they wish to observe.

### Key Methods:
- `void Widget::onCalculatePushButtonClicked()`: Handles the calculation and generates magnitude and phase graphs.
- `void Widget::onExportButtonClicked()`: Exports the graphs as a PDF using Qt libraries like QPainter and QPixmap.
- `void Widget::onClearButtonClicked()`: Clears the widget, allowing the user to input new data without restarting the application.

## Mathematical Part: `magnitudeandphase.h` and `magnitudeandphase.cpp`

In the mathematical part, I implemented methods to perform necessary calculations:

- **`processTransferFunction()`**: This method processes the transfer function, calculates the system's zeros and poles from the numerator and denominator coefficients, and returns the results in a human-readable format.

- **`parseCoefficients(const std::string& polynomial)`**: Processes a polynomial string and extracts the coefficients of each term using regular expressions, returning them as a vector.

- **`findRoots(const std::vector<double>& coefficients)`**: Solves the polynomial's roots using the companion matrix method, utilizing Eigen for matrix operations. It adjusts the roots based on a sigma value, if provided.

- **`frequencies()`**: Calculates logarithmic frequencies and corresponding angular frequencies, using a user-defined frequency range (_freqMin and _freqMax).

- **`translateFunction(double angularFrequency, bool isNumerator)`**: Evaluates the transfer function (numerator or denominator) at a specific angular frequency, handling terms using regular expressions and complex number calculations.

- **`calculateMagnitude(double angularFrequency)`**: Calculates the magnitude in dB of the transfer function evaluated at a specific angular frequency.

- **`calculatePhase(const std::complex<double>& transferFunction)`**: Calculates the phase of a complex transfer function using the `atan2` function.

- **`isStable()`**: Checks the stability of the system by examining the poles' position in the complex plane. If any pole has a real part greater than or equal to zero, the system is unstable.

These methods work together to provide the necessary calculations for the graphs and determine the stability of the system.

## Conclusion

This project represents a robust Bode diagram viewer that integrates both the user interface and complex mathematical calculations. It allows for dynamic interaction with the system's parameters, real-time graph updates, and stability analysis. The mathematical methods used provide accurate results, ensuring that the system's behavior is displayed correctly and that its stability is easily assessed.
