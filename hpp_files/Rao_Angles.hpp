#ifndef RAO_ANGLES_HPP
#define RAO_ANGLES_HPP

#include <vector>
#include <utility>

class RaoAngles {
public:
    // Compute starting and exit Rao optimal angles
    std::pair<double, double> compute(double ar_input, double l_input) const;

private:
    // Expansion ratio
    const std::vector<double> ar = {4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

    // fractional length of conical nozzle
    const std::vector<double> l_vals = {0.6, 0.7, 0.8, 0.9};

    // Starting angle nozzle (in degrees) for different expansion ratios and l_vals
    const std::vector<std::vector<double>> theta_n = {
        {26.5, 28.1, 29.2, 30.15, 30.9, 31.5, 32.0, 34.85, 36.3, 37.25, 37.9, 38.5, 39, 39.5, 39.95, 40.3},         // l = 0.6
        {23.7, 25.1, 26.0, 26.85, 27.45, 28, 28.45, 31.2, 32.7, 33.65, 34.4, 35, 35.5, 35.95, 36.3, 36.6},          // l = 0.7
        {21.6, 23.0, 23.9, 24.65, 25.3, 25.9, 26.3, 28.75, 30.1, 30.95, 31.65, 32.2, 32.65, 33.0, 33.3, 33.55},     // l = 0.8
        {19.9, 20.95, 21.85, 22.65, 23.25, 23.8, 24.25, 27.0, 28.4, 29.4, 30.2, 30.8, 31.35, 31.8, 32.15, 32.45}    // l = 0.9
    };

    // Exit angle nozzle (in degrees) for different expansion ratios and l_vals
    const std::vector<std::vector<double>> theta_e = {
        {20.75, 19.25, 18.3, 17.6, 17, 16.5, 16.15, 14.55, 13.9, 13.45, 13.15, 12.95, 12.755, 12.6, 12.45, 12.35},  // l = 0.6
        {17.2, 15.95, 15.1, 14.5, 14, 13.65, 13.35, 11.9, 11.25, 10.85, 10.5, 10.25, 10.05, 9.9, 9.8, 9.7},         // l = 0.7
        {13.95, 12.85, 12, 11.45, 11, 10.65, 10.4, 9, 8.35, 7.9, 7.55, 7.3, 7.15, 7.05, 7.0, 6.975},                // l = 0.8
        {11.3, 10.2, 9.5, 9, 8.6, 8.3, 8.1, 7, 6.55, 6.3, 6.15, 6.05, 6, 5.95, 5.9, 5.85}                           // l = 0.9
    };

    // Function for 1D linear interpolation
    double interpolate1D(const std::vector<double>& x, const std::vector<double>& y, double xi) const;
};

#endif