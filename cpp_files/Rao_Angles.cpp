#include "Rao_Angles.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>

// Helper function for linear interpolation
double RaoAngles::interpolate1D(const std::vector<double>& x, const std::vector<double>& y, double xi) const {
    auto it = std::lower_bound(x.begin(), x.end(), xi);
    if (it == x.begin()) return y.front();
    if (it == x.end()) return y.back();

    size_t i = std::distance(x.begin(), it);
    double x0 = x[i - 1], x1 = x[i];
    double y0 = y[i - 1], y1 = y[i];
    return y0 + (y1 - y0) / (x1 - x0) * (xi - x0);
}

std::pair<double, double> RaoAngles::compute(double ar_input, double l_input) const {

    // Check consistency data Input
    if (l_input < l_vals.front() || l_input > l_vals.back())
        throw std::out_of_range("Value of l must be within [0.6, 0.9]");
    if (ar_input < ar.front() || ar_input > ar.back())
        throw std::out_of_range("Value of ar must be within [4, 100]");

        auto it = std::upper_bound(l_vals.begin(), l_vals.end(), l_input);
        size_t i_high = std::distance(l_vals.begin(), it);
        size_t i_low = i_high - 1;

        // Lower and upper l within l_input
        double l_low = l_vals[i_low];
        double l_up = l_vals[i_high];

        // Get starting and exit angle of l_low and l_up
        double tn1 = interpolate1D(ar, theta_n[i_low], ar_input);
        double tn2 = interpolate1D(ar, theta_n[i_high], ar_input);
        double te1 = interpolate1D(ar, theta_e[i_low], ar_input);
        double te2 = interpolate1D(ar, theta_e[i_high], ar_input);

        // Interpolate starting and exit angle at l_input
        double theta_n_rad = (tn1 + (tn2 - tn1) * (l_input - l_low) / (l_up - l_low));
        double theta_e_rad = (te1 + (te2 - te1) * (l_input - l_low) / (l_up - l_low));

        return {theta_n_rad, theta_e_rad};
}