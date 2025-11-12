#include "Rao_Angles.hpp"
#include "Nozzle_Contour.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

int main() {
    // ---------------------------
    // Define nozzle geometry
    // ---------------------------
    double contraction_ratio = 1.6;      // Area contraction ratio (chamber to throat)
    double expansion_ratio   = 20.0;     // Area expansion ratio (throat to exit)
    double L_chamber         = 0.4613;   // Length of combustion chamber [m]
    double l_perc            = 0.80;     // Percentage of expansion cone length (used in divergent section)
    double r_chamber         = 0.4005;   // Radius of combustion chamber [m]
    double r_throat          = r_chamber / std::sqrt(contraction_ratio); // Radius of the throat [m]
    int N_r = 1;                         // Default resolution or number of radius samples
    
    // ---------------------------
    // Compute Rao optimal nozzle angles
    // ---------------------------
    RaoAngles rao;  
    auto [theta_n, theta_e] = rao.compute(expansion_ratio, l_perc);
    // theta_n: convergent (inlet) half-angle
    // theta_e: divergent  (exit)  half-angle

    // ---------------------------
    // Initialize nozzle contour generator
    // ---------------------------
    NozzleContour nozzle(
        contraction_ratio,       // Chamber to throat contraction ratio
        expansion_ratio,         // Throat to exit expansion ratio
        L_chamber,               // Chamber length
        l_perc,                  // Divergent cone length percentage
        r_chamber,               // Chamber radius

        // Angles (in degrees)
        {20.0},                  // Example angle for convergent section
        {theta_n},               // Rao optimal convergent half-angle
        {theta_e},               // Rao optimal divergent half-angle

        // Radii
        {},                           // r1_values: primary convergent radii, left empty (ignored)
        {1.5 * r_throat, r_throat},   // r2_values: secondary convergent radii
        {0.382 * r_throat, r_throat}, // r3_values: divergent section radii

        // Numerical settings (optional, default N_r = 1)
        // N_r,                   // Could be uncommented to set resolution explicitly

        // Output folder path (optional, default current directory)
        // "C:/Users/Name/Desktop/Nozzle_Contours"  // Path where to save counturs
    );

    // ---------------------------
    // Generate and draw the nozzle contour
    // ---------------------------
    nozzle.DrawContour();

    return 0;
}
