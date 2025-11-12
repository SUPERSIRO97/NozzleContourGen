#ifndef NOZZLE_CONTOUR_HPP
#define NOZZLE_CONTOUR_HPP

#include <vector>
#include <string>

class NozzleContour {
public:
    // Constructor 1
    NozzleContour(double contraction_ratio, double expansion_ratio, double L_chamber,
                  double l_perc, double r_chamber,
                  std::vector<double> theta_values,
                  std::vector<double> theta_n_values,
                  std::vector<double> theta_e_values,
                  std::vector<double> r1_values,
                  std::vector<double> r2_values,
                  std::vector<double> r3_values,
                  int Nr = 1,
                  const std::string& folderPath = ".");

    // Constructor 2
    NozzleContour(double contraction_ratio, double expansion_ratio, double L_chamber,
                  double l_perc, double r_chamber,
                  std::vector<double> theta_values,
                  std::vector<double> theta_n_values,
                  std::vector<double> theta_e_values,
                  std::vector<double> r1_values,
                  std::vector<double> r2_values,
                  std::vector<double> r3_values,
                  const std::string& folderPath);

    // Draw
    void DrawContour();
    void ExportContourToCSV(const std::string& base_filename,
                                       double r1, double r2, double r3,
                                       double theta_deg, double theta_n_deg, double theta_e_deg);
    
private:
    // Equally space grid creation
    std::vector<double> linspace(double start, double end, std::size_t num);
    
    // Maximum radius function
    double MaxRadius(double r_throat, double r_chamber, double theta, double r_input);
    
    // Draw chamber portion
    void DrawChamberSection();
    
    // Draw convergent portion
    std::pair<double, double> DrawConvergentSection(double r1_current, double r2_current, double theta_current);
    
    // Draw divergent portion
    void DrawDivergentSection(double x2, double r3_current, double theta_n_current, double theta_e_current);
    
    // Draw Bezier curve
    void DrawBezierCurve(double x3, double y3, double r3_current, double L_div, double theta_n, double theta_e);

    // Member variables
    double contraction_ratio;      ///< Chamber-to-throat area contraction ratio [-]
    double expansion_ratio;        ///< Throat-to-exit area expansion ratio [-]
    double L_chamber;              ///< Length of combustion chamber [m]
    double l_perc;                 ///< Percentage of divergent cone length to use (0-1) [-]
    double r_chamber;              ///< Radius of combustion chamber [m]

    double r_throat;               ///< Radius of nozzle throat [m]
    double r_ext;                  ///< Radius at nozzle exit   [m]

    std::vector<double> theta_values;   ///< Angles for convergent portion      [deg]
    std::vector<double> theta_n_values; ///< Rao optimal convergent half-angles [deg]
    std::vector<double> theta_e_values; ///< Rao optimal divergent half-angles  [deg]

    std::vector<double> r1_values;      ///< Radii for first circle in convergent section  [m]
    std::vector<double> r2_values;      ///< Radii for second circle in convergent section [m]
    std::vector<double> r3_values;      ///< Radii for divergent section / third circle    [m]

    int Nr;                             ///< Number of samples / resolution for dependent radius [-]

    int contour_counter = 1;            ///< Counter for exported contour files to avoid overwriting
    std::vector<double> x_contour;      ///< X coordinates of generated nozzle contour [m]
    std::vector<double> y_contour;      ///< Y coordinates of generated nozzle contour [m]
    std::string folderPath;             ///< Output folder path for saving CSV files
};

#endif