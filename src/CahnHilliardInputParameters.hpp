#ifndef CahnHilliardInputParameters_hpp
#define CahnHilliardInputParameters_hpp
#include <iostream>
#include <iomanip>
/**
 *\file 
 *\class CahnHilliardInputParameters
 *\brief Class for easily handling input parameters of Cahn-Hilliard differential equation solver.
 *
 * This class essentially just holds some values and has an operator to easily output
 * them to a stream.
 */
class CahnHilliardInputParameters 
{
public:

	/// Spatial discretisation step.
    double spaceStep;

    /// Temporal discretisation step.
    double timeStep;

    /// M positive constant from Cahn-Hilliard equation.
    double mConstant;

    /// a constant from the chemical potential.
    double aConstant;

    /// kapa constant from the chemical potential.
    double kConstant;

    /// \phi_0 initial value of order parameter.
    double initialValue;

    /// Maximum magnitude of initial noise.
    double noise;

    /// Number of steps to evolve equation for.
    int totalSteps;

    /// Number of rows in square lattice domain.
    int rowCount;

    /// Number of columns in square lattice domain.
    int colCount;

    /// Name of output directory to save any output into.
    std::string outputName;

    /** 
	 *\brief operator<< overload for outputting the results.
	 *\param out std::ostream reference that is the stream being outputted to.
	 *\param params constant CahnHilliardInputParameters instance to be output.
	 *\return std::ostream reference so the operator can be chained.
	 *
	 * Results will be output in a formatted table for easy viewing in the command line or a file.
	 */
    friend std::ostream& operator<<(std::ostream& out, const CahnHilliardInputParameters& params);

};
#endif /* CahnHilliardInputParameters_hpp */