#include <iostream> // For file IO.
#include <boost/filesystem.hpp> // For constructing directories for file IO.
#include <boost/program_options.hpp> // For command line arguments.
#include <fstream> // For file output.
#include <chrono> // For timing.
#include <ctime>  // For timing.
#include <random> // For generating random numbers.
#include <cmath> // For any maths functions.  
#include <iomanip> // For manipulating output.
#include <string> // For naming output directory.
#include <fstream> // For file output.
#include <algorithm> // For swapping the lattices.
#include "Timer.hpp" // For custom timer.
#include "makeDirectory.hpp" // For making directories.
#include "CahnHilliardInputParameters.hpp" // For neatly packaging together input parameters.
#include "CHLattice.hpp"


int main(int argc, char const *argv[])
{
/*************************************************************************************************************************
************************************************* Preparations **********************************************************
*************************************************************************************************************************/
    // Start the clock so execution time can be calculated. 
    Timer timer;

    // Seed the pseudo random number generator using the system clock.
    unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());

    // Create a generator that can be fed to any distribution to produce pseudo random numbers according to that distribution. 
    std::default_random_engine generator(seed);

/*************************************************************************************************************************
******************************************************** Input **********************************************************
*************************************************************************************************************************/
    // Input parameters for the simulation.

    // Spatial discretisation step.
    double spaceStep;

    // Temporal discretisation step.
    double timeStep;

    // M positive constant from Cahn-Hilliard equation.
    double mConstant;

    // a constant from the chemical potential.
    double aConstant;

    // kapa constant from the chemical potential.
    double kConstant;

    // \phi_0 initial value of order parameter.
    double initialValue;

    // Maximum magnitude of initial noise.
    double noise;

    // Number of steps to evolve equation for.
    int totalSteps;

    // Number of x values in square lattice domain.
    int xRange;

    // Number of y values in square lattice domain.
    int yRange;

    // Name of output directory to save any output into.
    std::string outputName;

    // Set up optional command line argument.
    boost::program_options::options_description desc("Options for Ising model simulation");

    desc.add_options()
        ("spatial-discretisation,x", boost::program_options::value<double>(&spaceStep)->default_value(1), "Spatial discretisation step size.")
        ("temporal-discretisation,t", boost::program_options::value<double>(&timeStep)->default_value(1), "Temporal discretisation step size.")
        ("M-Constant,M", boost::program_options::value<double>(&mConstant)->default_value(0.1), "M parameter from Cahn-Hilliard equation.")
        ("a-constant,a", boost::program_options::value<double>(&aConstant)->default_value(0.1),"a parameter from chemical potential.")
        ("k-constant,k", boost::program_options::value<double>(&kConstant)->default_value(0.1),"Kappa parameter from chemical potential.")
        ("initial-value,v", boost::program_options::value<double>(&initialValue)->default_value(0), "Initial value of order parameter.")
        ("noise,p",boost::program_options::value<double>(&noise)->default_value(0.1), "Maximum magnitude of initial noise.")
        ("steps,n", boost::program_options::value<int>(&totalSteps)->default_value(100000),"Total number of steps to evolve differential equation for.")
        ("x-range,r", boost::program_options::value<int>(&xRange)->default_value(100),"Total number of x points in domain of simulation domain.")
        ("y-range,c", boost::program_options::value<int>(&yRange)->default_value(100),"Total number of y points in domain of simulation domain.")
        ("output,o",boost::program_options::value<std::string>(&outputName)->default_value(getTimeStamp()), "Name of output directory to save output files into.")
        ("animate,a","Output the lattice after each update for animation.")
        ("help,h","Display help message.");
  

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc), vm);
    boost::program_options::notify(vm);

    // If the user asks for help display it then exit.
    if(vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    // Construct an input parameter object, this just makes printing a lot cleaner.
    CahnHilliardInputParameters inputParameters
    {
        spaceStep,
        timeStep,
        mConstant,
        aConstant,
        kConstant,
        initialValue,
        noise,
        totalSteps, 
        xRange,
        yRange,
        outputName
    };


/*************************************************************************************************************************
************************************************* Create Output Files ***************************************************
*************************************************************************************************************************/

    // Create an output directory from either the default time stamp or the user defined string.
    makeDirectory(outputName);

    // Create output file for the input parameters.
    std::fstream inputParameterOutput(outputName+"/input.txt",std::ios::out);

    // Create output file for lattice so we can animate the values.
    std::fstream latticeOutput(outputName+"/lattice.dat",std::ios::out);

    // Create an output file for the extensive free energy.
    std::fstream freeEnergy(outputName+"/freeEnergy.dat",std::ios::out);

    // Print input parameters to command line.
    std::cout << inputParameters << '\n';

    // Print input parameters to file. 
    inputParameterOutput << inputParameters << '\n';

/*************************************************************************************************************************
************************************************* The Simulation ********************************************************
*************************************************************************************************************************/

    // Create the lattice to be used in the simulation, we need two one to hold the current state and one to be updated
    // we can then swap them for performance.
    CHLattice currentLattice(xRange,yRange, mConstant, aConstant, kConstant, spaceStep);
    currentLattice.initialise(initialValue, noise, generator);
    CHLattice updatedLattice = currentLattice;

    // Print the initial lattice at t = 0.
    latticeOutput << currentLattice;

    // Variable to hold the current time step
    int t = 0;

    // Print the initial free energy at t = 0.
    freeEnergy << t << ' ' << updatedLattice.freeEnergy() << '\n';
    while(t < totalSteps)
    {
        // Update the lattice based on state at current time.
        update(currentLattice, updatedLattice, timeStep);

        if(vm.count("animate") && 0== t%1000)
        {
            // Order Parameter.
            // Move to the top of the file.
            latticeOutput.seekg(0,std::ios::beg);
            latticeOutput << updatedLattice << std::flush;

        }
        else
        {
            // Calculate and print the extensive free energy.
            freeEnergy << t << ' ' << updatedLattice.freeEnergy() << '\n';
        }   
        
        // Swap the current lattice and updated lattice so no unnecessary copying takes place.
        std::swap(currentLattice, updatedLattice);
        
        // Increment the loop variable.
        t++;

    }




/*************************************************************************************************************************
***********************************************  Output/Clean Up ********************************************************
*************************************************************************************************************************/


    // Report how long the program took to execute.
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Time take to execute(s):    " << std::right << timer.elapsed() << std::endl << std::endl;
    return 0;
}