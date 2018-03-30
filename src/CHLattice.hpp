#ifndef CHLattice_hpp
#define CHLattice_hpp

#include <vector> // For holding the values of the function.
#include <algorithm>
#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>

/**
 *\file
 *\class CHLattice
 *\brief Lattice for order parameter which can be evolved in time by the Cahn-Hilliard equation.
 *
 * 2D lattice consisiting of an array of floating points which represent the values of the order parameter
 * at some time t. The lattice can be evolved through time according to the Euler algorithm to give the 
 * evolution of the order parameter.
 */
 class CHLattice
 {
 private:

    /// x-range on lattice.
    int m_xRange;

    /// y-range on lattice.
    int m_yRange;

    /// spatial discretisation step size.
    double m_dx;

    /// M parameter from the Cahn-Hilliard equation.
    double m_M;

    /// ``a'' parameter from the chemical potential.
    double m_a;

    /// Kappa parameter from the chemical potential.
    double m_k;

    /// Vector to hold the values of the order parameter at each lattice site.
    std::vector<double> m_data;

 public:
    /**
     *\brief Creates a lattice of order parameter \phi
     *\param rows integer representing the number of rows in the lattice.
     *\param cols integer representing the number of columns in the lattice.
     *\param m floating point representing the M constant in CH equation.
     *\param a floating point representing the a in the chemical potential.
     *\param k floating point representing the kappa in the chemical potential.
     *\param initialValue floating point value representing the initial value of the order parameter.
     *\param dx floating point value representing temporal discretisation step size.
     */
    CHLattice(int xRange, int yRange, double m, double a, double k, double dx);

    /**
     *\brief Initializes lattice with some value at each site plus some noise
     *\param initialValue initial value at each lattice site.
     *\param maximum magnitude of initial noise which will be uniformly distributed. 
     *\param *\param generator reference to random engine to generate random noise.
     */
    void initialise(double initialValue, double noise, std::default_random_engine &generator);

    /**
     *\brief Calculates the chemical potential for a given point in the lattice.
     *\param i integer value representing x coordinate.
     *\param j integer value representing y coordinate.
     *\return floating point value representing the chemical potential at specified point at current time.
     */
    double chemicalPotential(int i, int j) const;

    /**
     *\brief Calculates the free energy on the lattice.
     *\param i x coordinate of the site at which we wish to calculate chemical potential.
     *\param j y coordinate of the site at which we wish to calculate the chemical potential.
     *\return floating point value representing the free energy at that point.
     */
    double freeEnergy(int i, int j) const;

    /**
     *\brief Prints the free energy to an output stream.
     *\param out ostream reference for data to be streamed to.
     */
    void printFreeEnergy(std::ostream &out) const;

    /**
     *\brief calculates the extensive free energy on the lattice according to the integral over all sites.
     *\return floating point representing the extensive free energy.
     */
    double freeEnergy() const;

    /**
     *\brief Calculates the next value for the order parameter at the (i,j)th lattice site in next step.
     *\param i integer value representing discretised x coordinate.
     *\param j integer value representing discretised y coordinate.
     *\param dt floating point representing discretised time step size.
     *\return floating point value representing the value of the order parameter at the next time step.
     */
    double nextValue(int i, int j, double dt) const;

    /**
     *\brief updates one lattice based on lattice state of other board.
     *\param currentLattice current lattice to update based on.
     *\param updateLattice lattice to be updated.
     *\param dt floating point value representing the value of the order parameter at the next time step.
     */
    friend void update(CHLattice &currentLattice, CHLattice &updateLattice, double dt);

    /**
     *\brief streams the lattice to an output stream in a nicely formatted way
     *\param out std::ostream reference that is being streamed to 
     *\param lattice CHLattice reference to be printed
     *\return std::ostream reference to output can be chained.
     */
     friend std::ostream& operator<<(std::ostream& out, const CHLattice &lattice);
     
     /**
     *\brief operator overload for getting the value of the order parameter at a site.
     *
     * This method is implemented values the spins are stored internally as a 1D vector, hence 
     * they need to be indexed in a special way in order to get the site that would correspond to 
     * the (i,j) site in coordinate notation. This function allows the caller to treat the lattice as a 
     * 2D coordinate system without having to worry about the internal implementation.
     *
     *\param x x index of site.
     *\param y y index of site.
     *\return reference to floating point stored at site so called can use it or set it.
     */
    double& operator()(int x, int y);

    /** 
     *\brief constant version of non-constant counterpart for use with constant ChLattice object.
     *
     * See non-constant version for description.
     *
     *\param x x index of site.
     *\param y y index of site.
     *\return constant reference to floating point stored at site so called can use it only.
     */
    const double& operator()(int x, int y) const;

 };

 #endif /* CHLattice_hpp */