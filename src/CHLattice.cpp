#include "CHLattice.hpp"


CHLattice::CHLattice(int xRange, int yRange, double m, double a, double k, double dx): m_xRange(xRange),
																			m_yRange(yRange),
																			m_M(m),
																			m_a(a),
																			m_k(k),
																			m_dx(dx),
																			m_data(xRange*yRange,0.0)
{



}

void CHLattice::initialise(double initialValue, std::default_random_engine &generator)
{
	static std::uniform_real_distribution<double> distribution(-0.01,0.01);

	for(auto &phi : m_data) 
	{
		phi = initialValue + distribution(generator);
	}
}

// The formula for chemical potential is calculated according to (25) in notes.
double CHLattice::chemicalPotential(int i, int j) const
{
	return (- m_a * (*this)(i,j) + m_a * pow((*this)(i,j), 3) 
			- m_k/(pow(m_dx,2)) * ((*this)(i+1,j)+(*this)(i-1,j) 
				+ (*this)(i,j+1)+(*this)(i,j-1)- 4 * (*this)(i,j)));


}

// Formula for free energy is computed according to (4) in notes.
double CHLattice::freeEnergy(int i, int j) const
{
	double gradSquaredTerm = pow(((*this)(i+1,j)-(*this)(i-1,j))/(2*m_dx),2) 
							+ pow(((*this)(i,j+1)-(*this)(i,j-1))/(2*m_dx),2);

	return (-m_a/2 * pow((*this)(i,j),2) + m_a/4 * pow((*this)(i,j),4) + m_k/2 * gradSquaredTerm);
}

void CHLattice::printFreeEnergy(std::ostream &out) const
{
	for(int j = m_yRange - 1; j >=0; --j)
 	{
 		for(int i = 0; i < m_xRange; ++i)
 		{
 			out << std::showpos << std::fixed << std::setprecision(6);
 			out << freeEnergy(i, j) << ' ';
 		}
 		out << '\n';
 	}

 	out << std::flush;

}

double CHLattice::nextValue(int i, int j, double dt) const
{
	return ((*this)(i,j) + m_M*dt/(pow(m_dx,2)) * (chemicalPotential(i+1, j)+chemicalPotential(i-1, j)
		+ chemicalPotential(i, j+1)+ chemicalPotential(i, j-1) - 4 * chemicalPotential(i, j)));

}


void update(CHLattice &currentLattice, CHLattice &updateLattice, double dt)
{
	for(int i = 0; i < currentLattice.m_xRange; ++i)
	{
		for(int j = 0; j < currentLattice.m_yRange; ++j)
		{
			updateLattice(i,j) = currentLattice.nextValue(i, j, dt);
		}
	}

}


 std::ostream& operator<<(std::ostream& out, const CHLattice &lattice)
 {
 	for(int j = lattice.m_yRange - 1; j >=0; --j)
 	{
 		for(int i = 0; i < lattice.m_xRange; ++i)
 		{
 			out << std::showpos << std::fixed << std::setprecision(6);
 			out << lattice(i,j) << ' ';
 		}
 		out << '\n';
 	}

 	return out;

 }
 
double& CHLattice::operator()(int x, int y)
{
	// Take into account periodic boundary conditions we add extra m_xRange and m_yRange
	// terms here to take into account the fact that the caller may be indexing with -1.
	x = (x + m_xRange) % m_xRange;
	y = (y + m_yRange) % m_yRange;

	return m_data[x + y * m_xRange];

}

const double& CHLattice::operator()(int x, int y) const
{
	// Take into account periodic boundary conditions we add extra m_xRange and m_yRange
	// terms here to take into account the fact that the caller may be indexing with -1.
	x = (x + m_xRange) % m_xRange;
	y = (y + m_yRange) % m_yRange;

	return m_data[x + y * m_xRange];

}

 
