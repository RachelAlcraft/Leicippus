// Leucippus.cpp : Defines the entry point for the application.
//

#include "Leucippus.h"
#include "Lattice.h"

using namespace std;

int main()
{
	// Diffraction Q 3
	Lattice* lat1 = new RealLattice(100.02, 90.57, 68.33, 90, 104.48, 90);
	Lattice* rec1 = lat1->makeInverseLattice();	
	string infoA1 = lat1->printLattice();	
	string infoB1 = rec1->printLattice();
	cout << "Answers to Question 3\n";
	cout << infoA1 << endl;
	cout << infoB1 << endl;

	// Diffraction Q 4
	Lattice* lat2 = new RealLattice(45.69, 150.22, 75.87, 90, 90, 90);
	Lattice* rec2 = lat2->makeInverseLattice();
	string infoA2 = lat2->printLattice();
	string infoB2 = rec2->printLattice();
	cout << "Answers to Question 4\n";
	cout << infoA2 << endl;
	cout << infoB2 << endl;

	//Add bragg's law in real and reciprocal space

		
	return 0;
}
