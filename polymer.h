#include<cmath>
#include<vector>
#include<string>
#include<iostream>

using namespace std;

class polymer
{
	public:
	// default constructor
	polymer(int i = 1000);	
	// methods
	bool	hasInitiated();
	
	int		endGroup(); // return integer
	
 	int		endGroupDyad(); // return two integers
	
	int		endGroupTriad(); // return three integers
	
	int		getDegreeOfPolymerization(); // return integer number of occupied 'sites'
	
	int		repeatUnitAtPosition(int loc);
	
	void	addMonomerToEnd(int type); // increment DP with monomer type
	
	void	removeMonomerFromEnd();
	
	void	printPolymer();
	
	double	compositionOf(int type);
	
	double	numberDyadsOfType(int i,int j);
	
	double	numberTriadsOfType(int i,int j,int k);

	double  numberTetradsOfType(int i, int j, int k, int l);

	double numberPentadsOfType(int i, int j, int k, int l, int m);

	
	vector<int> sequence;
	int	size;
	int	degreeOfPolymerization;
};

