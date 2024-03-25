#include "polymer.h"

// Implementation
/****************************************************************************************/
polymer::polymer(int length)
{
	size = length;
	degreeOfPolymerization = 0;
	sequence.resize(size);
	for (int i=0;i<=sequence.size()-1;i++)
	{
		sequence[i] = 0;
	}
}
/****************************************************************************************/
bool polymer::hasInitiated()
{
	if (degreeOfPolymerization == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}
/****************************************************************************************/
int	polymer::endGroup()
{
	if (degreeOfPolymerization >= 1)
	{
		return sequence[degreeOfPolymerization - 1];
	}
	else
	{
		return 123456 //just not 1 or 2
	}
}
/****************************************************************************************/
int	polymer::endGroupDyad()
{
	int value = 0;
	if (degreeOfPolymerization >= 2)
	{
		value=10*sequence[degreeOfPolymerization-2]
			 + sequence[degreeOfPolymerization-1];
	}
	return value;
}
/****************************************************************************************/
int	polymer::endGroupTriad()
{
	int value = 0;
	if (degreeOfPolymerization >= 3)
	{
		value=100*sequence[degreeOfPolymerization-3] 
				+ 10*sequence[degreeOfPolymerization-2]
				+ sequence[degreeOfPolymerization-1];
	}
	return value;
}
/****************************************************************************************/
int	polymer::getDegreeOfPolymerization()
{
	return degreeOfPolymerization;
}
/****************************************************************************************/
int	polymer::repeatUnitAtPosition(int loc)
{
	if (loc < degreeOfPolymerization)
	{
		return sequence[loc];
	}
	else
	{
		return 0;
	}
}
/****************************************************************************************/
void polymer::addMonomerToEnd(int type) // this needs work
{
	if (degreeOfPolymerization < sequence.size())
	{
		sequence[degreeOfPolymerization] = type;
	}
	else // add more size
	{
		sequence.push_back(100);
		// then put 
		sequence[degreeOfPolymerization] = type;
	}
	degreeOfPolymerization += 1;
	return;
}
/****************************************************************************************/
void polymer::removeMonomerFromEnd()
{
	sequence[degreeOfPolymerization-1] = 0;
	degreeOfPolymerization -= 1;
	return;
}
/****************************************************************************************/
void polymer::printPolymer()
{
	for (int i=0;i<=degreeOfPolymerization-1;i++)
	{
		cout << "(" << i <<")=" << sequence[i] << ", ";
	}
	cout << endl;
	return;
}
/****************************************************************************************/
double polymer::compositionOf(int type)
{
	double sum = 0.;
	for (int i=0;i<=sequence.size()-1;i++)
	{
		sum += double(sequence[i]);
	}
	sum /= double(degreeOfPolymerization);
	if (type == 1)
	{
		return (2.0-sum);
	}
	else if (type ==2)
	{
		return (sum - 1.0);
	}
	else
	{
		cout << "polymer::compositionOf(int type) - type unrecognized" << endl;
		return 0.0;
	}
}
/****************************************************************************************/
double polymer::numberDyadsOfType(int i, int j)
{
	double	dyads		= 0.;
	double	total_dyads	= 0.;
	if (degreeOfPolymerization >= 2 )
	{
		for (int k=0;k<=degreeOfPolymerization-2;k++)
		{
			if (sequence[k]==i && sequence[k+1]==j)
			{
				dyads		+= 1.0;
			}
			else if (k < degreeOfPolymerization)
			{
				total_dyads	+= 1.0;
			}
		}
		return dyads;
	}
	else
	{
		return 0.;
	}
}
/****************************************************************************************/
double polymer::numberTriadsOfType(int i, int j, int k)
{
	double	triads			= 0.;
	double	total_triads	= 0.;
	if (degreeOfPolymerization >= 3)
	{
		for (int l=0;l<=degreeOfPolymerization-3;l++)
		{
			if (sequence[l]==i && sequence[l+1]==j && sequence[l+2]==k)
			{
				triads			+= 1.0;
				total_triads	+= 1.0;
			}
			else if (k < degreeOfPolymerization)
			{
				total_triads	+= 1.0;
			}
		}
		return triads;
	}
	else
	{
		return 0.;
	}
}
/****************************************************************************************/
double polymer::numberTetradsOfType(int i, int j, int k, int l)
{
	double tetrads				= 0.;
	double total_tetrads		= 0.;
	if (degreeOfPolymerization >= 4)
	{
		for (int m = 0; m <= degreeOfPolymerization - 4; m++)
		{
			if (sequence[m] == i && sequence[m + 1] == j && sequence[m + 2] == k && sequence[m + 3] == l)
			{
				tetrads			+= 1.0;
				total_tetrads	+= 1.0;
			}
			else if (m < degreeOfPolymerization)
			{
				total_tetrads	+= 1.0;
			}
		}
		return tetrads;
	}
	else
	{
		return 0.;
	}
}

/****************************************************************************************/
double polymer::numberPentadsOfType(int i, int j, int k, int l, int m)
{
	double pentads			=0.;
	double total_pentads	=0.;
	if (degreeOfPolymerization >= 5)
	{
		for(int n=0;n<=degreeOfPolymerization-5;n++)
			if (sequence[n] == i && sequence[n + 1] == j && sequence[n + 2] == k && sequence[n + 3] == l && sequence[n + 4] == m)
			{
				pentads			+=1.0;
				total_pentads	+=1.0;
			}
			else if (n < degreeOfPolymerization)
			{
				total_pentads +=1.0;
			}
		return pentads;
	}
	else
	{
		return 0.;
	}
}

/****************************************************************************************/
