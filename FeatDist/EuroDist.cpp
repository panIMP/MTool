#include "EuroDist.h"
#include <math.h>

double getEuroDist(const double* vec1, const double* vec2, int dim)
{
	double dist = 0.0;
	double dif = 0.0;
	int d = 0;

	for (d = 0; d < dim; ++d)
	{
		dif = vec1[d] - vec2[d];
		dist += (dif * dif);
	}

	dist = sqrt(dist);

	return dist;
	
}


double getEuroDist2(const double* vec1, const double* vec2, int dim)
{
	double dist = 0.0;
	double dif = 0.0;
	int d = 0;

	for (d = 0; d < dim; ++d)
	{
		dif = vec1[d] - vec2[d];
		dist += (dif * dif);
	}

	return dist;

}