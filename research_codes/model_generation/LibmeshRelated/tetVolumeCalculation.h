#ifndef TETVOLUMECALCULATION_H
#define TETVOLUMECALCULATION_H

#include "tetgen.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>


double tetVolumeCalculation(tetgenio& out);
double tetVolumeCalculation(std::vector< std::vector<double> >& endo_points, int NoOfEndoNode);

#endif
