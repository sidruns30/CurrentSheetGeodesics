/*#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
*/
#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "defs.hpp"
#endif

const int NVARS = 20;

void InitializeArrays(std::string FILE_NAME);
extern const int NVARS;

extern void WriteVectorToFile(std::string fname, std::vector <std::vector <double>> data);
