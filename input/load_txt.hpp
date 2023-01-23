#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "../defs.hpp"
#endif


const int NVARS = 28;
void InitializeArrays(ARRAY2D &COORDS, ARRAY2D &PRIMS, std::string FILE_NAME);
void InitializeArraysMmap(std::string fname);
void InitializeNumpyArrays(ARRAY2D &COORDS, ARRAY2D &PRIMS, std::string fname);