#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "../defs.hpp"
#endif
#include "../fluid/getk.hpp"

void GetGeodesics(ARRAY3D &gdscs, const ARRAY2D &X_K);
void WriteGeodesics(const ARRAY3D &gdscs, const std::string fname);