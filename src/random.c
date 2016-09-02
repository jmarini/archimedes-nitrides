#include "random.h"

#include <math.h>

double rnd( ) {
   static double seed = 38467.;
   seed = fmod(1027. * seed, 1048576.);
   return seed / 1048576.;
}
