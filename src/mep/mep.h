/* mep.h -- This file is part of GNU archimedes

   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

   Copyright (C) 2004-2011 Jean Michel D. Sellier
   <jeanmichel.sellier@gmail.com>
   <jsellier@purdue.edu>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef ARCHIMEDES_MEP_H
#define ARCHIMEDES_MEP_H


#include "global_defines.h"
#include "constants.h"

#include "mep/constants.h"
#include "mep/extrema.h"
#include "mep/sign.h"
#include "mep/mm.h"


double bufx2d[NXM+1][NYM+1];
double bufy2d[NXM+1][NYM+1];
double ux2d[NXM+1][NYM+1][MN3+1];
double uy2d[NXM+1][NYM+1][MN3+1];
double f2d[NXM+1][NYM+1][MN3+1];
double g2d[NXM+1][NYM+1][MN3+1];
double fx2d[NXM+1][NYM+1][MN3+1];
double gy2d[NXM+1][NYM+1][MN3+1];
double c11[7],c12[7],c21[7],c22[7];
double u[7],f[7],g[7],cw[7];

double u2d[NXM+1][NYM+1][MN3+1];    // Hold summary values for electrons per cell, array indexed by mesh node and value type:
                                    //  type = 0: quantum effective potential
                                    //  type = 1: electron density
                                    //  type = 2: running sum of electron x-velocity (divide by MEDIA to get average)
                                    //  type = 3: running sum of electron y-velocity (divide by MEDIA to get average)
                                    //  type = 4: running sum of electron energy     (divide by MEDIA to get average)
double h2d[NXM+1][NYM+1][MN3+1];    // Hold summary values for holes per cell, array indexed by mesh node and value type:
                                    //  type = 0: quantum effective potential
                                    //  type = 1: hole density
                                    //  type = 2: running sum of hole x-velocity (divide by MEDIA to get average)
                                    //  type = 3: running sum of hole y-velocity (divide by MEDIA to get average)
                                    //  type = 4: running sum of hole energy     (divide by MEDIA to get average)


#include "mep/mep_interpolation.h"
#include "mep/electron_bcs.h"
#include "mep/electron_relaxation.h"
#include "mep/electron_mep.h"
#include "mep/hole_bcs.h"
#include "mep/hole_relaxation.h"
#include "mep/hole_mep.h"

extern inline real MM(real a, real b);
extern inline real MM2(real x, real a, real b);
extern inline real sign(real a, real b);
extern inline real minimus(real x, real y);
extern inline real maximus(real x, real y);

#endif
