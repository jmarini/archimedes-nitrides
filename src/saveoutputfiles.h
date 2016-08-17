/* saveoutputfiles.h -- This file is part of Archimedes release 0.0.3.
   Archimedes is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements the Monte Carlo method and Hybrid MEP model
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means 
   of effective potential method.

   Copyright (C) 2004, 2005 Jean Michel Sellier <sellier@dmi.unict.it>
 
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.  */


// ######################################################
// Created on 10 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 26 feb.2005, Siracusa, J.M.Sellier
// ######################################################

// Here we save the solutions in output files
// in a format according to the choose done.
// if File_Format = MESHFORMAT then output file in mesh format
// if File_Format = GNUPLOTFORMAT then output file in GNUPLOT format
// the input c is the final index for the output files,
// for example density00c.xyz

void
SaveOutputFiles(int File_Format,int c)
{
 if(File_Format==MESHFORMAT){
  if(g_config->simulation_model==MCE || g_config->simulation_model==MCEH || 
     g_config->simulation_model==MEPE || g_config->simulation_model==MEPEH)
       SaveOutput2D_MeshFormat(c);
  if(g_config->simulation_model==MCH || g_config->simulation_model==MCEH ||
     g_config->simulation_model==MEPH || g_config->simulation_model==MEPEH)
       SaveOutput2DHole_MeshFormat(c);
  return;
 }
 if(File_Format==GNUPLOTFORMAT){
  if(g_config->simulation_model==MCE || g_config->simulation_model==MCEH ||
     g_config->simulation_model==MEPE || g_config->simulation_model==MEPEH)
       SaveOutput2DGNUPLOT(c);
  if(g_config->simulation_model==MCH || g_config->simulation_model==MCEH ||
     g_config->simulation_model==MEPH || g_config->simulation_model==MEPEH) 
       SaveOutput2DHoleGNUPLOT(c);
  return;
 }
 printf("%s: Unknown output file format\n",progname);
 exit(EXIT_FAILURE);
}

// ====================================
