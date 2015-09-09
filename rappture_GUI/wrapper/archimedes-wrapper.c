/* archimedes-wrapper.c -- 
   This code reads the input of the rappture GUI, creates an input file and launch archimedes kernel.
   Archimedes code is a simulator for Submicron 2D III-V
   semiconductor Devices (along with SiO2). It implements the Monte Carlo method 
   and a  simplified version of MEP model for the simulation of the 
   semiclassical Boltzmann equation for both electrons and holes.
   It also includes the quantum effects by means of effective
   potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Jean Michel Sellier
   <jeanmichel.sellier@gmail.com>
 
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


// =================================================================
// File Name : archimedes-wrapper.c
// Version   : release 1.0.0
// Date of Creation : 06 Jan. 2009, West Lafayette, Indiana, USA, Jean Michel Sellier.
// Last Revision : 01 April 2009, West Lafayette, Indiana, USA, Jean Michel Sellier.
// =================================================================

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<memory.h>
#include "rappture.h"
#include<unistd.h>
#include<errno.h>
#include<sys/types.h>
#include<sys/stat.h>

FILE *fp;

int main(int argc,char* argv[])
{
 const char* text0=NULL;
 const char* text1=NULL;
 char buf[256];

// open the rappture library
 RpLibrary* lib=NULL;

 lib=rpLibrary(argv[1]);

 if(lib!=NULL) printf("Rappture Library loaded correctly.\n");
 else{
  printf("Unable to load the Rappture library.\n");
  return(0);
 }

// open the file.input
 fp=fopen("file.input","w");
 if(!fp){
  printf("not able to open the inputfile\n");
//  rpFreeLibrary(lib);
  exit(0);
 }

// #####################
// it reads all the GUI user inputs
// and put them in a file called file.input
 
 fprintf(fp,"# file automatically created by archimedes-gui\n\n");

 register int i;

// #######
// phase 2
// #######
 rpGetString(lib,"input.phase(two).group(time).number(finaltime).current",&text0);
 fprintf(fp,"FINALTIME %s\n\n",text0);
 rpGetString(lib,"input.phase(two).group(time).number(timestep).current",&text0);
 fprintf(fp,"TIMESTEP %s\n\n",text0);

 rpGetString(lib,"input.phase(two).group(lenghts).number(xlenght).current",&text0);
 fprintf(fp,"XLENGTH %s\n\n",text0);
 rpGetString(lib,"input.phase(two).group(lenghts).number(ylenght).current",&text0);
 fprintf(fp,"YLENGTH %s\n\n",text0);

 rpGetString(lib,"input.phase(two).group(lenghts).number(xspatialstep).current",&text0);
 fprintf(fp,"XSPATIALSTEP %s\n\n",text0);
 rpGetString(lib,"input.phase(two).group(lenghts).number(yspatialstep).current",&text0);
 fprintf(fp,"YSPATIALSTEP %s\n\n",text0);

 rpGetString(lib,"input.phase(two).group(temperature).number(latticetemperature).current",&text0);
 fprintf(fp,"LATTICETEMPERATURE ");
 for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
 fprintf(fp,"\n\n"); 

// #######
// phase 1
// #######
 rpGetString(lib,"input.phase(one).group(material).group(materiale1).boolean(materialflag1).current",&text0);
 if(strcmp(text0,"yes")==0){
  fprintf(fp,"MATERIAL X");
  rpGetString(lib,"input.phase(one).group(material).group(materiale1).number(materialxi1).current",&text0);
  fprintf(fp," %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale1).number(materialxf1).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale1).number(materialyi1).current",&text0);
  fprintf(fp,"Y %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale1).number(materialyf1).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale1).choice(material1).current",&text0); 
  for(i=0;i<=strlen(text0);i++){ if(text0[i]!='x') buf[i]=toupper(text0[i]); else buf[i]=text0[i];}
  fprintf(fp,"%s ",buf);
  rpGetString(lib,"input.phase(one).group(material).group(materiale1).number(materialx1).current",&text1);
  if(strstr(text0,"x")!=NULL) fprintf(fp,"%s\n",text1);
  else fprintf(fp,"\n\n");
 }
 rpGetString(lib,"input.phase(one).group(material).group(materiale2).boolean(materialflag2).current",&text0);
 if(strcmp(text0,"yes")==0){
  fprintf(fp,"MATERIAL X");
  rpGetString(lib,"input.phase(one).group(material).group(materiale2).number(materialxi2).current",&text0);
  fprintf(fp," %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale2).number(materialxf2).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale2).number(materialyi2).current",&text0);
  fprintf(fp,"Y %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale2).number(materialyf2).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale2).choice(material2).current",&text0); 
  for(i=0;i<=strlen(text0);i++){ if(text0[i]!='x') buf[i]=toupper(text0[i]); else buf[i]=text0[i];}
  fprintf(fp,"%s ",buf);
  rpGetString(lib,"input.phase(one).group(material).group(materiale2).number(materialx2).current",&text1);
  if(strstr(text0,"x")!=NULL) fprintf(fp,"%s\n",text1);
  else fprintf(fp,"\n\n");
 }
 rpGetString(lib,"input.phase(one).group(material).group(materiale3).boolean(materialflag3).current",&text0);
 if(strcmp(text0,"yes")==0){
  fprintf(fp,"MATERIAL X");
  rpGetString(lib,"input.phase(one).group(material).group(materiale3).number(materialxi3).current",&text0);
  fprintf(fp," %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale3).number(materialxf3).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale3).number(materialyi3).current",&text0);
  fprintf(fp,"Y %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale3).number(materialyf3).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale3).choice(material3).current",&text0); 
  for(i=0;i<=strlen(text0);i++){ if(text0[i]!='x') buf[i]=toupper(text0[i]); else buf[i]=text0[i];}
  fprintf(fp,"%s ",buf);
  rpGetString(lib,"input.phase(one).group(material).group(materiale3).number(materialx3).current",&text1);
  if(strstr(text0,"x")!=NULL) fprintf(fp,"%s\n",text1);
  else fprintf(fp,"\n\n");
 }
 rpGetString(lib,"input.phase(one).group(material).group(materiale4).boolean(materialflag4).current",&text0);
 if(strcmp(text0,"yes")==0){
  fprintf(fp,"MATERIAL X");
  rpGetString(lib,"input.phase(one).group(material).group(materiale4).number(materialxi4).current",&text0);
  fprintf(fp," %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale4).number(materialxf4).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale4).number(materialyi4).current",&text0);
  fprintf(fp,"Y %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale4).number(materialyf4).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale4).choice(material4).current",&text0); 
  for(i=0;i<=strlen(text0);i++){ if(text0[i]!='x') buf[i]=toupper(text0[i]); else buf[i]=text0[i];}
  fprintf(fp,"%s ",buf);
  rpGetString(lib,"input.phase(one).group(material).group(materiale4).number(materialx4).current",&text1);
  if(strstr(text0,"x")!=NULL) fprintf(fp,"%s\n",text1);
  else fprintf(fp,"\n\n");
 }
 rpGetString(lib,"input.phase(one).group(material).group(materiale5).boolean(materialflag5).current",&text0);
 if(strcmp(text0,"yes")==0){
  fprintf(fp,"MATERIAL X");
  rpGetString(lib,"input.phase(one).group(material).group(materiale5).number(materialxi5).current",&text0);
  fprintf(fp," %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale5).number(materialxf5).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale5).number(materialyi5).current",&text0);
  fprintf(fp,"Y %s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale5).number(materialyf5).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(one).group(material).group(materiale5).choice(material5).current",&text0); 
  for(i=0;i<=strlen(text0);i++){ if(text0[i]!='x') buf[i]=toupper(text0[i]); else buf[i]=text0[i];}
  fprintf(fp,"%s ",buf);
  rpGetString(lib,"input.phase(one).group(material).group(materiale5).number(materialx5).current",&text1);
  if(strstr(text0,"x")!=NULL) fprintf(fp,"%s\n",text1);
  else fprintf(fp,"\n\n");
 }

 rpGetString(lib,"input.phase(one).group(transport).choice(transportmodel).current",&text0);
 rpGetString(lib,"input.phase(one).group(transport).choice(transportparticles).current",&text1);
 fprintf(fp,"TRANSPORT %s %s\n\n",text0,text1);

 rpGetString(lib,"input.phase(one).group(leid).boolean(leid).current",&text0);
 if(strcmp(text0,"yes")==0) fprintf(fp,"LEID");

// #######
// phase 3
// #######
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity1).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity2).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity3).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity4).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity5).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity6).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity7).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity8).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity9).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(donors).string(donordensity10).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"DONORDENSITY %s\n\n",text0);

 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity1).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity2).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity3).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity4).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity5).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity6).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity7).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity8).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity9).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);
 rpGetString(lib,"input.phase(three).group(acceptors).string(acceptordensity10).current",&text0);
 if(strlen(text0)!=0) fprintf(fp,"ACCEPTORDENSITY %s\n\n",text0);

// #######
// phase 4
// #######

 rpGetString(lib,"input.phase(four).group(contact1).choice(contactposition1).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact1).string(xicontact1).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact1).string(xfcontact1).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact1).choice(contacttype1).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact1).number(contactvoltage1).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact1).string(contactdensity1).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact1).number(contactvoltage1).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }
 
 rpGetString(lib,"input.phase(four).group(contact2).choice(contactposition2).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact2).string(xicontact2).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact2).string(xfcontact2).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact2).choice(contacttype2).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact2).number(contactvoltage2).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact2).string(contactdensity2).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact2).number(contactvoltage2).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact3).choice(contactposition3).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact3).string(xicontact3).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact3).string(xfcontact3).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact3).choice(contacttype3).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact3).number(contactvoltage3).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact3).string(contactdensity3).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact3).number(contactvoltage3).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact4).choice(contactposition4).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact4).string(xicontact4).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact4).string(xfcontact4).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact4).choice(contacttype4).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact4).number(contactvoltage4).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact4).string(contactdensity4).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact4).number(contactvoltage4).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact5).choice(contactposition5).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact5).string(xicontact5).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact5).string(xfcontact5).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact5).choice(contacttype5).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact5).number(contactvoltage5).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact5).string(contactdensity5).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact5).number(contactvoltage5).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact6).choice(contactposition6).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact6).string(xicontact6).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact6).string(xfcontact6).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact6).choice(contacttype6).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact6).number(contactvoltage6).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact6).string(contactdensity6).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact6).number(contactvoltage6).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact7).choice(contactposition7).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact7).string(xicontact7).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact7).string(xfcontact7).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact7).choice(contacttype7).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact7).number(contactvoltage7).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact7).string(contactdensity7).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact7).number(contactvoltage7).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact8).choice(contactposition8).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact8).string(xicontact8).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact8).string(xfcontact8).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact8).choice(contacttype8).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact8).number(contactvoltage8).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact8).string(contactdensity8).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact8).number(contactvoltage8).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact9).choice(contactposition9).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact9).string(xicontact9).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact9).string(xfcontact9).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact9).choice(contacttype9).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact9).number(contactvoltage9).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact9).string(contactdensity9).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact9).number(contactvoltage9).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

 rpGetString(lib,"input.phase(four).group(contact10).choice(contactposition10).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"CONTACT %s ",text0);
  rpGetString(lib,"input.phase(four).group(contact10).string(xicontact10).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact10).string(xfcontact10).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(four).group(contact10).choice(contacttype10).current",&text0);
  fprintf(fp,"%s ",text0);
  if(strcmp(text0,"OHMIC")==0){
   rpGetString(lib,"input.phase(four).group(contact10).number(contactvoltage10).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
   rpGetString(lib,"input.phase(four).group(contact10).string(contactdensity10).current",&text0);
   fprintf(fp," %s\n\n",text0);
  }
  else{
   rpGetString(lib,"input.phase(four).group(contact10).number(contactvoltage10).current",&text0);
   for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]); 
   fprintf(fp,"\n\n");
  }
 }

// #######
// PHASE 5
// #######
 rpGetString(lib,"input.phase(fives).group(statistics).number(statisticalweight).current",&text0);
 fprintf(fp,"STATISTICALWEIGHT %s\n\n",text0);

 rpGetString(lib,"input.phase(fives).group(statistics).number(media).current",&text0);
 fprintf(fp,"MEDIA %s\n\n",text0);

 rpGetString(lib,"input.phase(fives).group(oxide).choice(oxideposition).current",&text0);
 if(strcmp(text0,"NONE")!=0){
  fprintf(fp,"OXIDE %s ",text0);
  rpGetString(lib,"input.phase(fives).group(oxide).string(oxidexi).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(fives).group(oxide).string(oxidexf).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(fives).group(oxide).string(oxidethickness).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(fives).group(oxide).number(oxidevoltage).current",&text0);
  for(i=0;i<strlen(text0)-1;i++) fprintf(fp,"%c",text0[i]);
 }

 rpGetString(lib,"input.phase(fives).group(Flags).group(uno).boolean(quantum).current",&text0); 
 if(strcmp(text0,"yes")==0) fprintf(fp,"QUANTUMEFFECTS\n\n");
 if(strcmp(text0,"no")==0) fprintf(fp,"NOQUANTUMEFFECTS\n\n");

 rpGetString(lib,"input.phase(fives).group(Flags).group(uno).boolean(faraday).current",&text0); 
 if(strcmp(text0,"yes")==0) fprintf(fp,"FARADAY ON\n\n");
 else fprintf(fp,"FARADAY OFF\n\n");

 rpGetString(lib,"input.phase(fives).group(Flags).group(uno).number(tauw).current",&text0); 
 fprintf(fp,"TAUW %s\n\n",text0);

 rpGetString(lib,"input.phase(fives).group(Flags).group(uno).number(cimp).current",&text0); 
 fprintf(fp,"CIMP %s\n\n",text0);

 rpGetString(lib,"input.phase(fives).group(Flags).group(due).boolean(magneticflag).current",&text0);
 if(strcmp(text0,"yes")==0){
  fprintf(fp,"CONSTANTMAGNETICFIELD ");
  rpGetString(lib,"input.phase(fives).group(Flags).group(due).number(xiconstantmagneticfield).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(fives).group(Flags).group(due).number(yiconstantmagneticfield).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(fives).group(Flags).group(due).number(xfconstantmagneticfield).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(fives).group(Flags).group(due).number(yfconstantmagneticfield).current",&text0);
  fprintf(fp,"%s ",text0);
  rpGetString(lib,"input.phase(fives).group(Flags).group(due).number(valueconstantmagneticfield).current",&text0);
  fprintf(fp,"%s\n\n",text0);
 }

 rpGetString(lib,"input.phase(fives).group(outputformats).boolean(maximini).current",&text0);
 if(strcmp(text0,"yes")==0) fprintf(fp,"MAXIMINI\n\n");
 else fprintf(fp,"NOMAXIMINI\n\n");

// rpGetString(lib,"input.phase(fives).group(outputformats).boolean(saveeachstep).current",&text0);
// if(strcmp(text0,"yes")==0) fprintf(fp,"SAVEEACHSTEP\n\n");

 rpGetString(lib,"input.phase(fives).group(outputformats).choice(outputformat).current",&text0);
 fprintf(fp,"OUTPUTFORMAT %s\n\n",text0);

 fprintf(fp,"# END OF FILE\n");

 fclose(fp);

// ############################
// we run the simulation now...
// ############################
 
 int err;
 err = system("archimedes-kernel-0.8.0 file.input");

// #######################
// we show the results now
// #######################
 rpPutString(lib,"output.cloud(points).about.label","Points",0);
 rpPutString(lib,"output.cloud(points).units","um",0);
 rpPutString(lib,"output.cloud(points).hide","yes",0);
 fp=fopen("./points.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.cloud(points).points",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(density).about.label","Electron Density",0);
 rpPutString(lib,"output.field(density).component.mesh","output.cloud(points)",0);
 fp=fopen("./density.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(density).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(x_velocity).about.label","X-velocity",0);
 rpPutString(lib,"output.field(x_velocity).component.mesh","output.cloud(points)",0);
 fp=fopen("./x_velocity.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(x_velocity).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(y_velocity).about.label","Y-velocity",0);
 rpPutString(lib,"output.field(y_velocity).component.mesh","output.cloud(points)",0);
 fp=fopen("./y_velocity.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(y_velocity).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(potential).about.label","Electrostatic Potential",0);
 rpPutString(lib,"output.field(potential).component.mesh","output.cloud(points)",0);
 fp=fopen("./potential.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(potential).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(quantumpotential).about.label","Electrostatic Potential + Bohm potential",0);
 rpPutString(lib,"output.field(quantumpotential).component.mesh","output.cloud(points)",0);
 fp=fopen("./quantum_potential.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(quantumpotential).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(x_Efield).about.label","X-Electric field",0);
 rpPutString(lib,"output.field(x_Efield).component.mesh","output.cloud(points)",0);
 fp=fopen("./x_Efield.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(x_Efield).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(y_Efield).about.label","Y-Electric field",0);
 rpPutString(lib,"output.field(y_Efield).component.mesh","output.cloud(points)",0);
 fp=fopen("./y_Efield.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(y_Efield).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(energy).about.label","Electron energy",0);
 rpPutString(lib,"output.field(energy).component.mesh","output.cloud(points)",0);
 fp=fopen("./energy.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(energy).component.values",buf,1);
 }
 fclose(fp);

 rpPutString(lib,"output.field(magfield).about.label","Magnetic Field",0);
 rpPutString(lib,"output.field(magfield).component.mesh","output.cloud(points)",0);
 fp=fopen("./magnetic_field.dat","r");
 while(!feof(fp)){
  fgets(buf,128,fp);
  rpPutString(lib,"output.field(magfield).component.values",buf,1);
 }
 fclose(fp);

 err=system("/bin/rm density.dat");
 err=system("/bin/rm x_velocity.dat");
 err=system("/bin/rm y_velocity.dat");
 err=system("/bin/rm potential.dat");
 err=system("/bin/rm quantum_potential.dat");
 err=system("/bin/rm energy.dat");
 err=system("/bin/rm x_Efield.dat");
 err=system("/bin/rm y_Efield.dat");
 err=system("/bin/rm points.dat");
 err=system("/bin/rm magnetic_field.dat");

 printf("All the output files have been saved on your local directory\n");


 rpResult(lib);
// rpFreeLibrary(lib);

 return(0);
}
