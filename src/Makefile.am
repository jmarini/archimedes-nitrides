bin_PROGRAMS = archimedes
archimedes_SOURCES = \
	materials/AlAs.h \
	materials/AlP.h \
	materials/AlSb.h \
	materials/GaAs.h \
	materials/GaN.h \
	materials/GaP.h \
	materials/GaSb.h \
	materials/Germanium.h \
	materials/InAs.h \
	materials/InP.h \
	materials/InSb.h \
	materials/Silicon.h \
	mep/HMEPbcs.h \
	mep/Hole_bcs.h \
	mep/MEP_interpolation.h \
	mep/ParabMEP2D.h \
	mep/electron_relaxation.h \
	mep/extrema.h \
	mep/hole_relaxation.h \
	mep/holemep2d.h \
	mep/mm.h \
	mep/mm2.h \
	mep/sign.h \
	archimedes.c \
	computecurrents.h \
	configuration.h \
	constants.h \
	deviceconfig.h \
	drift.h \
	electric_field.h \
	ensemblemontecarlo.h \
	faraday.h \
	faradaybcs.h \
	global_defines.h \
	mcparameters.h \
	media.h \
	mesh.c \
	mesh.h \
	particle.h \
	particle.c \
	particlecreation.h \
	particles_per_cell.h \
	poissonbcs.h \
	quantumeffectivepotential.h \
	random.h \
	readinputfile.h \
	saveoutput2dgnuplot.h \
	saveoutput2dholegnuplot.h \
	saveoutput2dholemeshformat.h \
	saveoutput2dmeshformat.h \
	saveoutputfiles.h \
	scattering.h \
	updating.h \
	utility.h \
	vec.h

archimedes_LDADD = -lm
archimedes_CFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-unused-result -pedantic -std=c11 -O3 -fms-extensions
