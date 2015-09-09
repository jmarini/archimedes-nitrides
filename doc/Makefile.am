## Process this file with automake to produce Makefile.in -*-Makefile-*-
#AUTOMAKE_OPTIONS = foreign 1.2h

EXTRA_DIST = archimedes.tex dblfont.sty fdl.tex diode_pictures mesfet_pictures

CLEANFILES = archimedes.aux archimedes.dvi archimedes.log archimedes.maf archimedes.pdf archimedes.ps \
archimedes.toc archimedes.bmt archimedes.mtc* 

DVIPS = @DVIPS@
DVIPDF = @DVIPDF@
LATEX = @LATEX@
#PDFLATEX = @PDFLATEX@

SUFFIXES = .dvi .pdf .tex .ps

# default target
all: @MANUAL@

manual: archimedes.pdf

nomanual:

# To touch it up after changes:
#remake: archimedes.dvi

archimedes.ps: archimedes.dvi

archimedes.dvi: archimedes.tex
	@echo "Building LaTeX tutorial (dvi file)"
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: $(LATEX) archimedes
	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: $(LATEX) archimedes

archimedes.pdf: archimedes.ps
	@echo "Building LaTeX tutorial (PDF version)"
	$(DVIPDF) archimedes.dvi
#	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: $(PDFLATEX) archimedes
#	TEXINPUTS=.:$(srcdir):${TEXINPUTS}: $(PDFLATEX) archimedes

.dvi.ps:
	$(DVIPS) -o $@ $<

dist-hook:
	rm -rf `find $(distdir) -name CVS`

#distclean-local:
#	@if test "$(top_srcdir)" != "$(top_builddir)" ; then \
#	  rm -f eg3.dat; \
#	fi
