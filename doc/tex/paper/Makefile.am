docfiles =  herafitter-epjc.tex 
if BUILD_DOCS

infiles = $(docfiles)
docfiles +=  herafitter-epjc.pdf
MANNAME =  herafitter-epjc
MANTEXSRC = $(MANNAME).tex intro.tex summary.tex method.tex theory.tex examples.tex alternative.tex techniques.tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf
#EXTRA_DIST = $(MANTEXSRC)

CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANAUX)  $(MANNAME).blg  $(MANNAME).bbl 

$(MANPDF): $(srcdir)/$(MANTEXSRC)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
endif
docdir = ${datadir}/doc/${PACKAGE}
# dist_doc_DATA = $(docfiles)
