docfiles =  herafitter-epjc.tex 
if HAVE_PDFLATEX
if HAVE_BIBTEX

infiles = $(docfiles)
docfiles +=  herafitter-epjc.pdf
MANNAME =  herafitter-epjc
MANTEXSRC = $(MANNAME).tex intro.tex summary.tex method.tex theory.tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf

CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANAUX)  $(MANNAME).blg  $(MANNAME).bbl

$(MANPDF): $(srcdir)/$(MANTEXSRC)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
endif
endif
docdir = ${datadir}/doc/${PACKAGE}
dist_doc_DATA = $(docfiles)
