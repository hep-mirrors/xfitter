docfiles = writeup.tex
if HAVE_PDFLATEX 
if HAVE_BIBTEX

infiles = $(docfiles)
docfiles += writeup.pdf
MANNAME = writeup
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf

CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANAUX)  $(MANNAME).blg  $(MANNAME).bbl

$(MANPDF): $(srcdir)/$(MANTEXSRC)
	$(PDFLATEX) $<
	$(BIBTEX) writeup
	$(PDFLATEX) $<
	$(PDFLATEX) $<
endif
endif
docdir = ${datadir}/doc/${PACKAGE} 
dist_doc_DATA = $(docfiles)

