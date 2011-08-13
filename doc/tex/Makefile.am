docfiles = writeup.tex
if HAVE_PDFLATEX

infiles = $(docfiles)
docfiles += writeup.pdf
MANNAME = writeup
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf

CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANAUX)

$(MANPDF): $(srcdir)/$(MANTEXSRC)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
	$(PDFLATEX) $<
endif
dist_doc_DATA = $(docfiles)

