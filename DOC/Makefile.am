all-local: $(DOCPDF)

doc:$(DOCPS) $(DOCPSGZ) $(DOCPDF)

%.dvi:%.tex
	latex $*
	latex $*
	makeindex $*
	latex $*
%.ps:%.dvi
	dvips -u +psfonts.map -u +pdftex.map -u +Ttbbold.map -z -K manual -o $@
%.pdf:%.tex %.dvi
	pdflatex $*
%.gz:%
	gzip -c  -9 $*  >$*.gz
clean:
	-rm *.dvi *.pdf *~ *.aux *.log *.ps *.toc *.ind *.glo *.out *.blg *ilg *.idx *.bbl *.tmp


EXTRA_DIST=plots/*.eps *.sty *.tex manual.ps.gz manual.pdf
