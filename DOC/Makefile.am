documentation:$(DOCPS) $(DOCPSGZ) $(DOCPDF)

%.dvi:%.tex
	-rm $*.dvi $*.aux $*.log \
                $*.toc $*.ind $*.glo $*.out $*.blg $*ilg $*.idx $*.bbl $*.tmp	
	latex $*
	latex $*
	makeindex $*
	latex $*
%.ps:%.dvi
	dvips -u +psfonts.map -u +pdftex.map -u +Ttbbold.map -z -K $* -o $@
%.pdf:%.tex
	-rm $*.pdf $*.aux $*.log \
                $*.toc $*.ind $*.glo $*.out $*.blg $*ilg $*.idx $*.bbl $*.tmp	
	pdflatex $*
	pdflatex $*
	makeindex $*
	pdflatex $*
%.gz:%
	gzip -c  -9 $*  >$*.gz
clean:
	-rm *.dvi *.pdf *~ *.aux *.log *.ps *.ps.gz \
		*.toc *.ind *.glo *.out *.blg *ilg *.idx *.bbl *.tmp

EXTRA_DIST=plots/*.eps figures/*.eps *.sty *.tex manual-full.pdf