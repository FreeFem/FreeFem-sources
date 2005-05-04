documentation:$(DOCPS) $(DOCPSGZ) $(DOCPDF)

%.dvi:%.tex $(FIGS_EPS)
	-rm $*.dvi $*.aux $*.log \
                $*.toc $*.ind $*.glo $*.out $*.blg $*ilg $*.idx $*.bbl $*.tmp	
	latex $*
	latex $*
	makeindex $*
	latex $*
%.ps:%.dvi
	dvips -u +psfonts.map -u +pdftex.map -u +Ttbbold.map -z -K $* -o $@
%.pdf:%.tex $(patsubst %.eps,%..pdf,$(FIGS_EPS))
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
%.pdf:%.eps
	epstopdf $*.eps --outfile=$*.pdf
## ls */*.eps|column -c 80|sed 's/$/ \\/' 
FIGS_EPS=  \
figures/1stCOD.eps              plots/Laplace.eps \
figures/1stCOD2.eps             plots/RT0.eps \
figures/1stPhoto.eps            plots/ThrhoP2.eps \
figures/1stPhoto2.eps           plots/Thwithhole.eps \
figures/BellInit.eps            plots/Thwithouthole.eps \
figures/BellLast.eps            plots/Thxy.eps \
figures/Bezier.eps              plots/dTh.eps \
figures/Cardioid.eps            plots/d_Thf.eps \
figures/Cassini.eps             plots/d_u.eps \
figures/Cycloid1.eps            plots/eigen11.eps \
figures/Cycloid2.eps            plots/eigen12.eps \
figures/Engine.eps              plots/emptymesh-1.eps \
figures/FreeFem++-cs.eps        plots/emptymesh-2.eps \
figures/L-shape.eps             plots/emptymesh.eps \
figures/L-shape2.eps            plots/err02.eps \
figures/LastCOD.eps             plots/fastInterpolat.eps \
figures/LastCOD2.eps            plots/firstTh.eps \
figures/LastPhoto.eps           plots/firstU.eps \
figures/LastPhoto2.eps          plots/gnuplot.eps \
figures/NACA0012.eps            plots/hat.eps \
figures/P0P1P2P1nc.eps          plots/laRTp.eps \
figures/P1P2.eps                plots/lapRTuv.eps \
figures/SmileFace.eps           plots/likegnu.eps \
figures/ThreePoint.eps          plots/logo.eps \
figures/TouchSide.eps           plots/lshape.eps \
figures/U-shape.eps             plots/lshapeSol.eps \
figures/V-shape.eps             plots/medit.eps \
figures/aTutorial.eps           plots/movemesh.eps \
figures/adaptmesh.eps           plots/nl-elas.eps \
figures/border.eps              plots/nosplitmesh.eps \
figures/cfunc1.eps              plots/onnewmesh.eps \
figures/cfunc2.eps              plots/onoldmesh.eps \
figures/electro.eps             plots/perio4.eps \
figures/electroMesh.eps         plots/period.eps \
figures/mesh_sample.eps         plots/region.eps \
figures/naca1.eps               plots/region_nu.eps \
figures/naca2.eps               plots/region_u.eps \
figures/projP0.eps              plots/rhoP2.eps \
figures/projP1.eps              plots/schwarz-no-th.eps \
figures/projP1b.eps             plots/schwarz-no-u.eps \
figures/projP1nc.eps            plots/schwarz-no-u0.eps \
figures/projP2.eps              plots/schwarz-th.eps \
figures/soapfilm.eps            plots/schwarz-u.eps \
figures/soapfilm3d.eps          plots/schwarz-u0.eps \
figures/square.eps              plots/secondT.eps \
plots/FIGII1.eps                plots/splitmesh.eps \
plots/FIGII2.eps                plots/squareb.eps \
plots/FIGII3.eps                plots/three.eps \
plots/FIGII4.eps                plots/threeg.eps \
plots/FIGII5.eps                plots/trunc0.eps \
plots/FIGII6.eps                plots/trunc6.eps \
plots/FIGII7.eps                plots/twosquare.eps \
plots/FIGII8.eps                plots/xyf.eps \
plots/FIGII9.eps 


EXTRA_DIST= addfe.tex docFFGUI.tex manual-full.tex manual.tex FF.sty FFF.sty dessin.sty pdfsync.sty styles.sty \
$(FIGS_EPS)