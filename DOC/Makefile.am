
documentation: $(DOCPDF) 
docps:  $(DOCPS) $(DOCPSGZ) 
docpdf: $(DOCPDF) 

%.dvi:%.tex cpfigs $(CPFIGS_EPS) 
	-rm $*.dvi $*.aux $*.log \
                $*.toc $*.ind $*.glo $*.out $*.blg $*ilg $*.idx $*.bbl $*.tmp	
	latex $*
	latex $*
	makeindex $*
	latex $*
%.ps:%.dvi
	dvips -u +psfonts.map -u +pdftex.map -u +Ttbbold.map -z -K $* -o $@

%.pdf:%.tex  cpfigs $(CPFIGS_PDF) 
	-rm $*.pdf $*.aux $*.log \
                $*.toc $*.ind $*.glo $*.out $*.blg $*ilg $*.idx $*.bbl $*.tmp	
	pdflatex $*
	pdflatex $*
	makeindex $*
	pdflatex $*
cpfigs: 
	-mkdir cpfigs

%.gz:%
	gzip -c  -9 $*  >$*.gz
clean:
	-rm *.dvi *.pdf *~ *.aux *.log *.ps *.ps.gz \
		*.toc *.ind *.glo *.out *.blg *ilg *.idx *.bbl *.tmp
	-rm -r cpfigs

cpfigs/%.pdf:plots/%.eps
	epstopdf $^ --outfile=$@
cpfigs/%.pdf:figures/%.eps
	epstopdf $^ --outfile=$@
cpfigs/%.eps:plots/%.jpg
	convert $^  $@
cpfigs/%.eps:plots/%.pdf
	pdf2ps $^  $@

## ls */*.eps|column -c 80|sed 's/$/ \\/' 
FIGS_EPS=  \
figures/1stCOD.eps              figures/LastCOD.eps             figures/cfunc1.eps \
figures/1stCOD2.eps             figures/LastCOD2.eps            figures/cfunc2.eps \
figures/1stPhoto.eps            figures/LastPhoto.eps           figures/electro.eps \
figures/1stPhoto2.eps           figures/LastPhoto2.eps          figures/electroMesh.eps \
figures/BellInit.eps            figures/NACA0012.eps            figures/mesh_sample.eps \
figures/BellLast.eps            figures/P0P1P2P1nc.eps          figures/naca1.eps \
figures/Bezier.eps              figures/P1P2.eps                figures/naca2.eps \
figures/Cardioid.eps            figures/SmileFace.eps           figures/projP0.eps \
figures/Cassini.eps             figures/ThreePoint.eps          figures/projP1.eps \
figures/Cycloid1.eps            figures/TouchSide.eps           figures/projP1b.eps \
figures/Cycloid2.eps            figures/U-shape.eps             figures/projP1nc.eps \
figures/Engine.eps              figures/V-shape.eps             figures/projP2.eps \
figures/FreeFem++-cs.eps        figures/aTutorial.eps           figures/soapfilm.eps \
figures/L-shape.eps             figures/adaptmesh.eps           figures/soapfilm3d.eps \
figures/L-shape2.eps            figures/border.eps              figures/square.eps \
 \
plots/BSth.eps                  plots/eigen12.eps               plots/onoldmesh.eps \
plots/BSval.eps                 plots/emptymesh-1.eps           plots/perio4.eps \
plots/FIGII1.eps                plots/emptymesh-2.eps           plots/period.eps \
plots/FIGII2.eps                plots/emptymesh.eps             plots/potential.eps \
plots/FIGII3.eps                plots/err02.eps                 plots/potheat.eps \
plots/FIGII4.eps                plots/fastInterpolat.eps        plots/region.eps \
plots/FIGII5.eps                plots/firstTh.eps               plots/region_nu.eps \
plots/FIGII6.eps                plots/firstU.eps                plots/region_u.eps \
plots/FIGII7.eps                plots/fluidstruct1.eps          plots/rhoP2.eps \
plots/FIGII8.eps                plots/fluidstruct2.eps          plots/rmuonde.eps \
plots/FIGII9.eps                plots/fluidstruct3.eps          plots/schwarz-no-th.eps \
plots/Laplace.eps               plots/gnumembrane.eps           plots/schwarz-no-u.eps \
plots/NScahouetChabart.eps      plots/gnuplot.eps               plots/schwarz-no-u0.eps \
plots/NSprojP.eps               plots/hat.eps                   plots/schwarz-th.eps \
plots/NSprojTh.eps              plots/heatex.eps                plots/schwarz-u.eps \
plots/NSprojU.eps               plots/heatexTh.eps              plots/schwarz-u0.eps \
plots/RT0.eps                   plots/imuonde.eps               plots/secondT.eps \
plots/ThrhoP2.eps               plots/laRTp.eps                 plots/sound.eps \
plots/Thwithhole.eps            plots/lamedeform.eps            plots/sound0.eps \
plots/Thwithouthole.eps         plots/lamevect.eps              plots/splitmesh.eps \
plots/Thxy.eps                  plots/lapRTuv.eps               plots/squareb.eps \
plots/condensor.eps             plots/likegnu.eps               plots/stokes.eps \
plots/condensorth.eps           plots/logo.eps                  plots/tempmuonde.eps \
plots/condersor.eps             plots/lshape.eps                plots/thermic.eps \
plots/convectCG.eps             plots/lshapeSol.eps             plots/thermicvst.eps \
plots/convectDG.eps             plots/medit.eps                 plots/three.eps \
plots/csSnap.eps                plots/membrane.eps              plots/threeg.eps \
plots/dTh.eps                   plots/membraneTh.eps            plots/trunc0.eps \
plots/d_Thf.eps                 plots/movemesh.eps              plots/trunc6.eps \
plots/d_u.eps                   plots/nl-elas.eps               plots/twosquare.eps \
plots/eigen.eps                 plots/nosplitmesh.eps           plots/xyf.eps \
plots/eigen11.eps               plots/onnewmesh.eps \
plots/medit2.eps  		plots/threehsv.eps		plots/hsv.eps \
plots/crimpson.eps              plots/J-bfgs.eps		plots/u-bfgs.eps \
plots/csSnapOld.eps \
plots/logo.eps  plots/arei-Thu.eps      plots/arei-etak.eps \
plots/overlapTh.eps plots/us-ug.eps plots/vs-vg.eps  \
plots/square-0.eps      plots/square-1.eps      plots/square-2.eps 
#  Pour la page de garde FH
FIGS_JPG=plots/fig-alh.jpg       plots/fig-fh.jpg        plots/fig-ko.jpg        plots/fig-op.jpg  \
figures/mi.jpg

FIGS_PDF= plots/LogoCNRS.pdf    plots/LogoLJLL.pdf    plots/LogoUPMC.pdf plots/ffauteurs4.pdf  \
plots/titre-ff.pdf   plots/BG-russe.pdf plots/sanskrit.pdf

EXTRA_DIST= addfe.tex docFFGUI.tex manual-full.tex freefem++doc.tex  manual.tex FF.sty FFF.sty dessin.sty pdfsync.sty styles.sty \
$(FIGS_EPS)    $(FIGS_JPG)    $(FIGS_PDF)  

##  la liste des figures a cree pour la doc (eps on pdf)
CPFIGS_PDF = $(addprefix cpfigs/,$(notdir $(patsubst %.eps,%.pdf,$(FIGS_EPS)) \
))

CPFIGS_EPS = $(addprefix cpfigs/,$(notdir $(patsubst %.pdf,%.eps,$(FIGS_PDF))) \
                                 $(notdir $(patsubst %.jpg,%.eps,$(FIGS_JPG))) \
)


# do not delete the copy figure to long to created 
.PRECIOUS: $(CPFIGS_PDF) $(CPFIGS_EPS)