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

EXTRA_DIST= 	\
A_tutorial.tex                  figures/Cassini.eps             plots/Thwithhole.eps	\
Adapt.tex                       figures/Cycloid1.eps            plots/Thwithouthole.eps	\
AdaptindicatorP2.tex            figures/Cycloid2.eps            plots/Thxy.eps	\
Beam.tex                        figures/Engine.eps              plots/dTh.eps	\
BisonGrammar.tex                figures/FreeFem++-cs.eps        plots/d_Thf.eps	\
figures/L-shape.eps             plots/d_u.eps	\
FF.sty                          figures/L-shape2.eps            plots/eigen11.eps	\
FFF.sty                         figures/LastCOD.eps             plots/eigen12.eps	\
FFGUI_cygwin_huge.tex           figures/LastCOD2.eps            plots/emptymesh-1.eps	\
FluidStruct.tex                 figures/LastPhoto.eps           plots/emptymesh-2.eps	\
FreeBoundary.tex                figures/LastPhoto2.eps          plots/emptymesh.eps	\
LB.tex                          figures/NACA0012.eps            plots/err02.eps	\
LaplaceP1.tex                   figures/P0P1P2P1nc.eps          plots/fastInterpolat.eps	\
LaplaceP1P2h.tex                figures/P1P2.eps                plots/firstTh.eps	\
LaplaceP1bis.tex                figures/SmileFace.eps           plots/firstU.eps	\
LaplaceP1huge.tex               figures/ThreePoint.eps          plots/gnuplot.eps	\
LaplaceRT.tex                   figures/TouchSide.eps           plots/hat.eps	\
ListAM.tex                      figures/U-shape.eps             plots/laRTp.eps	\
MetricFile.tex                  figures/V-shape.eps             plots/lapRTuv.eps	\
NSUawaCahouetChabart.tex        figures/aTutorial.eps           plots/likegnu.eps	\
NewtonRaphson.tex               figures/adaptmesh.eps           plots/logo.eps	\
No-linear-elas.tex              figures/border.eps              plots/lshape.eps	\
NonlinearCG.tex                 figures/cfunc1.eps              plots/lshapeSol.eps	\
Operators.tex                   figures/cfunc2.eps              plots/medit.eps	\
Parallelversionexperimental.tex figures/electro.eps             plots/movemesh.eps	\
Periodic.tex                    figures/electroMesh.eps         plots/nl-elas.eps	\
Readmesh.tex                    figures/mesh_sample.eps         plots/nosplitmesh.eps	\
Region.tex                      figures/naca1.eps               plots/onnewmesh.eps	\
Schwarz-GC.tex                  figures/naca2.eps               plots/onoldmesh.eps	\
Schwarz-no-overlap.tex          figures/projP0.eps              plots/perio4.eps	\
Schwarz-overlap.tex             figures/projP1.eps              plots/period.eps	\
StokesNavierStokes.tex          figures/projP1b.eps             plots/region.eps	\
TypesofLanguage.tex             figures/projP1nc.eps            plots/region_nu.eps	\
addfe-old.tex                   figures/projP2.eps              plots/region_u.eps	\
addfe.tex                       figures/soapfilm.eps            plots/rhoP2.eps	\
bbFiletype.tex                  figures/soapfilm3d.eps          plots/schwarz-no-th.eps	\
cavity.tex                      figures/square.eps              plots/schwarz-no-u.eps	\
dessin.sty                      mainaddfe.tex                   plots/schwarz-no-u0.eps	\
docFFGUI.tex                    manual-full.tex                 plots/schwarz-th.eps	\
docFFGUI_cygwin.tex             manual.tex                      plots/schwarz-u.eps	\
docFFGUI_cygwin_huge.tex        pdfsync.sty                     plots/schwarz-u0.eps	\
docFFGUI_cygwinhuge.tex         plots/FIGII1.eps                plots/secondT.eps	\
docFFGUI_unix.tex               plots/FIGII2.eps                plots/splitmesh.eps	\
docFFGUI_unix_huge.tex          plots/FIGII3.eps                plots/squareb.eps	\
docGUI.tex                      plots/FIGII4.eps                plots/three.eps	\
figures/1stCOD.eps              plots/FIGII5.eps                plots/threeg.eps	\
figures/1stCOD2.eps             plots/FIGII6.eps                plots/trunc0.eps	\
figures/1stPhoto.eps            plots/FIGII7.eps                plots/trunc6.eps	\
figures/1stPhoto2.eps           plots/FIGII8.eps                plots/twosquare.eps	\
figures/BellInit.eps            plots/FIGII9.eps                plots/xyf.eps	\
figures/BellLast.eps            plots/Laplace.eps               styles.sty	\
figures/Bezier.eps              plots/RT0.eps	\
figures/Cardioid.eps            plots/ThrhoP2.eps	


