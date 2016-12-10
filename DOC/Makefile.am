# ======================================================================
# Written by Antoine Le Hyaric
# http://www.ljll.math.upmc.fr/lehyaric
# ======================================================================
# headeralh default=0 make start=07/05/2014 written

# ALH - 7/5/14 - <<DOC_not_in_FFCS>> The DOC subdirectory is big (100Mb) and not useful for FFCS. So we just leave it
# out during FFCS compilation. Tip found at [[http://stackoverflow.com/questions/15874325/conditional-ac-output]]

## ls */*.eps|column -c 80|sed 's/$/ \\/' 
all-local:documentation
dist_pkgdata_DATA=$(DOCPDF)
documentation: $(DOCPDF) 
docps:  $(DOCPS) $(DOCPSGZ) 
docpdf: $(DOCPDF) 

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
figures/mach2r.eps \
figures/buillayermesh.eps figures/layer2D.eps figures/degenerate.eps \
figures/cube-bal-perio-medit.eps	figures/cube-bal-perio.eps \
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
plots/square-0.eps      plots/square-1.eps      plots/square-2.eps  \
plots/cube.eps plots/cone.eps \
plots/Hex-Sphere.eps plots/Cube-With-Ball.eps \
plots/multiendborder.eps plots/multiendmesh.eps \
plots/leman.eps 

freefem++doc.pdf:
	$(MAKE) -f gnu.mak $@ EPSTOPDF="$(EPSTOPDF)"
freefem++doc.ps:
	$(MAKE) -f gnu.mak $@ EPSTOPDF="$(EPSTOPDF)"
freefem++doc.ps.gz:
	$(MAKE) -f gnu.mak $@ EPSTOPDF="$(EPSTOPDF)"
clean:
	$(MAKE) -f gnu.mak $@

#  Pour la page de garde FH
FIGS_JPG=plots/fig-alh.jpg       plots/fig-fh.jpg        plots/fig-ko.jpg        plots/fig-op.jpg  \
figures/mi.jpg plots/LogoUPMC.jpg \
plots/beam-3d.jpg \
plots/logo-finance-par-anr.jpg \
plots/VarIneqFill.jpg plots/VarIneqIso.jpg \
plots/NSNewtonTh.jpg  plots/NSNewtonUP.jpg \
plots/chesapeake-2.jpg \
plots/minsurf3D.jpg \
plots/lg.jpg plots/INRIA-logo.jpg  \
figures/func-cube.jpg 

FIGS_PDF= plots/LogoCNRS.pdf    plots/LogoLJLL.pdf  plots/ffauteur.pdf  \
plots/titre-ff.pdf   plots/BG-russe.pdf plots/sanskrit.pdf \
plots/nlopttab.pdf

EXTRA_DIST= addfe.tex docFFGUI.tex manual-full.tex freefem++doc.tex  manual.tex FF.sty FFF.sty dessin.sty pdfsync.sty styles.sty petschpddm.tex 	macroListing.sty \
PARALLEL/FreefemProgram/diffusion-3d.edp \
PARALLEL/FreefemProgram/elasticity-3d.edp \
PARALLEL/FreefemProgram/Schwarz/additional_macro.idp \
PARALLEL/FreefemProgram/Schwarz/macro.idp \
PARALLEL/part.png \
$(FIGS_EPS)    $(FIGS_JPG)    $(FIGS_PDF)  

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
