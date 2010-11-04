#!/bin/sh -eu
#   the -e option to stop  on error
#   we are with awk to do real operation in the shell
#
bamg=../../bamg-g
NSC2KE=s/NSC2KE

# for awk because in french the number 1/1000 is written 0,001 not 0.001
# to be sure the RADIXCHAR  is '.'   (cf.  Native Language Support)
LANG=C
export LANG

#  some VAR
ifin=2
j=0
INIT=0
LastIteration=0
NBITER=50
HMIN=0.05
HMINGLOBAL=0.0005
HCOEF=0.8

#  clean of the output file 
rm -f [A-Z]*

#  create the geometry file 
awk -f naca.awk </dev/null >MESH_g.msh 

#  create the initial mesh  MESH_0.amdba
$bamg  -g MESH_g.msh -o  MESH_$j.msh -hmax 2 -oamdba MESH_$j.amdba

while [ $j -lt $ifin ] ; do
   #  i = j + 1
   i=`expr $j + 1`

   #  LastIteration = LastIteration + NBITER
   LastIteration=`expr $LastIteration + $NBITER`


   ## set the current MESH
   rm -f  MESH
   ln -s MESH_$j.amdba MESH

   ## create the DATA file for NSC2KE form file data
   ##  change  2 lines for initialisation 
   rm -f DATA
   sed -e "s/^INIT/$INIT/" -e "s/^LastIteration/$LastIteration/" <data >DATA

   echo "--------- NSC2KE iteration $j    -----------"
   $NSC2KE
   ## mv the GNU file 
   for file in GNU.MACH  GNU.MESH  GNU.PRES  GNU.TEMP  GNU.TURB  GNU.VECT;do
     mv -f $file $file.$j 
   done
   ##  find the nb of vertices in the file MESH 
   nbv=`head -1 MESH|awk '{print $1}'`

   ##  create the bb file for interpolation 
   echo "2 4 $nbv  2" > SOL_$j.bb
   cat  SOL_NS  >> SOL_$j.bb

   ## create the bb file for metric construction 
   ## in file SOL_NS  on each line i we have  ro ro*u ro*v  energy 
   ## at vertex i 
   ## + a last line with 2 number last iteration and last time
   echo "2 1 $nbv  2" > MACH.bb
   awk 'NF==4 { print sqrt($2*$2+$3*$3)/$1}' SOL_NS >>  MACH.bb

   ##  put all the residual in one file 
   cat RESIDUAL >>RESIDU

   ## set HMIN = MAX($HMINGLOBAL,$HMIN*$HCOEF) 
   HMIN=`awk "END {c=$HMIN*$HCOEF;c=c<$HMINGLOBAL ?$HMINGLOBAL:c; print c};" </dev/null`

   echo    -b MESH_$j.msh -err 0.001 -errg 0.05 -AbsError \
      -hmin $HMIN -hmax 3 -Mbb MACH.bb  -o  MESH_$i.msh \
      -oamdba MESH_$i.amdba  -raison 2 -rbb SOL_$j.bb -wbb INIT_$i.bb  \
      -maxsubdiv 2 > DATA_bamg
     $bamg  

   

   ##  creation of the INIT_NS for NSC2KE
   ##  remove fisrt line form bb file and  add the last line of  SOL_NS
   sed 1d <INIT_$i.bb >INIT_NS
   tail -1 SOL_NS >>INIT_NS

   # change i and not initialisation
   j=$i
   INIT=1
done
