END {
Pi=3.14159265358979;
i20=20;
i8=8;
r = 5;
c0x = 1;
c0y=0;
print "Dimension",2;
print "MaximalAngleOfCorner 46";


print "Vertices",i20+i20+i8;

# the vertex on naca012  wing ( clock wise)
 for (i=-i20;i<i20;i++) {
   x = i/i20;
   x = x^4; 
   t = 1.008930411365*x;
   y = 5*.12*(0.2969*sqrt(t) - 0.126*t - 0.3516*t^2 + 0.2843*t^3 - 0.1015*t^4);
   if(i<0) y=-y;
   print x,y,3;}
 
# vertex on circle  (counter clock wise)
 for (i=0;i<i8;i++) {
   t=i*Pi*2/i8;
   print c0x+r*(cos(t)),c0y+r*sin(t),5;}


 print "Edges",i20+i20+i8;

#  edge on wing  
 k = 1
 j = i20+i20-1; # previous points
 for (i=0;i<i20+i20;j=i++) 
   { print j+k,i+k,3;} # previous, current vertex

# edge on circle 
 k = i20+i20+1;
 j = i8-1;# previous points
 for (i=0;i<i8;j=i++) 
   { print k+j,k+i,5;} # previous, current vertex

#   one subdomain, region on left side of the wing
#   because clock wise sens.
 print "SubDomain",1;
 print 2,1,1,0;
}
