template<class T>
void  HeapSort(T *c,long n)
{
  long l,j,r,i;
  T crit;
  c--; // on decale de 1 pour que le tableau commence a 1
  if( n <= 1) return;
  l = n/2 + 1;
  r = n;
  while (1) { // label 2
    if(l <= 1 ) { // label 20
      crit = c[r];
      c[r--] = c[1];
    if ( r == 1 ) { c[1]=crit; return;}
    } else  crit = c[--l]; 
    j=l;
    while (1) {// label 4
      i=j;
      j=2*j;
      if  (j>r) {c[i]=crit;break;} // L8 -> G2
      if ((j<r) && (c[j] < c[j+1])) j++; // L5
      if (crit < c[j]) c[i]=c[j]; // L6+1 G4
      else {c[i]=crit;break;} //L8 -> G2
    }
  }
}
