#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cassert>

typedef int  Tab[1000];
typedef Tab TTab[100];

using namespace std;
int dump(int l,int *t, ostream & cc=cout) 
{
  cc << " { " ;
  for (int i=0;i<l;i++)
    cc << t[i] <<  " " << ( i < l-1 ? ',' : '}' ) ;
}

int dump(int l,int ll,TTab t, ostream & cc=cout) 
{
  cc << " { \n" ;
  for (int i=0;i<l;i++)
    { cc << "\t\t" ;
      dump(ll,t[i],cc);
      cc  <<  " " << ( i < l-1 ? ',' : '}' ) << endl;  ;
    }
}

int f(int k,int *nn, int * aa, int l0,int l1,int l2)
{
  int L[3]={l0,l1,l2};
  
  int i=1;
  for (int j=0;j<k;j++)
    {
    i*=(L[nn[j]]-aa[j]);
    cout << L[nn[j]] << " - " << aa[j]  << " = " << (L[nn[j]]-aa[j]) << "; ";
    }
  return i;
}
int main(int argc,char ** argv)
{
  if(argc<2) return  1;
  
  int k=atoi(argv[1]);
  char * prefix="";
  if(argc>2) prefix=argv[2];
  int i=0;
  Tab  num,num1,cc,ff;
  Tab il,jl,kl;

  //  fonction_i =  $\Pi_{j=0,k-1} (\Lambda_{nn[i][j]} - aa[i][j])  $
  TTab aa,nn;

  ostream * cf = 0;
   ofstream *  ccf=0;
  if(argc>3) cf = ccf = new ofstream(argv[3]);
  string s[1000];
  int e0=2;
  int e1=e0+k-1;
  int e2=e1+k-1;
  int t =e2+k;
 
  ostringstream si;
  for (int ii=0;ii<=k;ii++)
    { 
      int cc=1;
      ostringstream sj;
      for (int jj=0;jj+ii<=k;jj++)
	{
	
	  ostringstream sk;

          for(int kk=0;ii+jj+kk<k;kk++)
	    sk << "*(L2-"<<  kk << ")";
	  int kk = k-ii-jj;

          int l;
	  if      ( ii==k )  l=0;
	  else if ( jj==k )  l=1;
	  else if ( kk==k )  l=2;
	  else if ( ii==0 )  l=e0+kk;
	  else if ( jj==0 )  l=e1+ii;
	  else if ( kk==0 )  l=e2+jj;
          else l=t++;
          num1[l]=100*ii+10*jj+kk;
          num[l]=i++;
	  il[l]=ii;
	  jl[l]=jj;
	  kl[l]=kk;

          s[l]=si.str()+sj.str()+sk.str();
	  sj << "*(L1-"<<  jj << ")";
	  // l=i-1;
	  {
	    int i=0;
	    for (int j=0;j<ii;j++)
	      aa[l][i]=j,nn[l][i++]=0;
	    for (int j=0;j<jj;j++)
	      aa[l][i]=j,nn[l][i++]=1;
	    for (int j=0;j<kk;j++)
	      aa[l][i]=j,nn[l][i++]=2;
	    assert(i==k);
	  }

	}
      si << "*(L0-"<<  ii << ")";
    }
  cout << i << endl;
  for (int l=0;l<i;l++)
    {
      ff[l]=f(k,nn[l],aa[l],il[l],jl[l],kl[l]);
      cout <<  il[l] << " " << jl[l] << " " << kl[l] << " : " ;
      for(int j=0;j<k;j++)
	cout << "( L_" << nn[l][j]<< " - " << aa[l][j] << " ) ";
      cout << "/ " << ff[l]  << endl;
	}
      for (int l=0;l<i;l++)
    cout << setw(3) << num1[l] << " ->  " << num[l] << " " << s[l] << endl;
  for (int l=0;l<i;l++)
    cout << num[l] << ",  ";
  cout << endl;
 
  if(cf ==0) cf = & cout;
  *cf << prefix <<"nn[" << i << "][" << k << "] = " ;
  dump(i,k,nn,*cf);
  *cf << ";\n";
  *cf << prefix <<"aa[" << i << "][" << k << "] = " ;
  dump(i,k,aa,*cf);
  *cf << ";\n";
  *cf << prefix <<"ff[" << i << "] = " ;
  dump(i,ff,*cf);
  *cf << ";\n";
  *cf << prefix <<"il[" << i << "] = " ;
  dump(i,il,*cf);
  *cf << ";\n";
  *cf << prefix <<"jl[" << i << "] = " ;
  dump(i,jl,*cf);
  *cf << ";\n";
  *cf << prefix <<"kl[" << i << "] = " ;
  dump(i,kl,*cf);
  *cf << ";\n";


  if(ccf) delete ccf; // close file  

}
