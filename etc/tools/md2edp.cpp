#include <fstream>
#include <iostream>
#include <string>
using namespace std;
int  md2edp(const char *fname)
{
   // cout << fname << endl;
    int ok=0;
    string fn = fname;
    string ofn = fn;
    size_t lm=ofn.rfind(".md"),lofn=ofn.size();
    if ( lofn-lm == 3)
      ofn.replace(ofn.end()-3,ofn.end(), ".edp");
    cout << " md2edp " << fn << " -> " << ofn << endl;
    
    ifstream fin(fn.c_str());
    ofstream *fout=0;
    
    //
    long cas = 1;
    if( fin )
    {
        
        while (!fin.eof())
        {
            string ln ;
            getline(fin,ln);
            if(ln.find("~~~freefem")==0) cas++;
            else if ((cas%2 == 0 ) && (ln.find("~~~")==0))  cas++;
            else if( cas%2==0) { 
                if(!fout) { fout = new ofstream(ofn.c_str()); 
                if(!fout) break; 
                *fout << "//  created with md2edp " << fn << endl;
                }
                *fout << ln << endl;}

          //  cout << cas << " ..." << ln << endl;
            
        }
    } else {
        ok=1;cout << " error of file (skip) "<< fname << endl;
    }
    if( !fout) cout << " skip "<< fname << " no freefem code " << endl;
    if(fout) delete fout;
    return ok;
}
int main(int argc,const char ** argv)
{
    int ok =0;
    for(size_t i=1; i< argc;i++)
        if(md2edp(argv[i])) ok =1;
    return ok;
}
