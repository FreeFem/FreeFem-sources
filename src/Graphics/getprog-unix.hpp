
extern long verbosity;
extern FILE *ThePlotStream; //  Add for new plot. FH oct 2008
int getprog(char* fn,int argc, char **argv)
{
    
    int ret=0;
    *fn='\0';
    for (int i=1; i<argc;i++)
	if  (ret ==0 && strcmp(argv[i],"-f")==0 && i+1 < argc  ) 
	{
	    strcpy(fn,argv[i+1]);
	    i++;	
	    ret=1;
	}
	    else if  (strcmp(argv[i],"-v")==0 && i+1 < argc) 
	    {
		verbosity = atoi(argv[i+1]);
		i++;	
		if(verbosity>10) printf(" verbosity : %ld\n",verbosity);
	    }
            else if(strcmp(argv[i],"-fglut")==0 && i+1 < argc)
		    {
		      ThePlotStream = fopen(argv[i+1],"w");
		      i++;
		      printf(" save of all plot in file : %s\n",argv[i]);
		      if(!ThePlotStream) { cerr << "  Error fopen  "<< endl;exit(1);}
		    }
	    else if( strcmp(argv[i],"-glut")==0 && i+1 < argc) 
		 {
		    ThePlotStream = popen(argv[i+1],"w");
		   
		    i++;
		    printf(" EXEC of the plot in file : %s\n",argv[i]);
		    if(!ThePlotStream) { cerr << "  Error popen "<< endl;exit(1);}
		  }			    
	    else if(ret==0)
	    {
		strcpy(fn,argv[i]);
		ret=1;
	    }
	    if(ret==0) 
	    {
#ifndef PCRGRAPH_CPP
		if(argc>0)
		    cerr << " Syntaxe : " << argv[0] << "  -f filename  [-v verbosity] " << endl;
		else 
		    cerr << " Syntaxe : FreeFem++  -f filename  [-v verbosity] " << endl;
#endif
		return ret; 
	    }
	    if(verbosity>10) 
		cout << " file : " << fn << endl ; 
    return 1;
}
