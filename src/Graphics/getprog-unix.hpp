
extern long verbosity;

int getprog(char* fn,int argc, char **argv)
{
    int ret=0;
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
	    else if(ret==0)
	    {
		strcpy(fn,argv[i]);
		ret=1;
	    }
	    if(ret==0) 
	    {
		if(argc>0)
		    cerr << " Syntaxe : " << argv[0] << "  -f filename  [-v verbosity] " << endl;
		else 
		    cerr << " Syntaxe : FreeFem++  -f filename  [-v verbosity] " << endl;

		return ret; 
	    }
	    if(verbosity>10) 
		cout << " file : " << fn << endl ; 
    return 1;
}
