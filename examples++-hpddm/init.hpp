#ifdef DHYPRE
    const char ds = 'G';
#else
    const char ds = 'S';
#endif
    const char zs = 'G';
#ifdef SCHWARZ
    Schwarz::add<HpSchwarz, double, ds>();
    zzzfff->Add("dschwarz", atype<HpSchwarz<double, ds>*>());
#ifndef DHYPRE
    // Schwarz::add<HpSchwarz, float, ds>();
    // zzzfff->Add("sschwarz", atype<HpSchwarz<float, ds>*>());
    Schwarz::add<HpSchwarz, std::complex<double>, zs>();
    zzzfff->Add("zschwarz", atype<HpSchwarz<std::complex<double>, zs>*>());
    // Schwarz::add<HpSchwarz, std::complex<float>, zs>();
    // zzzfff->Add("cschwarz", atype<HpSchwarz<std::complex<float>, zs>*>());
#endif
#ifdef WITH_PETSC
    int argc = pkarg->n;
    char** argv = new char*[argc];
    for(int i = 0; i < argc; ++i)
        argv[i] = const_cast<char*>((*(*pkarg)[i].getap())->c_str());
    PetscInitialize(&argc, &argv, 0, "");
    delete [] argv;
    ff_atend(finalizePETSc);
    Dcl_Type<DistributedCSR*>(Initialize<DistributedCSR>, Delete<DistributedCSR>);
    Dcl_Type<GMV<DistributedCSR*, KN<double>*> >();
    zzzfff->Add("dmatrix", atype<DistributedCSR*>());
    TheOperators->Add("<-", new OneOperator1_<long, DistributedCSR*>(initEmptyCSR));
    TheOperators->Add("<-", new initCSR<double>);
    TheOperators->Add("*", new OneOperator2<GMV<DistributedCSR*, KN<double>*>, DistributedCSR*, KN<double>*>(Build));
    TheOperators->Add("=", new OneOperator2<KN<double>*, KN<double>*, GMV<DistributedCSR*, KN<double>*> >(GlobalMV));
    Global.Add("set", "(", new setOptions<double>());
    Dcl_Type<DistributedCSR_inv>();
    TheOperators->Add("^", new OneBinaryOperatorPETSc());
    Dcl_Type<Inv<DistributedCSR_inv, KN<double>*> >();
    TheOperators->Add("*", new OneOperator2<Inv<DistributedCSR_inv, KN<double>*>, DistributedCSR_inv, KN<double>*>(Build));
    TheOperators->Add("=", new OneOperator2<KN<double>*, KN<double>*, Inv<DistributedCSR_inv, KN<double>*> >(InvPETSc));
#endif
#endif
#if defined(BDD)
    Substructuring::add<HpBdd, double, ds>();
    zzzfff->Add("dbdd", atype<HpBdd<double, ds>*>());
#ifndef DHYPRE
    // Substructuring::add<HpBdd, float, ds>();
    // zzzfff->Add("sbdd", atype<HpBdd<float, ds>*>());
    Substructuring::add<HpBdd, std::complex<double>, zs>();
    zzzfff->Add("zbdd", atype<HpBdd<std::complex<double>, zs>*>());
    // Substructuring::add<HpBdd, std::complex<float>, zs>();
    // zzzfff->Add("cbdd", atype<HpBdd<std::complex<float>, zs>*>());
#endif
#endif
#if defined(FETI)
    Substructuring::add<HpFetiPrec, double, ds>();
    zzzfff->Add("dfeti", atype<HpFetiPrec<double, ds>*>());
#ifndef DHYPRE
    // Substructuring::add<HpFetiPrec, float, ds>();
    // zzzfff->Add("sfeti", atype<HpFetiPrec<float, ds>*>());
    Substructuring::add<HpFetiPrec, std::complex<double>, zs>();
    zzzfff->Add("zfeti", atype<HpFetiPrec<std::complex<double>, zs>*>());
    // Substructuring::add<HpFetiPrec, std::complex<float>, zs>();
    // zzzfff->Add("cfeti", atype<HpFetiPrec<std::complex<float>, zs>*>());
#endif
#endif
    // Dcl_Type<Pair<float>*>(InitP<Pair<float> >, Destroy<Pair<float> >);
    // zzzfff->Add("spair", atype<Pair<double>*>());
    Dcl_Type<Pair<double>*>(InitP<Pair<double> >, Destroy<Pair<double> >);
    zzzfff->Add("dpair", atype<Pair<double>*>());
    // Dcl_Type<Pair<std::complex<float> >*>(InitP<Pair<std::complex<float> > >, Destroy<Pair<std::complex<float> > >);
    // zzzfff->Add("cpair", atype<Pair<std::complex<float> >*>());
    Dcl_Type<Pair<std::complex<double> >*>(InitP<Pair<std::complex<double> > >, Destroy<Pair<std::complex<double> > >);
    zzzfff->Add("zpair", atype<Pair<std::complex<double> >*>());
