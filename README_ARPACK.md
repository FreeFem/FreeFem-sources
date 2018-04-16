# How to adapt ARPACK for compilation with FreeFem++

The simplest way is just 
```bash
configure --enable-download 
```

Or to compile arpack

Remark:

In arpack++, a lot of incoherance this moderne c++ (g++3 or better) so 
I write the driver by hand (from version 3). arpack++ in included in FreeFem++ now. 

The last one is in lapack lib second.f in an function not a procedure like in arpack.

Two weeks work to find this mistake.


Frederic Hecht


