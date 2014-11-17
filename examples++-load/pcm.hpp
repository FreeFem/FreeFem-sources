
/*
 * pcm.h
 *
 * Author: Dion Crannitch (adapted from Dave Mason's complexMap class)
 * Modified : 22 June 1999
 * Last Modified: 20 Oct. 2009  (F. Hecht,  Frederic.Hecht@upmc.fr)   for g++ v4 
 *
 * Declarations for PCM class, a data structure containing a 2D vector
 * image, with encapsulated operations for reading and writing PCM files.
 *
 */

#ifndef __PCM__
#define __PCM__

#include <cstdio>
using namespace std;
typedef struct { 
  float r,i;
} pcm_complex;

class PCM {
  void CalcMax();
 public:
  int  width,height;
  unsigned long pixels;
  float max;
  pcm_complex *image; 
  PCM(int w,int h);
  PCM(const char *filename);
  pcm_complex *Get(int x, int y);
  void Set(int x, int y, pcm_complex c);
  void Load(const char *filename);
  void Save(const char *filename);
  ~PCM();
private:
    PCM (const PCM &);
    void operator=(const PCM &);
};

#endif


