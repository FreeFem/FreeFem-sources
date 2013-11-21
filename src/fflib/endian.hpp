#ifndef ENDIAN_HPP__
#define  ENDIAN_HPP__
//  --------------------------------------------
//  read and write without endianess problem.
//  the choise is little endian , such than the 
//  the order of the bytes is the same order of
//  a shift operator <<  in a integer.   
//  F. Hecht.
//  --------------------------------------------
template<int L>
inline void w_endian(const unsigned char * c,unsigned char * k)
{cerr<< " L = "<< L  << endl; assert(0); }
template<int L>
inline void r_endian(const unsigned char *c,unsigned char * k)
{assert(0); }

template<>
inline  void w_endian<1>(const unsigned char *c,unsigned char* k)
{*k=*c;}

template<>
inline  void w_endian<8>(const unsigned char *c,unsigned char* k)
{
  static  const unsigned long long  ull =( 0ULL + (1ULL << 8) + (2ULL << 16) + (3ULL << 24)
					   + (4ULL << 32) + (5ULL << 40) + (6ULL << 48 ) + (7ULL << 56 ));
  static  const unsigned char *o = reinterpret_cast<const unsigned char *>(&ull);
  //  cout << " ++++ " << ull << " .. "   ;
  // for(int i=0;i<8;++i) cout <<(int) o[i] ;
  assert(8==sizeof(ull));
  k[o[0]]=c[0];
  k[o[1]]=c[1];
  k[o[2]]=c[2];
  k[o[3]]=c[3];  
  k[o[4]]=c[4];
  k[o[5]]=c[5];
  k[o[6]]=c[6];
  k[o[7]]=c[7];  
  // cout << " ---  " ;
  //  for (int i=0;i<8;++i) 
  //  cout << k[i] ;
  // cout << k << endl;
}

/*
template<class T> std::ostream &dump(std::ostream & f, const T &t)
{
  const unsigned char *o = reinterpret_cast<const unsigned char *>(&t);
  //  f << " " << t <<  " " << & t << " == " << (void *) o << " -> ";
  for( int i=0;i< sizeof(t); ++i)
    f  << (unsigned int) o[i] ;
  return f;
}
*/

template<>
inline  void w_endian<4>(const unsigned char *c,unsigned char * k)
{
  static  const unsigned int u0123 = 0U + (1U << 8) + (2U << 16) + (3U << 24);
  static  const unsigned char *o = reinterpret_cast<const unsigned char *>(&u0123);
  //   ordre   c[o^-1] :   si o:   0123 -> 2130  
  //   si  c == 1U + (2U << 8) + (3U << 16) + (4U << 24);   
  //   alors  k[0] = 1, k[1]=2, k[2] = 3, k[3]=4 
  //  ------------------
  assert(4==sizeof(u0123));
  k[o[0]]=c[0];
  k[o[1]]=c[1];
  k[o[2]]=c[2];
  k[o[3]]=c[3];  
}

template<>
inline  void w_endian<2>(const unsigned char *c,unsigned char *k)
{
  static  const unsigned short int  u01 = 0U + (1U << 8);
  static  const unsigned char *o = reinterpret_cast<const unsigned char *>(&u01);
  assert(2==sizeof(u01));
  k[o[0]]=c[0];
  k[o[1]]=c[1];
}

template<>
inline  void r_endian<1>(const unsigned char *c,unsigned char *k)
{  *k=*c; }

template<>
inline  void r_endian<8>(const unsigned char *c,unsigned char *k)
{
  static  const unsigned long long  ull =( 0ULL + (1ULL << 8) + (2ULL << 16) + (3ULL << 24)
					   + (4ULL << 32) + (5ULL << 40) + (6ULL << 48 ) + (7ULL << 56 ));
  static  const unsigned char *o = reinterpret_cast<const unsigned char *>(&ull);
  assert(8==sizeof(ull));
  
  k[0]=c[o[0]];
  k[1]=c[o[1]];
  k[2]=c[o[2]];
  k[3]=c[o[3]];  
  k[4]=c[o[4]];
  k[5]=c[o[5]];
  k[6]=c[o[6]];
  k[7]=c[o[7]];  
  
}


template<>
inline  void r_endian<4>(const unsigned char c[4],unsigned char k[4])
{
  static  const unsigned int u0123 = 0U + (1U << 8) + (2U << 16) + (3U << 24);
  static  const unsigned char *o = reinterpret_cast<const unsigned char *>(&u0123);
  assert(4==sizeof(u0123));
  
  k[0]=c[o[0]];
  k[1]=c[o[1]];
  k[2]=c[o[2]];
  k[3]=c[o[3]];  
}

template<>
inline  void r_endian<2>(const unsigned char c[2],unsigned char k[2])
{
  static  const unsigned short u01 = 0U + (1U << 8) ;
  static  const unsigned char *o = reinterpret_cast<const unsigned char *>(&u01);
  assert(2==sizeof(u01));
  k[0]=c[o[0]];
  k[1]=c[o[1]];
}


template<class T>
inline  T r_endian(const T & t )
{
  T r;
  r_endian<sizeof(t)>(reinterpret_cast<const unsigned char *>(&t),reinterpret_cast<unsigned char *>(&r));
  return r;
}

template<class T>
inline  T w_endian(const T & t )
{
  T r;
  w_endian<sizeof(t)>(reinterpret_cast<const unsigned char *>(&t),reinterpret_cast<unsigned char *>(&r));
  return r;
}
#endif
