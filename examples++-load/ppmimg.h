#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char ubyte;

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

enum imgtyp {DEFAULT=0, P2,P3,P4,P5,P6,
             GREY,RGB,RED,GREEN,BLUE,COLOR};
                     
typedef struct {
  short    sizeX,sizeY;
  ubyte   *data;
} PPMimage;
typedef PPMimage * pPPMimage;

/* prototypes */
PPMimage *loadPPM(const char *imgname,ubyte *type,ubyte quiet);
int       savePPM(const char *imgname,pPPMimage img,int typimg);
pPPMimage diffImg(pPPMimage bits,pPPMimage img,ubyte ityp);

#ifdef __cplusplus
}
#endif

