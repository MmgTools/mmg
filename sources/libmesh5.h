/*----------------------------------------------------------*/
/*                                                                                                                      */
/*                                              LIBMESH V 5.3                                           */
/*                                                                                                                      */
/*----------------------------------------------------------*/
/*                                                                                                                      */
/*      Description:            handle .meshb file format I/O           */
/*      Author:                         Loic MARECHAL                                           */
/*      Creation date:          feb 16 2007                                                     */
/*      Last modification:      dec 12 2008                                                     */
/*                                                                                                                      */
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Defines                                                                                                      */
/*----------------------------------------------------------*/

#define GmfStrSiz 1024
#define GmfMaxTyp 20
#define GmfMaxKwd 79
#define GmfMshVer 1
#define GmfRead 1
#define GmfWrite 2
#define GmfSca 1
#define GmfVec 2
#define GmfSymMat 3
#define GmfMat 4
#define GmfFloat 1
#define GmfDouble 2

enum GmfKwdCod
  {
    GmfReserved1, \
    GmfVersionFormatted, \
    GmfReserved2, \
    GmfDimension, \
    GmfVertices, \
    GmfEdges, \
    GmfTriangles, \
    GmfQuadrilaterals, \
    GmfTetrahedra, \
    GmfPentahedra, \
    GmfHexahedra, \
    GmfReserved3, \
    GmfReserved4, \
    GmfCorners, \
    GmfRidges, \
    GmfRequiredVertices, \
    GmfRequiredEdges, \
    GmfRequiredTriangles, \
    GmfRequiredQuadrilaterals, \
    GmfTangentAtEdgeVertices, \
    GmfNormalAtVertices, \
    GmfNormalAtTriangleVertices, \
    GmfNormalAtQuadrilateralVertices, \
    GmfAngleOfCornerBound, \
    GmfTrianglesP2, \
    GmfTrianglesP3, \
    GmfTrianglesP4, \
    GmfQuadrilateralsP2, \
    GmfQuadrilateralsP3, \
    GmfQuadrilateralsP4, \
    GmfTetrahedraP2, \
    GmfTetrahedraP3, \
    GmfTetrahedraP4, \
    GmfHexahedraP2, \
    GmfHexahedraP3, \
    GmfHexahedraP4, \
    GmfReserved17, \
    GmfReserved18, \
    GmfReserved19, \
    GmfReserved20, \
    GmfReserved21, \
    GmfReserved22, \
    GmfReserved23, \
    GmfReserved24, \
    GmfReserved25, \
    GmfReserved26, \
    GmfReserved27, \
    GmfReserved28, \
    GmfReserved29, \
    GmfReserved30, \
    GmfBoundingBox, \
    GmfReserved31, \
    GmfReserved32, \
    GmfReserved33, \
    GmfEnd, \
    GmfReserved34, \
    GmfReserved35, \
    GmfReserved36, \
    GmfReserved37, \
    GmfTangents, \
    GmfNormals, \
    GmfTangentAtVertices, \
    GmfSolAtVertices, \
    GmfSolAtEdges, \
    GmfSolAtTriangles, \
    GmfSolAtQuadrilaterals, \
    GmfSolAtTetrahedra, \
    GmfSolAtPentahedra, \
    GmfSolAtHexahedra, \
    GmfDSolAtVertices, \
    GmfISolAtVertices, \
    GmfISolAtEdges, \
    GmfISolAtTriangles, \
    GmfISolAtQuadrilaterals, \
    GmfISolAtTetrahedra, \
    GmfISolAtPentahedra, \
    GmfISolAtHexahedra, \
    GmfIterations, \
    GmfTime, \
    GmfReserved38
  };


/*----------------------------------------------------------*/
/* External procedures                                                                          */
/*----------------------------------------------------------*/

extern int GmfOpenMesh(char *, int, ...);
extern int GmfCloseMesh(int);
extern int GmfStatKwd(int, int, ...);
extern int GmfGotoKwd(int, int);
extern int GmfSetKwd(int, int, ...);
extern void GmfGetLin(int, int, ...);
extern void GmfSetLin(int, int, ...);


/*----------------------------------------------------------*/
/* Fortran 77 API                                                                                       */
/*----------------------------------------------------------*/

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif


/*----------------------------------------------------------*/
/* Transmesh private API                                                                        */
/*----------------------------------------------------------*/

#ifdef TRANSMESH

extern char *GmfKwdFmt[ GmfMaxKwd + 1 ][4];
extern int GmfCpyLin(int, int, int);

#endif
