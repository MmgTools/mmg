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
/* Includes                                                                                                     */
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "libmesh5.h"


/*----------------------------------------------------------*/
/* Structures                                                                                           */
/*----------------------------------------------------------*/

typedef struct
{
  int typ, SolSiz, NmbLin, NmbTyp, TypTab[ GmfMaxTyp ];
  long pos;
  char fmt[ GmfMaxTyp ];
}KwdSct;

typedef struct
{
  int dim, ver, iter, mod, typ, cod;
  long NexKwdPos;
  double angle, bbox[3][2], time;
  KwdSct KwdTab[ GmfMaxKwd + 1 ];
  FILE *hdl;
  char FilNam[ GmfStrSiz ];
}GmfMshSct;


/*----------------------------------------------------------*/
/* Defines                                                                                                      */
/*----------------------------------------------------------*/

#define Asc 1
#define Bin 2
#define MshFil 4
#define SolFil 8
#define MaxMsh 100
#define InfKwd 1
#define RegKwd 2
#define SolKwd 3
#define WrdSiz 4


/*----------------------------------------------------------*/
/* Global variables                                                                                     */
/*----------------------------------------------------------*/

int GmfIniFlg=0;
GmfMshSct *GmfMshTab[ MaxMsh + 1 ];
char *GmfKwdFmt[ GmfMaxKwd + 1 ][4] =
  {     {"Reserved", "", "", ""},
        {"MeshVersionFormatted", "", "", "i"},
        {"Reserved", "", "", ""},
        {"Dimension", "", "", "i"},
        {"Vertices", "Vertex", "i", "dri"},
        {"Edges", "Edge", "i", "iii"},
        {"Triangles", "Triangle", "i", "iiii"},
        {"Quadrilaterals", "Quadrilateral", "i", "iiiii"},
        {"Tetrahedra", "Tetrahedron", "i", "iiiii"},
        {"Pentahedra", "Pentahedron", "i", "iiiiiii"},
        {"Hexahedra", "Hexahedron", "i", "iiiiiiiii"},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Corners", "Corner", "i", "i"},
        {"Ridges", "Ridge", "i", "i"},
        {"RequiredVertices", "RequiredVertex", "i", "i"},
        {"RequiredEdges", "RequiredEdge", "i", "i"},
        {"RequiredTriangles", "RequiredTriangle", "i", "i"},
        {"RequiredQuadrilaterals", "RequiredQuadrilateral", "i", "i"},
        {"TangentAtEdgeVertices", "TangentAtEdgeVertex", "i", "iii"},
        {"NormalAtVertices", "NormalAtVertex", "i", "ii"},
        {"NormalAtTriangleVertices", "NormalAtTriangleVertex", "i", "iii"},
        {"NormalAtQuadrilateralVertices", "NormalAtQuadrilateralVertex", "i", "iiii"},
        {"AngleOfCornerBound", "", "", "r"},
        {"TrianglesP2", "TriangleP2", "i", "iiiiiii"},
        {"TrianglesP3", "TriangleP3", "i", "iiiiiiiiii"},
        {"TrianglesP4", "TriangleP4", "i", "iiiiiiiiiiiii"},
        {"QuadrilateralsP2", "QuadrilateralP2", "i", "iiiiiiiii"},
        {"QuadrilateralsP3", "QuadrilateralP3", "i", "iiiiiiiiiiiii"},
        {"QuadrilateralsP4", "QuadrilateralP4", "i", "iiiiiiiiiiiiiiiii"},
        {"TetrahedraP2", "TetrahedronP2", "i", "iiiiiiiiii"},
        {"TetrahedraP3", "TetrahedronP3", "i", "iiiiiiiiiiiiiiii"},
        {"TetrahedraP4", "TetrahedronP4", "i", "iiiiiiiiiiiiiiiiiiiiii"},
        {"HexahedraP2", "HexahedronP2", "i", "iiiiiiiiiiiiiiiiiiiii"},
        {"HexahedraP3", "HexahedronP3", "i", "iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii"},
        {"HexahedraP4", "HexahedronP4", "i", "iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii"},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"BoundingBox", "", "", "drdr"},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"End", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Reserved", "", "", ""},
        {"Tangents", "Tangent", "i", "dr"},
        {"Normals", "Normal", "i", "dr"},
        {"TangentAtVertices", "TangentAtVertex", "i", "ii"},
        {"SolAtVertices", "SolAtVertex", "i", "sr"},
        {"SolAtEdges", "SolAtEdge", "i", "sr"},
        {"SolAtTriangles", "SolAtTriangle", "i", "sr"},
        {"SolAtQuadrilaterals", "SolAtQuadrilateral", "i", "sr"},
        {"SolAtTetrahedra", "SolAtTetrahedron", "i", "sr"},
        {"SolAtPentahedra", "SolAtPentahedron", "i", "sr"},
        {"SolAtHexahedra", "SolAtHexahedron", "i", "sr"},
        {"DSolAtVertices", "DSolAtVertex", "i", "sr"},
        {"ISolAtVertices", "ISolAtVertex", "i", "i"},
        {"ISolAtEdges", "ISolAtEdge", "i", "ii"},
        {"ISolAtTriangles", "ISolAtTriangle", "i", "iii"},
        {"ISolAtQuadrilaterals", "ISolAtQuadrilateral", "i", "iiii"},
        {"ISolAtTetrahedra", "ISolAtTetrahedron", "i", "iiii"},
        {"ISolAtPentahedra", "ISolAtPentahedron", "i", "iiiiii"},
        {"ISolAtHexahedra", "ISolAtHexahedron", "i", "iiiiiiii"},
        {"Iterations", "","","i"},
        {"Time", "","","r"},
        {"Reserved", "","",""}
  };


/*----------------------------------------------------------*/
/* Prototypes of local procedures                                                       */
/*----------------------------------------------------------*/

static void ScaWrd(GmfMshSct *, unsigned char *);
static void ScaDblWrd(GmfMshSct *, unsigned char *);
static void ScaBlk(GmfMshSct *, unsigned char *, int);
static long GetPos(GmfMshSct *);
static void RecWrd(GmfMshSct *, unsigned char *);
static void RecDblWrd(GmfMshSct *, unsigned char *);
static void RecBlk(GmfMshSct *, unsigned char *, int);
static void SetPos(GmfMshSct *, long);
static int ScaKwdTab(GmfMshSct *);
static void ExpFmt(GmfMshSct *, int);
static void ScaKwdHdr(GmfMshSct *, int);


/*----------------------------------------------------------*/
/* Open a mesh file in read or write mod                                        */
/*----------------------------------------------------------*/

int GmfOpenMesh(char *FilNam, int mod, ...)
{
  int i, KwdCod, res, *PtrVer, *PtrDim, MshIdx=0;
  char str[ GmfStrSiz ];
  va_list VarArg;
  GmfMshSct *msh;

  if(!GmfIniFlg)
    {
      for(i=0;i<MaxMsh;i++)
        GmfMshTab[i] = NULL;

      GmfIniFlg = 1;
    }

  /*---------------------*/
  /* MESH STRUCTURE INIT */
  /*---------------------*/

  for(i=1;i<MaxMsh;i++)
    if(!GmfMshTab[i])
      {
        MshIdx = i;
        break;
      }

  if( !MshIdx || !(msh = calloc(1, sizeof(GmfMshSct))) )
    return(0);

  /* Copy the FilNam into the structure */

  if(strlen(FilNam) + 7 >= GmfStrSiz)
    return(0);

  strcpy(msh->FilNam, FilNam);

  /* Store the opening mod (read or write) and guess the filetype (binary or ascii) depending on the extension */

  msh->mod = mod;

  if(strstr(msh->FilNam, ".meshb"))
    msh->typ |= (Bin | MshFil);
  else if(strstr(msh->FilNam, ".mesh"))
    msh->typ |= (Asc | MshFil);
  else if(strstr(msh->FilNam, ".solb"))
    msh->typ |= (Bin | SolFil);
  else if(strstr(msh->FilNam, ".sol"))
    msh->typ |= (Asc | SolFil);
  else
    return(0);

  /* Open the file in the required mod and initialyse the mesh structure */

  if(msh->mod == GmfRead)
    {

      /*-----------------------*/
      /* OPEN FILE FOR READING */
      /*-----------------------*/

      va_start(VarArg, mod);
      PtrVer = va_arg(VarArg, int *);
      PtrDim = va_arg(VarArg, int *);
      va_end(VarArg);

      /* Create the name string and open the file */

      if(!(msh->hdl = fopen(msh->FilNam, "rb")))
        return(0);

      /* Read the endian coding tag, the mesh version and the mesh dimension (mandatory kwd) */

      if(msh->typ & Bin)
        {
          fread((unsigned char *)&msh->cod, WrdSiz, 1, msh->hdl);

          if( (msh->cod != 1) && (msh->cod != 16777216) )
            return(0);

          ScaWrd(msh, (unsigned char *)&msh->ver);

          if( (msh->ver < 1) || (msh->ver > 3) )
            return(0);

          if( (msh->ver == 3) && (sizeof(long) == 4) )
            return(0);

          ScaWrd(msh, (unsigned char *)&KwdCod);

          if(KwdCod != GmfDimension)
            return(0);

          GetPos(msh);
          ScaWrd(msh, (unsigned char *)&msh->dim);
        }
      else
        {
          do
            {
              res = fscanf(msh->hdl, "%s", str);
            }while( (res != EOF) && strcmp(str, "MeshVersionFormatted") );

          if(res == EOF)
            return(0);

          fscanf(msh->hdl, "%d", &msh->ver);

          if( (msh->ver < 1) || (msh->ver > 3) )
            return(0);

          do
            {
              res = fscanf(msh->hdl, "%s", str);
            }while( (res != EOF) && strcmp(str, "Dimension") );

          if(res == EOF)
            return(0);

          fscanf(msh->hdl, "%d", &msh->dim);
        }

      if( (msh->dim != 2) && (msh->dim != 3) )
        return(0);

      (*PtrVer) = msh->ver;
      (*PtrDim) = msh->dim;

      /*------------*/
      /* KW READING */
      /*------------*/

      /* Read the list of kw present in the file */

      if(!ScaKwdTab(msh))
        return(0);

      GmfMshTab[ MshIdx ] = msh;

      return(MshIdx);
    }
  else if(msh->mod == GmfWrite)
    {

      /*-----------------------*/
      /* OPEN FILE FOR WRITING */
      /*-----------------------*/

      msh->cod = 1;

      /* Check if the user provided a valid version number and dimension */

      va_start(VarArg, mod);
      msh->ver = va_arg(VarArg, int);
      msh->dim = va_arg(VarArg, int);
      va_end(VarArg);

      if( (msh->ver < 1) || (msh->ver > 3) )
        return(0);

      if( (msh->ver == 3) && (sizeof(long) == 4) )
        return(0);

      if( (msh->dim != 2) && (msh->dim != 3) )
        return(0);

      /* Create the mesh file */

      if(!(msh->hdl = fopen(msh->FilNam, "wb")))
        return(0);

      GmfMshTab[ MshIdx ] = msh;


      /*------------*/
      /* KW WRITING */
      /*------------*/

      /* Write the mesh version and dimension */

      if(msh->typ & Asc)
        {
          fprintf(msh->hdl, "%s %d\n\n", GmfKwdFmt[ GmfVersionFormatted ][0], msh->ver);
          fprintf(msh->hdl, "%s %d\n", GmfKwdFmt[ GmfDimension ][0], msh->dim);
        }
      else
        {
          RecWrd(msh, (unsigned char *)&msh->cod);
          RecWrd(msh, (unsigned char *)&msh->ver);
          GmfSetKwd(MshIdx, GmfDimension, 0);
          RecWrd(msh, (unsigned char *)&msh->dim);
        }

      return(MshIdx);
    }
  else
    return(0);
}


/*----------------------------------------------------------*/
/* Close a meshfile in the right way                                            */
/*----------------------------------------------------------*/

int GmfCloseMesh(int MshIdx)
{
  int res = 1;
  GmfMshSct *msh;

  if( (MshIdx < 1) || (MshIdx > MaxMsh) )
    return(0);

  msh = GmfMshTab[ MshIdx ];

  /* In write down the "End" kw in write mode */

  if(msh->mod == GmfWrite){
    if(msh->typ & Asc)
      fprintf(msh->hdl, "\n%s\n", GmfKwdFmt[ GmfEnd ][0]);
    else
      GmfSetKwd(MshIdx, GmfEnd, 0);
  }
  /* Close the file and free the mesh structure */

  if(fclose(msh->hdl))
    res = 0;

  free(msh);
  GmfMshTab[ MshIdx ] = NULL;

  return(res);
}


/*----------------------------------------------------------*/
/* Read the number of lines and set the position to this kwd*/
/*----------------------------------------------------------*/

int GmfStatKwd(int MshIdx, int KwdCod, ...)
{
  int i, *PtrNmbTyp, *PtrSolSiz, *TypTab;
  GmfMshSct *msh;
  KwdSct *kwd;
  va_list VarArg;

  if( (MshIdx < 1) || (MshIdx > MaxMsh) )
    return(0);

  msh = GmfMshTab[ MshIdx ];

  if( (KwdCod < 1) || (KwdCod > GmfMaxKwd) )
    return(0);

  kwd = &msh->KwdTab[ KwdCod ];

  if(!kwd->NmbLin)
    return(0);

  /* Read further arguments if this kw is a sol */

  if(kwd->typ == SolKwd)
    {
      va_start(VarArg, KwdCod);

      PtrNmbTyp = va_arg(VarArg, int *);
      *PtrNmbTyp = kwd->NmbTyp;

      PtrSolSiz = va_arg(VarArg, int *);
      *PtrSolSiz = kwd->SolSiz;

      TypTab = va_arg(VarArg, int *);

      for(i=0;i<kwd->NmbTyp;i++)
        TypTab[i] = kwd->TypTab[i];

      va_end(VarArg);
    }

  return(kwd->NmbLin);
}


/*----------------------------------------------------------*/
/* Set the current file position to a given kwd                         */
/*----------------------------------------------------------*/

int GmfGotoKwd(int MshIdx, int KwdCod)
{
  GmfMshSct *msh;
  KwdSct *kwd;

  if( (MshIdx < 1) || (MshIdx > MaxMsh) )
    return(0);

  msh = GmfMshTab[ MshIdx ];

  if( (KwdCod < 1) || (KwdCod > GmfMaxKwd) )
    return(0);

  kwd = &msh->KwdTab[ KwdCod ];

  if(!kwd->NmbLin)
    return(0);

  return(fseek(msh->hdl, kwd->pos, SEEK_SET));
}


/*----------------------------------------------------------*/
/* Write the kwd and set the number of lines                            */
/*----------------------------------------------------------*/

int GmfSetKwd(int MshIdx, int KwdCod, ...)
{
  int i, NmbLin=0, *TypTab;
  long CurPos;
  va_list VarArg;
  GmfMshSct *msh;
  KwdSct *kwd;

  if( (MshIdx < 1) || (MshIdx > MaxMsh) )
    return(0);

  msh = GmfMshTab[ MshIdx ];

  if( (KwdCod < 1) || (KwdCod > GmfMaxKwd) )
    return(0);

  kwd = &msh->KwdTab[ KwdCod ];

  /* Read further arguments if this kw has a header */

  if(strlen(GmfKwdFmt[ KwdCod ][2]))
    {
      va_start(VarArg, KwdCod);
      NmbLin = va_arg(VarArg, int);

      if(!strcmp(GmfKwdFmt[ KwdCod ][3], "sr"))
        {
          kwd->NmbTyp = va_arg(VarArg, int);
          TypTab = va_arg(VarArg, int *);

          for(i=0;i<kwd->NmbTyp;i++)
            kwd->TypTab[i] = TypTab[i];
        }

      va_end(VarArg);
    }

  /* Setup the kwd info */

  ExpFmt(msh, KwdCod);

  if(!kwd->typ)
    return(0);
  else if(kwd->typ == InfKwd)
    kwd->NmbLin = 1;
  else
    kwd->NmbLin = NmbLin;

  /* Store the next kwd position in binary file */

  if( (msh->typ & Bin) && msh->NexKwdPos )
    {
      CurPos = ftell(msh->hdl);
      fseek(msh->hdl, msh->NexKwdPos, SEEK_SET);
      SetPos(msh, CurPos);
      fseek(msh->hdl, CurPos, SEEK_SET);
    }

  /* Write the header */

  if(msh->typ & Asc)
    {
      fprintf(msh->hdl, "\n%s\n", GmfKwdFmt[ KwdCod ][0]);

      if(kwd->typ != InfKwd)
        fprintf(msh->hdl, "%d\n", kwd->NmbLin);

      /* In case of solution field, write the extended header */

      if(kwd->typ == SolKwd)
        {
          fprintf(msh->hdl, "%d ", kwd->NmbTyp);

          for(i=0;i<kwd->NmbTyp;i++)
            fprintf(msh->hdl, "%d ", kwd->TypTab[i]);

          fprintf(msh->hdl, "\n\n");
        }
    }
  else
    {
      RecWrd(msh, (unsigned char *)&KwdCod);
      msh->NexKwdPos = ftell(msh->hdl);
      SetPos(msh, 0);

      if(kwd->typ != InfKwd)
        RecWrd(msh, (unsigned char *)&kwd->NmbLin);

      /* In case of solution field, write the extended header at once */

      if(kwd->typ == SolKwd)
        {
          RecWrd(msh, (unsigned char *)&kwd->NmbTyp);

          for(i=0;i<kwd->NmbTyp;i++)
            RecWrd(msh, (unsigned char *)&kwd->TypTab[i]);
        }
    }

  return(kwd->NmbLin);
}


/*----------------------------------------------------------*/
/* Read a full line from the current kwd                                        */
/*----------------------------------------------------------*/

void GmfGetLin(int MshIdx, int KwdCod, ...)
{
  unsigned char buf[ GmfMaxTyp * WrdSiz ];
  int i, j, *IntPtr, *IntBuf = (int *)buf;
  float *FltPtr, *FltSolTab, *FltBuf = (float *)buf;
  double *DblPtr, *DblSolTab;
  va_list VarArg;
  GmfMshSct *msh = GmfMshTab[ MshIdx ];
  KwdSct *kwd = &msh->KwdTab[ KwdCod ];

  /* Start decoding the arguments */

  va_start(VarArg, KwdCod);

  if(kwd->typ != SolKwd)
    {
      if(msh->ver == 1)
        {
          if(msh->typ & Asc)
            {
              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      FltPtr = va_arg(VarArg, float *);
                      fscanf(msh->hdl, "%f", FltPtr);
                    }
                  else
                    {
                      IntPtr = va_arg(VarArg, int *);
                      fscanf(msh->hdl, "%d", IntPtr);
                    }
                }
            }
          else
            {
              ScaBlk(msh, buf, kwd->SolSiz);

              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      FltPtr = va_arg(VarArg, float *);
                      *FltPtr = FltBuf[i];
                    }
                  else
                    {
                      IntPtr = va_arg(VarArg, int *);
                      *IntPtr = IntBuf[i];
                    }
                }
            }
        }
      else
        {
          if(msh->typ & Asc)
            {
              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      DblPtr = va_arg(VarArg, double *);
                      fscanf(msh->hdl, "%lf", DblPtr);
                    }
                  else
                    {
                      IntPtr = va_arg(VarArg, int *);
                      fscanf(msh->hdl, "%d", IntPtr);
                    }
                }
            }
          else
            {
              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      DblPtr = va_arg(VarArg, double *);
                      ScaDblWrd(msh, (unsigned char *)DblPtr);
                    }
                  else
                    {
                      IntPtr = va_arg(VarArg, int *);
                      ScaWrd(msh, (unsigned char *)IntPtr);
                    }
                }
            }
        }
    }
  else
    {
      if(msh->ver == 1)
        {
          FltSolTab = va_arg(VarArg, float *);

          if(msh->typ & Asc)
            for(j=0;j<kwd->SolSiz;j++)
              fscanf(msh->hdl, "%f", &FltSolTab[j]);
          else
            for(j=0;j<kwd->SolSiz;j++)
              ScaWrd(msh, (unsigned char *)&FltSolTab[j]);
        }
      else if(msh->ver == 2)
        {
          DblSolTab = va_arg(VarArg, double *);

          if(msh->typ & Asc)
            for(j=0;j<kwd->SolSiz;j++)
              fscanf(msh->hdl, "%lf", &DblSolTab[j]);
          else
            for(j=0;j<kwd->SolSiz;j++)
              ScaDblWrd(msh, (unsigned char *)&DblSolTab[j]);
        }
    }

  va_end(VarArg);
}


/*----------------------------------------------------------*/
/* Write a full line from the current kwd                                       */
/*----------------------------------------------------------*/

void GmfSetLin(int MshIdx, int KwdCod, ...)
{
  unsigned char buf[ GmfMaxTyp * WrdSiz ];
  int i, j, *IntBuf = (int *)buf;
  float *FltSolTab, *FltBuf = (float *)buf;
  double d, *DblSolTab;
  va_list VarArg;
  GmfMshSct *msh = GmfMshTab[ MshIdx ];
  KwdSct *kwd = &msh->KwdTab[ KwdCod ];

  /* Start decoding the arguments */

  va_start(VarArg, KwdCod);

  if(kwd->typ != SolKwd)
    {
      if(msh->ver == 1)
        {
          if(msh->typ & Asc)
            {
              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      d = va_arg(VarArg, double);
                      fprintf(msh->hdl, "%g ", (float)d);
                    }
                  else
                    {
                      j = va_arg(VarArg, int);
                      fprintf(msh->hdl, "%d ", j);
                    }
                }
            }
          else
            {
              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      d = va_arg(VarArg, double);
                      FltBuf[i] = d;
                    }
                  else
                    {
                      j = va_arg(VarArg, int);
                      IntBuf[i] = j;
                    }
                }

              RecBlk(msh, buf, kwd->SolSiz);
            }
        }
      else
        {
          if(msh->typ & Asc)
            {
              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      d = va_arg(VarArg, double);
                      fprintf(msh->hdl, "%.15lg ", d);
                    }
                  else
                    {
                      j = va_arg(VarArg, int);
                      fprintf(msh->hdl, "%d ", j);
                    }
                }
            }
          else
            {
              for(i=0;i<kwd->SolSiz;i++)
                {
                  if(kwd->fmt[i] == 'r')
                    {
                      d = va_arg(VarArg, double);
                      RecDblWrd(msh, (unsigned char *)&d);
                    }
                  else
                    {
                      j = va_arg(VarArg, int);
                      RecWrd(msh, (unsigned char *)&j);
                    }
                }
            }
        }
    }
  else
    {
      if(msh->ver == 1)
        {
          FltSolTab = va_arg(VarArg, float *);

          if(msh->typ & Asc)
            for(j=0;j<kwd->SolSiz;j++)
              fprintf(msh->hdl, "%g ", FltSolTab[j]);
          else
            for(j=0;j<kwd->SolSiz;j++)
              RecWrd(msh, (unsigned char *)&FltSolTab[j]);
        }
      else if(msh->ver == 2)
        {
          DblSolTab = va_arg(VarArg, double *);

          if(msh->typ & Asc)
            for(j=0;j<kwd->SolSiz;j++)
              fprintf(msh->hdl, "%.15lg ", DblSolTab[j]);
          else
            for(j=0;j<kwd->SolSiz;j++)
              RecDblWrd(msh, (unsigned char *)&DblSolTab[j]);
        }
    }

  va_end(VarArg);

  if(msh->typ & Asc)
    fprintf(msh->hdl, "\n");
}


/*----------------------------------------------------------*/
/* Private procedure for transmesh : copy a whole line          */
/*----------------------------------------------------------*/

void GmfCpyLin(int InpIdx, int OutIdx, int KwdCod)
{
  double d;
  float f;
  int i, a;
  GmfMshSct *InpMsh = GmfMshTab[ InpIdx ], *OutMsh = GmfMshTab[ OutIdx ];
  KwdSct *kwd = &InpMsh->KwdTab[ KwdCod ];

  for(i=0;i<kwd->SolSiz;i++)
    {
      if(kwd->fmt[i] == 'r')
        {
          if(InpMsh->ver == 1)
            {
              if(InpMsh->typ & Asc)
                fscanf(InpMsh->hdl, "%f", &f);
              else
                ScaWrd(InpMsh, (unsigned char *)&f);

              d = f;
            }
          else
            {
              if(InpMsh->typ & Asc)
                fscanf(InpMsh->hdl, "%lf", &d);
              else
                ScaDblWrd(InpMsh, (unsigned char *)&d);

              f = (float)d;
            }

          if(OutMsh->ver == 1)
            if(OutMsh->typ & Asc)
              fprintf(OutMsh->hdl, "%g ", f);
            else
              RecWrd(OutMsh, (unsigned char *)&f);
          else
            if(OutMsh->typ & Asc)
              fprintf(OutMsh->hdl, "%.15g ", d);
            else
              RecDblWrd(OutMsh, (unsigned char *)&d);
        }
      else
        {
          if(InpMsh->typ & Asc)
            fscanf(InpMsh->hdl, "%d", &a);
          else
            ScaWrd(InpMsh, (unsigned char *)&a);

          if(OutMsh->typ & Asc)
            fprintf(OutMsh->hdl, "%d ", a);
          else
            RecWrd(OutMsh, (unsigned char *)&a);
        }
    }

  if(OutMsh->typ & Asc)
    fprintf(OutMsh->hdl, "\n");
}


/*----------------------------------------------------------*/
/* Find every kw present in a meshfile                                          */
/*----------------------------------------------------------*/

static int ScaKwdTab(GmfMshSct *msh)
{
  int KwdCod;
  long  NexPos, CurPos, EndPos;
  char str[ GmfStrSiz ];

  if(msh->typ & Asc)
    {
      /* Scan each string in the file until the end */

      while(fscanf(msh->hdl, "%s", str) != EOF)
        {
          /* Fast test in order to reject quickly the numeric values */

          if(isalpha(str[0]))
            {
              /* Search which kwd code this string is associated with,
                 then get its header and save the curent position in file (just before the data) */

              for(KwdCod=1; KwdCod<= GmfMaxKwd; KwdCod++)
                if(!strcmp(str, GmfKwdFmt[ KwdCod ][0]))
                  {
                    ScaKwdHdr(msh, KwdCod);
                    break;
                  }
            }
          else if(str[0] == '#')
            while(fgetc(msh->hdl) != '\n');
        }
    }
  else
    {
      /* Get file size */

      CurPos = ftell(msh->hdl);
      fseek(msh->hdl, 0, SEEK_END);
      EndPos = ftell(msh->hdl);
      fseek(msh->hdl, CurPos, SEEK_SET);

      /* Jump through kwd positions in the file */

      do
        {
          /* Get the kwd code and the next kwd position */

          ScaWrd(msh, (unsigned char *)&KwdCod);
          NexPos = GetPos(msh);

          if(NexPos > EndPos)
            return(0);

          /* Check if this kwd belongs to this mesh version */

          if( (KwdCod >= 1) && (KwdCod <= GmfMaxKwd) )
            ScaKwdHdr(msh, KwdCod);

          /* Go to the next kwd */

          if(NexPos)
            fseek(msh->hdl, NexPos, SEEK_SET);
        }while(NexPos && (KwdCod != GmfEnd));
    }

  return(1);
}


/*----------------------------------------------------------*/
/* Read and setup the keyword's header                                          */
/*----------------------------------------------------------*/

static void ScaKwdHdr(GmfMshSct *msh, int KwdCod)
{
  int i;
  KwdSct *kwd = &msh->KwdTab[ KwdCod ];

  if(!strcmp("i", GmfKwdFmt[ KwdCod ][2]))
    {
      if(msh->typ & Asc)
        fscanf(msh->hdl, "%d", &kwd->NmbLin);
      else
        ScaWrd(msh, (unsigned char *)&kwd->NmbLin);
    }
  else
    kwd->NmbLin = 1;

  if(!strcmp("sr", GmfKwdFmt[ KwdCod ][3]))
    {
      if(msh->typ & Asc)
        {
          fscanf(msh->hdl, "%d", &kwd->NmbTyp);

          for(i=0;i<kwd->NmbTyp;i++)
            fscanf(msh->hdl, "%d", &kwd->TypTab[i]);
        }
      else
        {
          ScaWrd(msh, (unsigned char *)&kwd->NmbTyp);

          for(i=0;i<kwd->NmbTyp;i++)
            ScaWrd(msh, (unsigned char *)&kwd->TypTab[i]);
        }
    }

  ExpFmt(msh, KwdCod);
  kwd->pos = ftell(msh->hdl);
}


/*----------------------------------------------------------*/
/* Expand the compacted format and compute the line size        */
/*----------------------------------------------------------*/

static void ExpFmt(GmfMshSct *msh, int KwdCod)
{
  int i, j, TmpSiz=0;
  char chr, *InpFmt = GmfKwdFmt[ KwdCod ][3];
  KwdSct *kwd = &msh->KwdTab[ KwdCod ];

  /* Set the kwd's type */

  if(!strlen(GmfKwdFmt[ KwdCod ][2]))
    kwd->typ = InfKwd;
  else if(!strcmp(InpFmt, "sr"))
    kwd->typ = SolKwd;
  else
    kwd->typ = RegKwd;

  /* Get the solution-field's size */

  if(kwd->typ == SolKwd)
    for(i=0;i<kwd->NmbTyp;i++)
      switch(kwd->TypTab[i])
        {
        case GmfSca    : TmpSiz += 1; break;
        case GmfVec    : TmpSiz += msh->dim; break;
        case GmfSymMat : TmpSiz += (msh->dim * (msh->dim+1)) / 2; break;
        case GmfMat    : TmpSiz += msh->dim * msh->dim; break;
        }

  /* Scan each character from the format string */

  i = kwd->SolSiz = 0;

  while(i < strlen(InpFmt))
    {
      chr = InpFmt[ i++ ];

      if(chr == 'd')
        {
          chr = InpFmt[i++];

          for(j=0;j<msh->dim;j++)
            kwd->fmt[ kwd->SolSiz++ ] = chr;
        }
      else if(chr == 's')
        {
          chr = InpFmt[i++];

          for(j=0;j<TmpSiz;j++)
            kwd->fmt[ kwd->SolSiz++ ] = chr;
        }
      else
        kwd->fmt[ kwd->SolSiz++ ] = chr;
    }
}


/*----------------------------------------------------------*/
/* Generate fortran api automatically                                           */
/*----------------------------------------------------------*/

#ifdef f77api
int main()
{
  int i, j, k, l, m, DimMod, FltMod;
  char *ActStr[2][2] = { {"Get","Set"}, {"", "*"} };
  char *DimStr[2][2] = { {"", ""}, {"2d", "3d"} };
  char *FltStr[4][2] = { {"", ""}, {"R4", "R8"}, {"float", "double"}, {"real(4)", "real(8)"} };
  char PrcNam[256], PrcPar[256], F77Nam[256], F77Par[256];
  KwdSct *kwd;
  FILE *hdl, *F90Hdl, *F77Hdl;
  GmfMshSct msh;

  /* Create ".c" and a ".h" files */

  hdl = fopen("libmesh5_fortran_api.c", "w");

  system("cp head_f90 M_libmesh5_api.f90");
  F90Hdl = fopen("M_libmesh5_api.f90", "a");
  F77Hdl = fopen("libmesh5_api.ins", "w");

  /* Write headers */

  fprintf(hdl, "#include \"libmesh5.h\"\n\n");
  fprintf(hdl, "/* Generated automatically by libmesh5.2 */\n\n");

  fprintf(F77Hdl, "c Generated automatically by libmesh5.2 \n\n");

  fprintf(F77Hdl, "      external GmfOpenMeshF77\n");
  fprintf(F77Hdl, "      external GmfCloseMeshF77\n");
  fprintf(F77Hdl, "      external GmfStatKwdF77\n");
  fprintf(F77Hdl, "      external GmfSetKwdF77\n");
  fprintf(F77Hdl, "      external GmfGotoKwdF77\n");
  fprintf(F77Hdl, "\n");

  fprintf(F77Hdl, "      integer GmfOpenMeshF77\n");
  fprintf(F77Hdl, "      integer GmfCloseMeshF77\n");
  fprintf(F77Hdl, "      integer GmfStatKwdF77\n");
  fprintf(F77Hdl, "      integer GmfSetKwdF77\n");
  fprintf(F77Hdl, "      integer GmfGotoKwdF77\n");
  fprintf(F77Hdl, "\n");

  /* Start generating keyword's prototypes */

  for(i=1;i<=GmfMaxKwd;i++)
    {
      if(strcmp(GmfKwdFmt[i][2], "i"))
        continue;

      kwd = &msh.KwdTab[i];

      if(strchr(GmfKwdFmt[i][3], 'd'))
        DimMod = 1;
      else
        DimMod = 0;

      if(strchr(GmfKwdFmt[i][3], 'r'))
        FltMod = 1;
      else
        FltMod = 0;

      for(j=0;j<2;j++)
        for(k=0;k<=DimMod;k++)
          for(l=0;l<=FltMod;l++)
            {
              /* Generate function name's strings */

              sprintf(F77Nam, "Gmf%s%s%s%s",ActStr[0][j], GmfKwdFmt[i][1], DimStr[ DimMod ][k], FltStr[ FltMod ][l]);
              sprintf(PrcNam, "Gmf%sLin",ActStr[0][j]);

              for(m=0;m<strlen(F77Nam);m++)
                F77Nam[m] = tolower(F77Nam[m]);


              /* Write f77api.c */

              fprintf(hdl, "void call(%s)(int *MshIdx", F77Nam);

              msh.dim = k+2;
              ExpFmt(&msh, i);

              if(kwd->typ == RegKwd)
                {
                  for(m=0;m<kwd->SolSiz;m++)
                    if(kwd->fmt[m] == 'i')
                      fprintf(hdl, ", int *i%d", m);
                    else if(kwd->fmt[m] == 'r')
                      fprintf(hdl, ", %s *r%d", FltStr[2][l], m);
                }
              else if(kwd->typ == SolKwd)
                {
                  fprintf(hdl, ", %s *SolTab", FltStr[2][l]);
                }

              fprintf(hdl, ")\n");

              fprintf(hdl, "{\n %s(*MshIdx, Gmf%s",PrcNam,GmfKwdFmt[i][0]);

              if(kwd->typ == RegKwd)
                {
                  for(m=0;m<kwd->SolSiz;m++)
                    if(kwd->fmt[m] == 'i')
                      fprintf(hdl, ", %si%d", ActStr[1][j], m);
                    else if(kwd->fmt[m] == 'r')
                      fprintf(hdl, ", %sr%d", ActStr[1][j], m);
                }
              else if(kwd->typ == SolKwd)
                {
                  fprintf(hdl, ", SolTab", FltStr[2][l]);
                }

              fprintf(hdl, ");\n}\n\n");


              /* Write f77api.f90 */

              fprintf(F90Hdl, "\ninterface\n");
              fprintf(F90Hdl, "  subroutine %s(MshIdx", F77Nam);

              msh.dim = k+2;
              ExpFmt(&msh, i);

              if(kwd->typ == RegKwd)
                {
                  for(m=0;m<kwd->SolSiz;m++)
                    if(kwd->fmt[m] == 'i')
                      fprintf(F90Hdl, ", i%d",m);
                    else if(kwd->fmt[m] == 'r')
                      fprintf(F90Hdl, ", r%d", m);
                }
              else if(kwd->typ == SolKwd)
                {
                  fprintf(F90Hdl, ", SolTab");
                }

              fprintf(F90Hdl, ")\n");

              fprintf(F90Hdl, "  integer :: MshIdx\n");

              if(kwd->typ == RegKwd)
                {
                  for(m=0;m<kwd->SolSiz;m++)
                    if(kwd->fmt[m] == 'i')
                      fprintf(F90Hdl, "  integer :: i%d\n",m);
                    else if(kwd->fmt[m] == 'r')
                      fprintf(F90Hdl, "  %s :: r%d\n", FltStr[3][l],m);
                }
              else if(kwd->typ == SolKwd)
                fprintf(F90Hdl, "  %s :: SolTab(*)\n", FltStr[3][l]);

              fprintf(F90Hdl, "  end subroutine %s\n", F77Nam);
              fprintf(F90Hdl, "end interface\n");


              /* Write f77api.ins */

              fprintf(F77Hdl, "      external %s\n", F77Nam);
            }
    }

  /* Generate f90 keywords */

  fprintf(F90Hdl, "\n");
  fprintf(F90Hdl, "integer,parameter :: GmfMaxTyp=%d\n", GmfMaxTyp);
  fprintf(F90Hdl, "integer,parameter :: GmfMaxKwd=%d\n", GmfMaxKwd);
  fprintf(F90Hdl, "integer,parameter :: GmfRead=%d\n", GmfRead);
  fprintf(F90Hdl, "integer,parameter :: GmfWrite=%d\n", GmfWrite);
  fprintf(F90Hdl, "integer,parameter :: GmfSca=%d\n", GmfSca);
  fprintf(F90Hdl, "integer,parameter :: GmfVec=%d\n", GmfVec);
  fprintf(F90Hdl, "integer,parameter :: GmfSymMat=%d\n", GmfSymMat);
  fprintf(F90Hdl, "integer,parameter :: GmfMat=%d\n", GmfMat);
  fprintf(F90Hdl, "\n");

  for(i=1;i<=GmfMaxKwd;i++)
    if(strcmp(GmfKwdFmt[i][0], "Reserved"))
      fprintf(F90Hdl, "integer,parameter :: Gmf%s=%d\n", GmfKwdFmt[i][0], i);

  fprintf(F90Hdl, "\nend module M_libmesh5_api\n\n");

  /* Generate f77 keywords */

  fprintf(F77Hdl, "\n");
  fprintf(F77Hdl, "      integer GmfMaxTyp\n");
  fprintf(F77Hdl, "      integer GmfMaxKwd\n");
  fprintf(F77Hdl, "      integer GmfRead\n");
  fprintf(F77Hdl, "      integer GmfWrite\n");
  fprintf(F77Hdl, "      integer GmfSca\n");
  fprintf(F77Hdl, "      integer GmfVec\n");
  fprintf(F77Hdl, "      integer GmfSymMat\n");
  fprintf(F77Hdl, "      integer GmfMat\n");
  fprintf(F77Hdl, "\n");

  fprintf(F77Hdl, "      parameter GmfMaxTyp=%d\n", GmfMaxTyp);
  fprintf(F77Hdl, "      parameter GmfMaxKwd=%d\n", GmfMaxKwd);
  fprintf(F77Hdl, "      parameter GmfRead=%d\n", GmfRead);
  fprintf(F77Hdl, "      parameter GmfWrite=%d\n", GmfWrite);
  fprintf(F77Hdl, "      parameter GmfSca=%d\n", GmfSca);
  fprintf(F77Hdl, "      parameter GmfVec=%d\n", GmfVec);
  fprintf(F77Hdl, "      parameter GmfSymMat=%d\n", GmfSymMat);
  fprintf(F77Hdl, "      parameter GmfMat=%d\n", GmfMat);
  fprintf(F77Hdl, "\n");

  for(i=1;i<=GmfMaxKwd;i++)
    if(strcmp(GmfKwdFmt[i][0], "Reserved"))
      fprintf(F77Hdl, "      integer Gmf%s\n", GmfKwdFmt[i][0]);

  fprintf(F77Hdl, "\n");

  for(i=1;i<=GmfMaxKwd;i++)
    if(strcmp(GmfKwdFmt[i][0], "Reserved"))
      fprintf(F77Hdl, "      parameter Gmf%s=%d\n", GmfKwdFmt[i][0], i);

  fclose(hdl);
  fclose(F90Hdl);
  fclose(F77Hdl);
}
#endif


/*----------------------------------------------------------*/
/* Read a four bytes word from a mesh file                                      */
/*----------------------------------------------------------*/

static void ScaWrd(GmfMshSct *msh, unsigned char *wrd)
{
  unsigned char swp;

  fread(wrd, WrdSiz, 1, msh->hdl);

  if(msh->cod == 1)
    return;

  swp = wrd[3];
  wrd[3] = wrd[0];
  wrd[0] = swp;

  swp = wrd[2];
  wrd[2] = wrd[1];
  wrd[1] = swp;
}


/*----------------------------------------------------------*/
/* Read an eight bytes word from a mesh file                            */
/*----------------------------------------------------------*/

static void ScaDblWrd(GmfMshSct *msh, unsigned char *wrd)
{
  int i;
  unsigned char swp;

  fread(wrd, WrdSiz, 2, msh->hdl);

  if(msh->cod == 1)
    return;

  for(i=0;i<4;i++)
    {
      swp = wrd[7-i];
      wrd[7-i] = wrd[i];
      wrd[i] = swp;
    }
}


/*----------------------------------------------------------*/
/* Read ablock of four bytes word from a mesh file                      */
/*----------------------------------------------------------*/

static void ScaBlk(GmfMshSct *msh, unsigned char *blk, int siz)
{
  int i, j;
  unsigned char swp, *wrd;

  fread(blk, WrdSiz, siz, msh->hdl);

  if(msh->cod == 1)
    return;

  for(i=0;i<siz;i++)
    {
      wrd = &blk[ i * 4 ];

      for(j=0;j<2;j++)
        {
          swp = wrd[ 3-j ];
          wrd[ 3-j ] = wrd[j];
          wrd[j] = swp;
        }
    }
}


/*----------------------------------------------------------*/
/* Read a 4 or 8 bytes position in mesh file                            */
/*----------------------------------------------------------*/

static long GetPos(GmfMshSct *msh)
{
  int IntVal;
  long pos;

  if(msh->ver >= 3)
    ScaDblWrd(msh, (unsigned char*)&pos);
  else
    {
      ScaWrd(msh, (unsigned char*)&IntVal);
      pos = IntVal;
    }

  return(pos);
}


/*----------------------------------------------------------*/
/* Write a four bytes word to a mesh file                                       */
/*----------------------------------------------------------*/

static void RecWrd(GmfMshSct *msh, unsigned char *wrd)
{
  fwrite(wrd, WrdSiz, 1, msh->hdl);
}


/*----------------------------------------------------------*/
/* Write an eight bytes word to a mesh file                                     */
/*----------------------------------------------------------*/

static void RecDblWrd(GmfMshSct *msh, unsigned char *wrd)
{
  fwrite(wrd, WrdSiz, 2, msh->hdl);
}


/*----------------------------------------------------------*/
/* Write a block of four bytes word to a mesh file                      */
/*----------------------------------------------------------*/

static void RecBlk(GmfMshSct *msh, unsigned char *blk, int siz)
{
  fwrite(blk, WrdSiz, siz, msh->hdl);
}


/*----------------------------------------------------------*/
/* Read a 4 or 8 bytes position in mesh file                            */
/*----------------------------------------------------------*/

static void SetPos(GmfMshSct *msh, long pos)
{
  int IntVal;

  if(msh->ver >= 3)
    RecDblWrd(msh, (unsigned char*)&pos);
  else
    {
      IntVal = pos;
      RecWrd(msh, (unsigned char*)&IntVal);
    }
}


/*----------------------------------------------------------*/
/* Fortran 77 API                                                                                       */
/*----------------------------------------------------------*/

int call(gmfopenmeshf77)(char *FilNam, int *mod, int *ver, int *dim, int StrSiz)
{
  int i;
  char TmpNam[ GmfStrSiz ];

  for(i=0;i<StrSiz;i++)
    TmpNam[i] = FilNam[i];

  TmpNam[ StrSiz ] = 0;

  if(*mod == GmfRead)
    return(GmfOpenMesh(TmpNam, *mod, ver, dim));
  else
    return(GmfOpenMesh(TmpNam, *mod, *ver, *dim));
}

int call(gmfclosemeshf77)(int *idx)
{
  return(GmfCloseMesh(*idx));
}

int call(gmfstatkwdf77)(int *MshIdx, int *KwdIdx, int *NmbTyp, int *SolSiz, int *TypTab)
{
  if(!strcmp(GmfKwdFmt[ *KwdIdx ][3], "sr"))
    return(GmfStatKwd(*MshIdx, *KwdIdx, NmbTyp, SolSiz, TypTab));
  else
    return(GmfStatKwd(*MshIdx, *KwdIdx));
}

int call(gmfgotokwdf77)(int *MshIdx, int *KwdIdx)
{
  return(GmfGotoKwd(*MshIdx, *KwdIdx));
}

int call(gmfsetkwdf77)(int *MshIdx, int *KwdIdx, int *NmbLin, int *NmbTyp, int *TypTab)
{
  if(!strcmp(GmfKwdFmt[ *KwdIdx ][3], "sr"))
    return(GmfSetKwd(*MshIdx, *KwdIdx, *NmbLin, *NmbTyp, TypTab));
  else if(strlen(GmfKwdFmt[ *KwdIdx ][2]))
    return(GmfSetKwd(*MshIdx, *KwdIdx, *NmbLin));
  else
    return(GmfSetKwd(*MshIdx, *KwdIdx));
}
