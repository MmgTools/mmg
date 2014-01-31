#include "mmg3d.h"

#define KA     7
#define KB    11
#define KC    13

extern char  ddb;

/** tetra packing */
static void paktet(pMesh mesh) {
  pTetra   pt,pt1;
  int      k;

  k = 1;
  do {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) {
      pt1 = &mesh->tetra[mesh->ne];
      assert(MG_EOK(pt1));
      memcpy(pt,pt1,sizeof(Tetra));
      delElt(mesh,mesh->ne);
    }
  }
  while ( ++k < mesh->ne );

  /* Recreate nil chain */
  mesh->nenil = mesh->ne + 1;

  for(k=mesh->nenil; k<=mesh->nemax-1; k++){
    mesh->tetra[k].v[3] = k+1;
  }
}

/** Create table of adjacency. *
 *  Set pack variable to 0 for a compact mesh and to 1 for *
 *  a mesh that need to be packed */
int hashTetra(pMesh mesh, int pack) {
  pTetra         pt,pt1;
  int            k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
  int           *hcode,*link,hsize,inival;
  unsigned char  i,ii,i1,i2,i3;
  unsigned int   key;

  /* default */
  if ( mesh->adja ) {
    if( !mesh->info.sing ) {
      if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
        fprintf(stdout,"  ## Warning: no re-build of adjacencies of mesh. ");
        fprintf(stdout,"mesh->adja must be freed to enforce analysis.\n");
      }
    }
    return(1);
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING STRUCTURE\n");

  /* packing : if not hash does not work */
  if ( pack )  paktet(mesh);

  /* memory alloc */
  ADD_MEM(mesh,(4*mesh->nemax+5)*sizeof(int),"adjacency table",
          printf("  Exit program.\n");
          exit(EXIT_FAILURE));
  SAFE_CALLOC(mesh->adja,4*mesh->nemax+5,int);
  SAFE_CALLOC(hcode,mesh->ne+5,int);

  link  = mesh->adja;
  hsize = mesh->ne;

  /* init */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");
  inival = 2147483647;
  iadr   = 0;
  for (k=0; k<=mesh->ne; k++)
    hcode[k] = -inival;

  /* hash tetras */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<4; i++) {
      i1 = idir[i][0];
      i2 = idir[i][1];
      i3 = idir[i][2];
      mins = MG_MIN(pt->v[i1],MG_MIN(pt->v[i2],pt->v[i3]));
      maxs = MG_MAX(pt->v[i1],MG_MAX(pt->v[i2],pt->v[i3]));

      /* compute key and insert */
      sum = pt->v[i1] + pt->v[i2] + pt->v[i3];
      key = KA*mins + KB*maxs + KC*sum;
      key = key % hsize + 1;
      iadr++;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 2: adjacencies\n");
  for (l=iadr; l>0; l--) {
    if ( link[l] >= 0 )  continue;

    /* current element */
    k = (l-1) / 4 + 1;
    i = (l-1) % 4;
    i1 = idir[i][0];
    i2 = idir[i][1];
    i3 = idir[i][2];
    pt = &mesh->tetra[k];
    mins = MG_MIN(pt->v[i1],MG_MIN(pt->v[i2],pt->v[i3]));
    maxs = MG_MAX(pt->v[i1],MG_MAX(pt->v[i2],pt->v[i3]));
    sum  = pt->v[i1] + pt->v[i2] + pt->v[i3];

    /* accross link */
    ll      = -link[l];
    pp      = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 4 + 1;
      ii = (ll-1) % 4;
      i1 = idir[ii][0];
      i2 = idir[ii][1];
      i3 = idir[ii][2];
      pt1  = &mesh->tetra[kk];
      sum1 = pt1->v[i1] + pt1->v[i2] + pt1->v[i3];
      if ( sum1 == sum ) {
        mins1 = MG_MIN(pt1->v[i1],MG_MIN(pt1->v[i2],pt1->v[i3]));
        maxs1 = MG_MAX(pt1->v[i1],MG_MAX(pt1->v[i2],pt1->v[i3]));

        /* adjacent found */
        if ( mins1 == mins && maxs1 == maxs ) {
          if ( pp != 0 )  link[pp] = link[ll];
          link[l]  = 4*kk + ii;
          link[ll] = 4*k + i;
          break;
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  SAFE_FREE(hcode);
  return(1);
}

/** Create surface adjacency */
int hashTria(pMesh mesh) {
  pTria     pt,pt1;
  Hash      hash;
  hedge    *ph;
  int      *adja,k,jel,lel,hmax,dup,nmf,ia,ib;
  char      i,i1,i2,j,l,ok;
  unsigned int key;

  ADD_MEM(mesh,(3*mesh->nt+4)*sizeof(int),"surfacic adjacency table",return(0));
  SAFE_CALLOC(mesh->adjt,3*mesh->nt+4,int);

  /* adjust hash table params */
  hmax = 3.71*mesh->np;
  hash.siz  = mesh->np;
  hash.max  = hmax + 1;
  hash.nxt  = hash.siz;
  ADD_MEM(mesh,(hash.max+1)*sizeof(hedge),"hash table",return(0));
  SAFE_CALLOC(hash.item,hash.max+1,hedge);

  for (k=hash.siz; k<hash.max; k++)
    hash.item[k].nxt = k+1;

  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");

  /* hash triangles */
  mesh->base = 1;
  dup = nmf = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    pt->flag = 0;
    pt->base = mesh->base;
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      i1 = inxt2[i];
      i2 = iprv2[i];

      /* compute key */
      ia  = MG_MIN(pt->v[i1],pt->v[i2]);
      ib  = MG_MAX(pt->v[i1],pt->v[i2]);
      key = (KA*ia + KB*ib) % hash.siz;
      ph  = &hash.item[key];

      /* store edge */
      if ( ph->a == 0 ) {
        ph->a = ia;
        ph->b = ib;
        ph->k = 3*k + i;
        ph->nxt = 0;
        continue;
      }
      /* update info about adjacent */
      ok = 0;
      while ( ph->a ) {
        if ( ph->a == ia && ph->b == ib ) {
          jel = ph->k / 3;
          j   = ph->k % 3;
          pt1 = &mesh->tria[jel];
          /* discard duplicate face */
          if ( pt1->v[j] == pt->v[i] ) {
            pt1->v[0] = 0;
            dup++;
          }
          /* update adjacent */
          else if ( !mesh->adjt[3*(jel-1)+1+j] ) {
            adja[i] = 3*jel + j;
            mesh->adjt[3*(jel-1)+1+j] = 3*k + i;
          }
          /* non-manifold case */
          else if ( adja[i] != 3*jel+j ) {
            if ( pt->ref == MG_ISO ) {
              lel = mesh->adjt[3*(jel-1)+1+j]/3;
              l   = mesh->adjt[3*(jel-1)+1+j]%3;
              mesh->adjt[3*(lel-1)+1+l] = 0;
              adja[i] = 3*jel+j;
              mesh->adjt[3*(jel-1)+1+j] = 3*k + i;
              (mesh->tria[lel]).tag[l] |= MG_GEO + MG_NOM;
            }
            else {
              pt1->tag[j] |= MG_GEO + MG_NOM;
            }
            pt->tag[i] |= MG_GEO + MG_NOM;
            nmf++;
          }
          ok = 1;
          break;
        }
        else if ( !ph->nxt ) {
          ph->nxt = hash.nxt;
          ph = &hash.item[ph->nxt];
          assert(ph);
          hash.nxt = ph->nxt;
          ph->a = ia;
          ph->b = ib;
          ph->k = 3*k + i;
          ph->nxt = 0;
          ok = 1;
          break;
        }
        else
          ph = &hash.item[ph->nxt];
      }
      if ( !ok ) {
        ph->a = ia;
        ph->b = ib;
        ph->k = 3*k + i;
        ph->nxt = 0;
      }
    }
  }
  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));

  /* set tag */
  for (k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_NOM ) {
        mesh->point[pt->v[inxt2[i]]].tag |= MG_NOM;
        mesh->point[pt->v[iprv2[i]]].tag |= MG_NOM;
      }
    }
  }

  if ( abs(mesh->info.imprim) > 4 && dup+nmf > 0 ) {
    fprintf(stdout,"  ## ");  fflush(stdout);
    if ( nmf > 0 )  fprintf(stdout,"[non-manifold model]  ");
    if ( dup > 0 )  fprintf(stdout," %d duplicate removed",dup);
    fprintf(stdout,"\n");
  }
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- completed.\n");
  return(1);
}

int hashEdge(pMesh mesh,Hash *hash, int a,int b,int k) {
  hedge  *ph;
  int     key,ia,ib,j;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a == ia && ph->b == ib )
    return(1);
  else if ( ph->a ) {
    while ( ph->nxt && ph->nxt < hash->max ) {
      ph = &hash->item[ph->nxt];
      if ( ph->a == ia && ph->b == ib )  return(1);
    }
    ph->nxt   = hash->nxt;
    ph        = &hash->item[hash->nxt];
    ph->a     = ia;  ph->b   = ib;
    ph->k     = k;
    hash->nxt = ph->nxt;
    ph->nxt   = 0;
    if ( hash->nxt >= hash->max ) {
      if ( mesh->info.ddebug )
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",hash->max);
      TAB_RECALLOC(mesh,hash->item,hash->max,0.2,hedge,"edge",return(0));
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    return(1);
  }
  /* insert new edge */
  ph->a = ia;  ph->b = ib;
  ph->k = k;
  ph->nxt = 0;
  return(1);
}

/** return index of point stored along (ia,ib) */
int hashGet(Hash *hash,int a,int b) {
  hedge  *ph;
  int     key,ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->a )  return(0);
  if ( ph->a == ia && ph->b == ib )  return(ph->k);
  while ( ph->nxt ) {
    ph = &hash->item[ph->nxt];
    if ( ph->a == ia && ph->b == ib )  return(ph->k);
  }
  return(0);
}

/** remove edge from hash table */
int hashPop(Hash *hash,int a,int b) {
  hedge  *ph,*php;
  int     key,ia,ib,iph,iphp;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->a ) return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    if ( !ph->nxt ) {
      memset(ph,0,sizeof(hedge));
      return(1);
    }
    else {
      iph = ph->nxt;
      php = ph;
      ph  = &hash->item[ph->nxt];
      memcpy(php,ph,sizeof(hedge));
      memset(ph,0,sizeof(hedge));
      ph->nxt   = hash->nxt;
      hash->nxt = iph;
      return(1);
    }
  }
  while ( ph->nxt ) {
    php = ph;
    ph  = &hash->item[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      if ( !ph->nxt ) {
        memset(ph,0,sizeof(hedge));
        ph->nxt   = hash->nxt;
        hash->nxt = php->nxt;
        php->nxt  = 0;
      }
      else {
        iph  = ph->nxt;
        iphp = php->nxt;
        php->nxt = iph;
        memset(ph,0,sizeof(hedge));
        ph->nxt   = hash->nxt;
        hash->nxt = iphp;
      }
      return(1);
    }
  }
  return(0);
}

/** used to hash edges or faces */
int hashNew(pMesh mesh,Hash *hash,int hsiz,int hmax) {
  int   k;

  /* adjust hash table params */
  ADD_MEM(mesh,(hmax+2)*sizeof(hedge),"hash table",
          return(0));
  SAFE_CALLOC(hash->item,hmax+2,hedge);

  hash->siz  = hsiz;
  hash->max  = hmax + 1;
  hash->nxt  = hsiz;
  for (k=hsiz; k<hash->max; k++)
    hash->item[k].nxt = k+1;

  return(1);
}

/** set tag to edge on geometry */
int hTag(HGeom *hash,int a,int b,int ref,char tag) {
  hgeom  *ph;
  int     key,ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )
    return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    ph->tag |= tag;
    ph->ref  = ref;
    return(1);
  }
  while ( ph->nxt ) {
    ph = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      ph->tag |= tag;
      ph->ref  = ref;
      return(1);
    }
  }
  return(0);
}

/** remove edge from hash table */
int hPop(HGeom *hash,int a,int b,int *ref,char *tag) {
  hgeom  *ph,*php;
  int     key,ia,ib,iph,iphp;

  *ref = 0;
  *tag = 0;
  if ( !hash->siz )  return(0);

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )  return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    *ref = ph->ref;
    *tag = ph->tag;
    if ( !ph->nxt ) {
      memset(ph,0,sizeof(hgeom));
    }
    else {
      iph = ph->nxt;
      php = ph;
      ph  = &hash->geom[ph->nxt];
      memcpy(php,ph,sizeof(hgeom));
      memset(ph,0,sizeof(hgeom));
      ph->nxt   = hash->nxt;
      hash->nxt = iph;
    }
    return(1);
  }
  while ( ph->nxt ) {
    php = ph;
    ph  = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      *ref = ph->ref;
      *tag = ph->tag;
      if ( !ph->nxt ) {
        memset(ph,0,sizeof(hgeom));
        ph->nxt   = hash->nxt;
        hash->nxt = php->nxt;
        php->nxt  = 0;
      }
      else {
        iph  = ph->nxt;
        iphp = php->nxt;
        php->nxt = iph;
        memset(ph,0,sizeof(hgeom));
        ph->nxt   = hash->nxt;
        hash->nxt = iphp;
      }
      return(1);
    }
  }
  return(0);
}

/** get ref and tag to edge on geometry */
int hGet(HGeom *hash,int a,int b,int *ref,char *tag) {
  hgeom  *ph;
  int     key,ia,ib;

  *tag = 0;
  *ref = 0;
  if ( !hash->siz )  return(0);
  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )  return(0);
  else if ( ph->a == ia && ph->b == ib ) {
    *ref = ph->ref;
    *tag = ph->tag;
    return(1);
  }
  while ( ph->nxt ) {
    ph = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      *ref = ph->ref;
      *tag = ph->tag;
      return(1);
    }
  }
  return(0);
}

/** store edge on geometry */
void hEdge(pMesh mesh,int a,int b,int ref,char tag) {
  hgeom  *ph;
  int     key,ia,ib,j;

  if ( !mesh->htab.siz )  return;
  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % mesh->htab.siz;
  ph  = &mesh->htab.geom[key];

  if ( ph->a == ia && ph->b == ib )
    return;
  else if ( ph->a ) {
    while ( ph->nxt ) {
      ph = &mesh->htab.geom[ph->nxt];
      if ( ph->a == ia && ph->b == ib )  return;
    }
    ph->nxt = mesh->htab.nxt;
    ph      = &mesh->htab.geom[mesh->htab.nxt];
    ph->a   = ia;   ph->b   = ib;
    ph->ref = ref;  ph->tag = tag;
    mesh->htab.nxt = ph->nxt;
    ph->nxt = 0;
    if ( mesh->htab.nxt >= mesh->htab.max ) {
      if ( mesh->info.ddebug )
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",mesh->htab.max);
      TAB_RECALLOC(mesh,mesh->htab.geom,mesh->htab.max,0.2,hgeom,
                   "larger htab table",
                   printf("  Exit program.\n");
                   exit(EXIT_FAILURE));
      for (j=mesh->htab.nxt; j<mesh->htab.max; j++)  mesh->htab.geom[j].nxt = j+1;
    }
    return;
  }
  /* insert new edge */
  ph->a   = ia;   ph->b   = ib;
  ph->ref = ref;  ph->tag = tag;
  ph->nxt = 0;
  return;
}

/** to store edge on geometry */
int hNew(HGeom *hash,int hsiz,int hmax,int secure) {
  int   k;

  /* adjust hash table params */
  hash->geom = (hgeom*)calloc(hmax+2,sizeof(hgeom));
  if ( !hash->geom ) {
    perror("  ## Memory problem: calloc");
    if ( !secure )  return(0);
    else  exit(EXIT_FAILURE);
  }
  hash->siz  = hsiz;
  hash->max  = hmax + 1;
  hash->nxt  = hsiz;
  for (k=hsiz; k<hash->max; k++)
    hash->geom[k].nxt = k+1;
  return 1;
}

int hGeom(pMesh mesh) {
  pTria   pt;
  pEdge   pa;
  int    *adja,k,kk,edg;
  char    i,i1,i2,tag;

  /* if edges exist in mesh, hash special edges from existing field */
  if ( mesh->na ) {
    if ( !mesh->htab.geom ) {
      mesh->namax = MG_MAX(1.5*mesh->na,NAMAX);
      ADD_MEM(mesh,(3*mesh->namax+2)*sizeof(hgeom),"htab",return(0));
      hNew(&mesh->htab,mesh->na,3*mesh->namax,1);
    }
    else {
#ifdef SINGUL
      if ( !mesh->info.sing ) {
#endif
        if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
          fprintf(stdout,"  ## Warning: no re-hash of edges of mesh. ");
          fprintf(stdout,"mesh->htab.geom must be freed to enforce analysis.\n");
        }
        DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(Edge));
        mesh->na   = 0;
        return(1);
#ifdef SINGUL
      }
#endif
    }

    /* store initial edges */
    for (k=1; k<=mesh->na; k++) {
      pa = &mesh->edge[k];
      hEdge(mesh,pa->a,pa->b,pa->ref,pa->tag);
    }

    /* now check triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        i1 = inxt2[i];
        i2 = iprv2[i];
        /* transfer non manifold tag to edges */
        if ( pt->tag[i] & MG_NOM )
          hTag(&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);

        hGet(&mesh->htab,pt->v[i1],pt->v[i2],&edg,&tag);
        pt->edg[i]  = edg;
        /* Mark edges as boundary edges */
        pt->tag[i] |= (tag | MG_BDY);
        hTag(&mesh->htab,pt->v[i1],pt->v[i2],edg,pt->tag[i]);
      }
    }
    DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(Edge));
    mesh->na   = 0;
  }
  /* else, infer special edges from information carried by triangles */
  else {
    if ( !mesh->adjt && !hashTria(mesh) )  return(0);
#ifdef SINGUL
    if ( !mesh->info.sing || !mesh->htab.geom ) {
      for (k=1; k<=mesh->nt; k++) {
        pt   = &mesh->tria[k];
        adja = &mesh->adjt[3*(k-1)+1];
        for (i=0; i<3; i++) {
          i1  = inxt2[i];
          i2  = iprv2[i];
          kk  = adja[i] / 3;
          if ( !kk || pt->tag[i] & MG_NOM )
            mesh->na++;
          else if ( (k < kk) && ( pt->edg[i] || pt->tag[i] ) )  mesh->na++;
        }
      }

      if ( mesh->htab.geom )
        DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(hgeom));

      mesh->namax = MG_MAX(1.5*mesh->na,NAMAX);
      ADD_MEM(mesh,(3*mesh->namax+2)*sizeof(hgeom),"htab",return(0));
      hNew(&mesh->htab,mesh->na,3*mesh->namax,1);
      mesh->na = 0;
    }
#else
    for (k=1; k<=mesh->nt; k++) {
      pt   = &mesh->tria[k];
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = inxt2[i];
        i2  = iprv2[i];
        kk  = adja[i] / 3;
        if ( !kk || pt->tag[i] & MG_NOM )
          mesh->na++;
        else if ( (k < kk) && ( pt->edg[i] || pt->tag[i] ) )  mesh->na++;
      }
    }

    if ( mesh->htab.geom )
      DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(hgeom));

    mesh->namax = MG_MAX(1.5*mesh->na,NAMAX);
    ADD_MEM(mesh,(3*mesh->namax+2)*sizeof(hgeom),"htab",return(0));
    hNew(&mesh->htab,mesh->na,3*mesh->namax,1);
    mesh->na = 0;
#endif

    /* build hash for edges */
    for (k=1; k<=mesh->nt; k++) {
      pt   = &mesh->tria[k];
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = inxt2[i];
        i2  = iprv2[i];
        kk  = adja[i] / 3;
        if ( !kk || pt->tag[i] & MG_NOM ) {
          if ( pt->tag[i] & MG_NOM )
            pt->edg[i] = ( mesh->info.iso && pt->edg[i] != 0 ) ?  -abs(pt->edg[i]) : MG_ISO;
          hEdge(mesh,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
        }
        else if ( k < kk && ( pt->edg[i] || pt->tag[i] ) )
          hEdge(mesh,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
      }
    }
    /* now check triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        i1 = inxt2[i];
        i2 = iprv2[i];
        hGet(&mesh->htab,pt->v[i1],pt->v[i2],&edg,&tag);
        pt->edg[i]  = edg;
        pt->tag[i] |= tag;
      }
    }
  }
  return(1);
}

/** Count the number of faces in mesh and compare this number    *
 *  to the number of given triangles. Delete the given triangles *
 *  if they didn't correspond to the founded faces and fill      *
 *  mesh->nt with the real number of faces.                      *
 *  Retrun 1 if numbers of given and founded faces correspond.   */
int chkNumberOfTri(pMesh mesh) {
  pTetra    pt,pt1;
  int      *adja,adj,k,i,nttmp;

  nttmp = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( !adj || ( pt->ref > pt1->ref) )
        nttmp++;
    }
  }
  if ( mesh->nt == nttmp ) return(1);
  else if ( mesh->nt ){
    if ( !mesh->info.iso && (mesh->info.imprim > 4 || mesh->info.ddebug) ) {
      fprintf(stdout,"  ## WARNING: INITIAL TRIANGLES ARE DELETED.\n");
      fprintf(stdout,"  Not enough or too much triangles for geometry (maybe");
      fprintf(stdout," you have 2 domains but only boundary/interface triangles).\n");
      fprintf(stdout," %d given triangles and %d counted triangles.\n",mesh->nt,nttmp);
    }
    DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(Tria));
  }
  mesh->nt = nttmp;
  return(0);
}

/** identify boundary triangles */
int bdryTria(pMesh mesh) {
  pTetra    pt,pt1;
  pTria     ptt;
  pPoint    ppt;
  pxTetra   pxt;
  int      *adja,adj,k;
  char      i;

  /* create triangles */
  ADD_MEM(mesh,(mesh->nt+1)*sizeof(Tria),"triangles",return(0));
  SAFE_CALLOC(mesh->tria,mesh->nt+1,Tria);

  mesh->nt = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    pxt = 0;
    if ( pt->xt )  pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( adj && ( pt->ref <= pt1->ref) )  continue;
      mesh->nt++;
      ptt = &mesh->tria[mesh->nt];
      ptt->v[0] = pt->v[idir[i][0]];
      ptt->v[1] = pt->v[idir[i][1]];
      ptt->v[2] = pt->v[idir[i][2]];
      if ( !adj ) {
        if ( pxt ) {
          if ( pxt->tag[iarf[i][0]] )  ptt->tag[0] = pxt->tag[iarf[i][0]];
          if ( pxt->tag[iarf[i][1]] )  ptt->tag[1] = pxt->tag[iarf[i][1]];
          if ( pxt->tag[iarf[i][2]] )  ptt->tag[2] = pxt->tag[iarf[i][2]];
          /* useful only when saving mesh */
          ptt->ref = pxt->ref[i];
        }
      }
      else {
        if ( pxt ) {
          if ( pxt->tag[iarf[i][0]] )  ptt->tag[0] = pxt->tag[iarf[i][0]];
          if ( pxt->tag[iarf[i][1]] )  ptt->tag[1] = pxt->tag[iarf[i][1]];
          if ( pxt->tag[iarf[i][2]] )  ptt->tag[2] = pxt->tag[iarf[i][2]];
          /* useful only when saving mesh */
        }
        ptt->ref = mesh->info.iso ? MG_ISO : 0;
      }
    }
  }

  /* set point tag */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      ppt->tag |= MG_BDY;
    }
  }

  return(1);
}

/** identify boundary triangles for implicit surface */
/* On ne passe jamais ici non? (pour y passer, il faudrait que l'on ait lu les
   triangles de peau et les triangles ISO or on zappe ces derniers dans loadMesh) */
int bdryIso(pMesh mesh) {
  pTetra    pt,pt1;
  pTria     ptt;
  pPoint    ppt;
  int      *adja,adj,k,nt;
  char      i;

  /* step 1: count interface faces */
  nt = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      if ( adj ) {
        pt1 = &mesh->tetra[adj];
        if ( pt->ref > pt1->ref )  nt++;
      }
    }
  }
  if ( !nt )  return(1);

  /* step 2 : create triangles */
  ADD_MEM(mesh,nt*sizeof(Tria),"triangles",return(0));
  if ( mesh->ntmax < (mesh->nt+nt ) )  mesh->ntmax = mesh->nt+nt;
  SAFE_RECALLOC(mesh->tria,mesh->nt+1,(mesh->nt+nt+1),Tria);

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      if ( adj ) {
        pt1 = &mesh->tetra[adj];
        if ( pt->ref > pt1->ref ) {
          mesh->nt++;
          ptt = &mesh->tria[mesh->nt];
          ptt->v[0] = pt->v[idir[i][0]];
          ptt->v[1] = pt->v[idir[i][1]];
          ptt->v[2] = pt->v[idir[i][2]];
          ptt->ref  = mesh->info.iso ? 100 : 0;  /* useful only when saving mesh */
        }
      }
    }
  }

  /* set point tag */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      ppt->tag |= MG_BDY;
    }
  }

  return(1);
}

static int hashFace(pMesh mesh,Hash *hash,int ia,int ib,int ic,int k) {
  hedge     *ph;
  int        key,mins,maxs,sum,j;

  mins = MG_MIN(ia,MG_MIN(ib,ic));
  maxs = MG_MAX(ia,MG_MAX(ib,ic));

  /* compute key */
  sum = ia + ib + ic;
  key = (KA*mins + KB*maxs) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == mins && ph->b == maxs && ph->s == sum )
      return(ph->k);
    else {
      while ( ph->nxt && ph->nxt < hash->max ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return(ph->k);
      }
    }
    ph->nxt = hash->nxt;
    ph      = &hash->item[hash->nxt];
    ph->a   = mins;
    ph->b   = maxs;
    ph->s   = sum;
    ph->k   = k;
    hash->nxt = ph->nxt;
    ph->nxt = 0;

    if ( hash->nxt >= hash->max ) {
      TAB_RECALLOC(mesh,hash->item,hash->max,0.2,hedge,"face",return(0));
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    return(1);
  }

  /* insert new face */
  ph->a = mins;
  ph->b = maxs;
  ph->s = sum;
  ph->k = k;
  ph->nxt = 0;

  return(1);
}

/** return index of triangle ia ib ic */
static int hashGetFace(Hash *hash,int ia,int ib,int ic) {
  hedge  *ph;
  int     key,mins,maxs,sum;

  if ( !hash->item )  return(0);
  mins = MG_MIN(ia,MG_MIN(ib,ic));
  maxs = MG_MAX(ia,MG_MAX(ib,ic));

  /* compute key */
  sum = ia + ib + ic;
  key = (KA*mins + KB*maxs) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == mins && ph->b == maxs && ph->s == sum )
      return(ph->k);
    else {
      while ( ph->nxt ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return(ph->k);
      }
    }
  }

  return(0);
}

/** if triangles, set ref to tetra faces and edges */
int bdrySet(pMesh mesh) {
  pTetra   pt,pt1;
  pTria    ptt;
  pxTetra  pxt;
  Hash     hash;
  int      *adja,adj,k,kt,ia,ib,ic,j,na;
  char     i,tag;
#ifdef SINGUL
  hgeom    *ph;
  int      ref;

  if ( (!mesh->info.sing) && (!mesh->nt) )  return(1);
#else
  if ( !mesh->nt )  return(1);
#endif

  if ( mesh->xtetra ) {
    if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: no re-build of boundary tetras. ");
      fprintf(stdout,"mesh->xtetra must be freed to enforce analysis.\n");
    }
    return(1);
  }

  if ( ! hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return(0);
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k);
  }
  na = 0;
#ifdef SINGUL
  if ( mesh->info.sing ) {
    for (k=0; k<=mesh->htab.max; k++) {
      ph = &mesh->htab.geom[k];
      if ( !ph->a )  continue;
      na++;
    }
  }
  if ( !mesh->nt ) {
    if ( !na )  return(1);
    else {
      fprintf(stdout,"%s:%d: Error: we should not pass in",__FILE__,__LINE__);
      fprintf(stdout," this routine without triangles and with");
      fprintf(stdout," singular edges.\n");
      return(0);
    }
  }
#endif

  mesh->xt     = 0;
  mesh->xtmax  = mesh->ntmax + 2*na;

  ADD_MEM(mesh,(mesh->xtmax+1)*sizeof(xTetra),"boundary tetrahedra",
          printf("  Exit program.\n");
          exit(EXIT_FAILURE));
  SAFE_CALLOC(mesh->xtetra,mesh->xtmax+1,xTetra);

  /* assign references to tetras faces */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( !adj || (adj > 0 && pt->ref != pt1->ref) ) {
        ia = pt->v[idir[i][0]];
        ib = pt->v[idir[i][1]];
        ic = pt->v[idir[i][2]];
        kt = hashGetFace(&hash,ia,ib,ic);
        assert(kt);
        if ( !pt->xt ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,xTetra,
                         "larger xtetra table",
                         mesh->xt--;
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
          }
          pt->xt = mesh->xt;
        }
        ptt = &mesh->tria[kt];
        pxt = &mesh->xtetra[mesh->xt];
        pxt->ref[i]   = ptt->ref;
        pxt->ftag[i] |= MG_BDY;
        pxt->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);
#ifdef SINGUL
        if ( mesh->info.sing ) {
          pxt->tag[iarf[i][0]] |= MG_BDY;
          pxt->tag[iarf[i][1]] |= MG_BDY;
          pxt->tag[iarf[i][2]] |= MG_BDY;
        }
#endif
      }
    }
  }

#ifdef SINGUL
  /* Add xtetras for singularities */
  if ( mesh->info.sing ) {
    for (k=0; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      for (i=0; i<6; i++) {
        ia = iare[i][0];
        ib = iare[i][1];
        hGet(&mesh->htab,pt->v[ia],pt->v[ib],&ref,&tag);
        if ( pt->xt )
          tag |= mesh->xtetra[pt->xt].tag[i];
        if ( !(tag & MG_BDY) ) {
          if ( tag || ref ) {
            if ( !pt->xt ) {
              mesh->xt++;
              if ( mesh->xt > mesh->xtmax ) {
                TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             printf("  Exit program.\n");
                             exit(EXIT_FAILURE));
              }
              pt->xt = mesh->xt;
            }
            pxt = &mesh->xtetra[pt->xt];
            pxt->edg[i]  = ref;
            pxt->tag[i] |= tag | MG_SGL;
            mesh->point[pt->v[ia]].tag |= MG_SGL;
            mesh->point[pt->v[ib]].tag |= MG_SGL;
          }
        }
      }
    }
  }
#endif

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt )  continue;
    pxt = &mesh->xtetra[pt->xt];
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      /* Set flag to know if tetra has the same orientation than the triangle */
      if ( adj && pt->ref < pt1->ref )  MG_CLR(pxt->ori,i);
      else  MG_SET(pxt->ori,i);
      /* Set edge tag */
      if ( pxt->ftag[i] ) {
        if ( adj && (pt->ref <= pt1->ref || (pt->ref == MG_PLUS)) ) {
          continue;
        }
        else {
          ia = pt->v[idir[i][0]];
          ib = pt->v[idir[i][1]];
          ic = pt->v[idir[i][2]];
          kt = hashGetFace(&hash,ia,ib,ic);
          ptt = &mesh->tria[kt];
          for (j=0; j<3; j++) {
            tag = pxt->ftag[i] | ptt->tag[j];
            if ( tag )
              settag(mesh,k,iarf[i][j],tag,ptt->edg[j]);
          }
        }
      }
    }
  }
  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
  return(1);
}

/** Update tag and refs of tetra edges.
    If tetra is required, set the faces/edges to required */
int bdryUpdate(pMesh mesh) {
  pTetra   pt;
  pTria    ptt;
  pxTetra  pxt;
  Hash     hash;
  int      k,kt,ia,ib,ic,j;
  char     i,tag;

  if ( !mesh->nt )  return(1);
  if ( !hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) )  return(0);
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k);
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( pt->tag & MG_REQ ) {
      mesh->point[mesh->tetra[k].v[0]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[1]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[2]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[3]].tag |= MG_REQ;
      settag(mesh,k,0,MG_REQ,0);
      settag(mesh,k,1,MG_REQ,0);
      settag(mesh,k,2,MG_REQ,0);
      settag(mesh,k,3,MG_REQ,0);
      settag(mesh,k,4,MG_REQ,0);
      settag(mesh,k,5,MG_REQ,0);
    }

    if ( !pt->xt )  continue;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++) {
      /* Set edge tag */
      if ( ! MG_GET(pxt->ori,i) ) continue;
      if ( pxt->ftag[i] & MG_BDY ) {
        ia = pt->v[idir[i][0]];
        ib = pt->v[idir[i][1]];
        ic = pt->v[idir[i][2]];
        kt = hashGetFace(&hash,ia,ib,ic);
        assert(kt);
        ptt = &mesh->tria[kt];
        if ( pt->tag & MG_REQ ) {
          pxt->ftag[i] |= MG_REQ;
          ptt->tag[0]   = MG_REQ;
          ptt->tag[1]   = MG_REQ;
          ptt->tag[2]   = MG_REQ;
        }
        for ( j=0; j<3; j++ ) {
          tag = ptt->tag[j];
          if ( tag || ptt->edg[j] )
            settag(mesh,k,iarf[i][j],tag,ptt->edg[j]);
        }
      }
    }
  }
  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
  return(1);
}

/** make orientation of triangles compatible with tetra faces */
int bdryPerm(pMesh mesh) {
  pTetra   pt,pt1;
  pTria    ptt;
  pPoint   ppt;
  Hash     hash;
  int     *adja,adj,k,kt,ia,ib,ic,nf;
  char     i;

  assert(mesh->nt);

  /* store triangles temporarily */
  if ( !hashNew(mesh,&hash,MG_MAX(0.51*mesh->nt,100),MG_MAX(1.51*mesh->nt,300)) )
    return(0);
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) )  return(0);
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      if ( !mesh->info.iso ) ppt->tag |= MG_BDY;
    }
  }

  /* check orientation */
  nf = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( adj && (pt->ref <= pt1->ref || pt->ref == MG_PLUS) )
        continue;
      ia = pt->v[idir[i][0]];
      ib = pt->v[idir[i][1]];
      ic = pt->v[idir[i][2]];
      kt = hashGetFace(&hash,ia,ib,ic);
      if ( !kt ) {
        fprintf(stdout,"%s:%d: Error: function hashGetFace return 0.\n",__FILE__,__LINE__);
        fprintf(stdout," Maybe you have non-boundary triangles.");
        fprintf(stdout," Check triangle of vertices %d %d %d.\n",ia,ib,ic);
        exit(EXIT_FAILURE);
      }

      /* check orientation */
      ptt = &mesh->tria[kt];
      if ( ptt->v[0] == ia && ptt->v[1] == ib && ptt->v[2] == ic )
        continue;
      else {
        ptt->v[0] = ia;
        ptt->v[1] = ib;
        ptt->v[2] = ic;
        nf++;
      }
    }
  }
  if ( mesh->info.ddebug && nf > 0 )
    fprintf(stdout,"  ## %d faces reoriented\n",nf);

  DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(hedge));
  return(1);
}
