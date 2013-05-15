#include "mmg3d.h"

#define KA     7
#define KB    11
#define KC    13

extern Info  info;
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

int hashTetra(pMesh mesh) {
  pTetra    pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
  int      *hcode,*link,hsize,inival;
  unsigned char  i,ii,i1,i2,i3;
  unsigned int   key;

  /* default */
  if ( mesh->adja )  return(1);
  if ( abs(info.imprim) > 5 || info.ddebug )
    fprintf(stdout,"  ** SETTING STRUCTURE\n");

  /* packing : if not hash does not work */
  paktet(mesh);

  /* memory alloc */
  mesh->adja = (int*)calloc(4*mesh->nemax+5,sizeof(int));
  assert(mesh->adja);
  hcode = (int*)calloc(mesh->ne+5,sizeof(int));
  assert(hcode);
  link  = mesh->adja;
  hsize = mesh->ne;

  /* init */
  if ( info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");
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
  if ( info.ddebug )  fprintf(stdout,"  h- stage 2: adjacencies\n");
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
  free(hcode);
  hcode=NULL;
  return(1);
}

/** Create surface adjacency */
int hashTria(pMesh mesh) {
  pTria     pt,pt1;
  Hash      hash;
  hedge    *ph;
  int      *adja,k,jel,hmax,dup,nmf,ia,ib;
  char      i,i1,i2,j,ok;
  unsigned int key;

  mesh->adjt = (int*)calloc(3*mesh->nt+5,sizeof(int));
  assert(mesh->adjt);

  /* adjust hash table params */
  hmax = 3.71*mesh->np;
  hash.siz  = mesh->np;
  hash.max  = hmax + 1;
  hash.nxt  = hash.siz;
  hash.item = (hedge*)calloc(hmax+1,sizeof(hedge));
  assert(hash.item);

  for (k=hash.siz; k<hash.max; k++)
    hash.item[k].nxt = k+1;

  if ( info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");

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
            pt->tag[i] |= MG_GEO + MG_NOM;
            pt1->tag[j]|= MG_GEO + MG_NOM;
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
  free(hash.item);
  hash.item = NULL;

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

  if ( abs(info.imprim) > 4 && dup+nmf > 0 ) {
    fprintf(stdout,"  ## ");  fflush(stdout);
    if ( nmf > 0 )  fprintf(stdout,"[non-manifold model]  ");
    if ( dup > 0 )  fprintf(stdout," %d duplicate removed",dup);
    fprintf(stdout,"\n");
  }
  if ( info.ddebug )  fprintf(stdout,"  h- completed.\n");
  return(1);
}

int hashEdge(Hash *hash,int a,int b,int k) {
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
    ph->nxt = hash->nxt;
    ph      = &hash->item[hash->nxt];
    ph->a   = ia;  ph->b   = ib;
    ph->k   = k;
    hash->nxt = ph->nxt;
    ph->nxt = 0;
    if ( hash->nxt >= hash->max ) {
      if ( info.ddebug )  fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",hash->max);
      hash->max *= 1.2;
      hash->item = (hedge*)realloc(hash->item,hash->max*sizeof(hedge));
      assert(hash->item);
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
int hashNew(Hash *hash,int hsiz,int hmax) {
  int   k;

  /* adjust hash table params */
  hash->item = (hedge*)calloc(hmax+10,sizeof(hedge));
  assert(hash->item);
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
void hEdge(HGeom *hash,int a,int b,int ref,char tag) {
  hgeom  *ph;
  int     key,ia,ib,j;

  if ( !hash->siz )  return;
  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (KA*ia + KB*ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( ph->a == ia && ph->b == ib )
    return;
  else if ( ph->a ) {
    while ( ph->nxt ) {
      ph = &hash->geom[ph->nxt];
      if ( ph->a == ia && ph->b == ib )  return;
    }
    ph->nxt = hash->nxt;
    ph      = &hash->geom[hash->nxt];
    ph->a   = ia;   ph->b   = ib;
    ph->ref = ref;  ph->tag = tag;
    hash->nxt = ph->nxt;
    ph->nxt = 0;
    if ( hash->nxt >= hash->max ) {
      if ( info.ddebug )  fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",hash->max);
      hash->max *= 1.2;
      hash->geom = (hgeom*)realloc(hash->geom,hash->max*sizeof(hgeom));
      assert(hash->geom);
      for (j=hash->nxt; j<hash->max; j++)  hash->geom[j].nxt = j+1;
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
void hNew(HGeom *hash,int hsiz,int hmax) {
  int   k;

  /* adjust hash table params */
  hash->geom = (hgeom*)calloc(hmax+10,sizeof(hgeom));
  assert(hash->geom);
  hash->siz  = hsiz;
  hash->max  = hmax + 1;
  hash->nxt  = hsiz;
  for (k=hsiz; k<hash->max; k++)
    hash->geom[k].nxt = k+1;
  return;
}

int hGeom(pMesh mesh) {
  pTria   pt;
  pEdge   pa;
  int    *adja,k,kk,edg;
  char    i,i1,i2,tag;

  /* if edges exist in mesh, hash special edges from existing field */
  if ( mesh->na ) {
    mesh->namax = MG_MAX(1.5*mesh->na,NAMAX);
    hNew(&mesh->htab,mesh->na,3*mesh->namax);

    /* store initial edges */
    for (k=1; k<=mesh->na; k++) {
      pa = &mesh->edge[k];
      hEdge(&mesh->htab,pa->a,pa->b,pa->ref,pa->tag);
    }
    if ( mesh->na )  free(mesh->edge);
    mesh->edge = 0;
    mesh->na = 0;

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
        pt->tag[i] |= tag;
      }
    }
  }
  /* else, infer special edges from information carried by triangles */
  else {
		if ( !mesh->adjt && !hashTria(mesh) )  return(0);
		for (k=1; k<=mesh->nt; k++) {
      pt   = &mesh->tria[k];
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = inxt2[i];
        i2  = iprv2[i];
        kk  = adja[i] / 3;
        if ( !kk || pt->tag[i] & MG_NOM )
          mesh->na++;
        else if ( (k < kk) && pt->edg[i]+pt->tag[i] )  mesh->na++;
      }
    }
    mesh->namax = MG_MAX(1.5*mesh->na,NAMAX);
    if(mesh->htab.geom){
      free(mesh->htab.geom);
      mesh->htab.geom=NULL;
    }
    hNew(&mesh->htab,mesh->na,3*mesh->namax);
    mesh->na = 0;

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
            pt->edg[i] = ( info.iso && pt->edg[i] != 0 ) ?  -abs(pt->edg[i]) : MG_ISO;
          hEdge(&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
        }
        else if ( k < kk && pt->edg[i]+pt->tag[i] )
          hEdge(&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
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

/** identify boundary triangles */
int bdryTria(pMesh mesh) {
  pTetra    pt,pt1;
  pTria     ptt;
  pPoint    ppt;
  pxTetra   pxt;
  int      *adja,adj,k,edg;
  char      i,i1,i2,tag;

  if ( mesh->nt )  return(1);

  /* step 1: count external faces */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( !adj || (k < adj && pt->ref != pt1->ref) )
        mesh->nt++;
    }
  }
  if ( !mesh->nt )  return(1);

  /* step 2 : create triangles */
  mesh->tria = (pTria)calloc(mesh->nt+1,sizeof(Tria));
  assert(mesh->tria);

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
      if ( !adj || (k < adj && pt->ref != pt1->ref) ) {
        mesh->nt++;
        ptt = &mesh->tria[mesh->nt];
        ptt->v[0] = pt->v[idir[i][0]];
        ptt->v[1] = pt->v[idir[i][1]];
        ptt->v[2] = pt->v[idir[i][2]];
        if ( !adj ) {
          if ( pxt ) {
            if ( pxt->tag[iarf[i][0]] & MG_REQ ) ptt->tag[0] = MG_REQ;
            if ( pxt->tag[iarf[i][1]] & MG_REQ ) ptt->tag[1] = MG_REQ;
            if ( pxt->tag[iarf[i][2]] & MG_REQ ) ptt->tag[2] = MG_REQ;
            /* useful only when saving mesh */
            ptt->ref = pxt->ref[i];
          }
        }
        else
					ptt->ref = info.iso ? MG_ISO : 0;
      }
    }
  }

  /* set point tag */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      ppt->tag |= MG_BDY;

      i1 = inxt2[i];
      i2 = iprv2[i];
      hGet(&mesh->htab,ptt->v[i1],ptt->v[i2],&edg,&tag);
      ptt->edg[i]  = edg;
      ptt->tag[i] |= tag;
    }
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    pt->xt = 0;
  }
  free(mesh->xtetra);
  mesh->xtetra = 0;
  mesh->xt = 0;
  return(1);
}

/** identify boundary triangles for implicit surface */
int bdryIso(pMesh mesh) {
  pTetra    pt,pt1;
  pTria     ptt;
  pPoint    ppt;
  int      *adja,adj,k,nt;
  char      i,alloc;

  alloc = mesh->nt == 0;
  /* step 1: count external faces */
  nt = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( k < adj && pt->ref != pt1->ref )  nt++;
    }
  }
  if ( !nt )  return(1);

  /* step 2 : create triangles */
  if ( alloc )
    mesh->tria = (pTria)calloc(mesh->nt+nt+1,sizeof(Tria));
  else if ( mesh->nt+nt >= mesh->ntmax )
    mesh->tria = (pTria)realloc(mesh->tria,(mesh->nt+nt+1)*sizeof(Tria));
  assert(mesh->tria);

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];
      if ( k < adj && pt->ref != pt1->ref ) {
        mesh->nt++;
        ptt = &mesh->tria[mesh->nt];
        ptt->v[0] = pt->v[idir[i][0]];
        ptt->v[1] = pt->v[idir[i][1]];
        ptt->v[2] = pt->v[idir[i][2]];
        ptt->ref  = info.iso ? 100 : 0;  /* useful only when saving mesh */
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

static int hashFace(Hash *hash,int ia,int ib,int ic,int k) {
  hedge     *ph;
  int        key,mins,maxs,sum;

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
    ++hash->nxt;
    if ( hash->nxt == hash->max ) {
      fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",hash->max);
      return(0);
    }
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
  int     *adja,adj,k,kt,ia,ib,ic,j;
  char     i,tag;

  if ( !mesh->nt )  return(1);

  hashNew(&hash,0.51*mesh->nt,1.51*mesh->nt);
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    hashFace(&hash,ptt->v[0],ptt->v[1],ptt->v[2],k);
  }
  mesh->xt     = 0;
  mesh->xtmax  = mesh->ntmax; //10 * NTMAX;
  mesh->xtetra = (pxTetra)calloc(mesh->xtmax+1,sizeof(xTetra));
  assert(mesh->xtetra);

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
          pt->xt = mesh->xt;
        }
        ptt = &mesh->tria[kt];
        pxt = &mesh->xtetra[mesh->xt];
        pxt->ref[i]   = ptt->ref;
        pxt->ftag[i] |= MG_BDY;
        pxt->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);
      }
    }
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt )  continue;
    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] ) {
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
  free(hash.item);
  hash.item=NULL;
  return(1);
}

/** make orientation of triangles compatible with tetra faces */
int bdryPerm(pMesh mesh,int iso) {
  pTetra   pt,pt1;
  pTria    ptt;
  pPoint   ppt;
  Hash     hash;
  int     *adja,adj,k,kt,ia,ib,ic,nf;
  char     i;

  if ( !mesh->nt )  return(0);

  /* store triangles temporarily */
  hashNew(&hash,MG_MAX(0.51*mesh->nt,100),MG_MAX(1.51*mesh->nt,300));
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !hashFace(&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) )  return(0);
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      if (!iso) ppt->tag |= MG_BDY;
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
      if ( adj && (pt->ref == pt1->ref || pt->ref == MG_PLUS) )
        continue;
      else {
        ia = pt->v[idir[i][0]];
        ib = pt->v[idir[i][1]];
        ic = pt->v[idir[i][2]];
        kt = hashGetFace(&hash,ia,ib,ic);
        if ( !kt ) {
	        fprintf(stdout,"%s:%d: Error: function hashGetFace return 0\n",__FILE__,__LINE__);
	        fprintf(stdout," Maybe you have non-boundary triangles.");
	        fprintf(stdout," Check triangle %d %d %d\n",ia,ib,ic);
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
  }
  if ( info.ddebug && nf > 0 )
    fprintf(stdout,"  ## %d faces reoriented\n",nf);

  free(hash.item);
  hash.item=NULL;
  return(1);
}
