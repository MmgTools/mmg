//
//  cenrad.c
//  MMG3D5delone
//
//  Created by Dobrzynski on 08/07/13.
//  Copyright (c) 2013 Dobrzynski. All rights reserved.
//

#include "mmg3d.h"
/* compute circumradius and center */
int cenrad_iso(pMesh mesh,double *ct,double *c,double *rad) {
    double      dd,ux,uy,uz,n1[3],n2[3],n3[3],*c1,*c2,*c3,*c4,pl1,pl2,pl3;
    double      cc1,cc2,cc3;
    
    c1 = &ct[0];
    c2 = &ct[3];
    c3 = &ct[6];
    c4 = &ct[9];
    
    ux = c4[0] - c1[0];
    uy = c4[1] - c1[1];
    uz = c4[2] - c1[2];
    //if(fabs(sqrt(ux*ux + uy*uy + uz*uz))<1e-12) printf("garg 1 %e\n",fabs(sqrt(ux*ux + uy*uy + uz*uz)));
    dd = 1.0 / sqrt(ux*ux + uy*uy + uz*uz);
    n1[0] = ux*dd;
    n1[1] = uy*dd;
    n1[2] = uz*dd;
    
    /* plan: vecteur directeur passant par milieu(1,4) */
    pl1 = n1[0]*(c4[0]+c1[0]) \
    + n1[1]*(c4[1]+c1[1]) + n1[2]*(c4[2]+c1[2]);
    
    ux = c4[0] - c2[0];
    uy = c4[1] - c2[1];
    uz = c4[2] - c2[2];
    //if(fabs(sqrt(ux*ux + uy*uy + uz*uz))<1e-12) printf("garg 2 %e\n",fabs(sqrt(ux*ux + uy*uy + uz*uz)));
    dd = 1.0 / sqrt(ux*ux + uy*uy + uz*uz);
    n2[0] = ux*dd;
    n2[1] = uy*dd;
    n2[2] = uz*dd;
    pl2 = n2[0]*(c4[0]+c2[0]) \
    + n2[1]*(c4[1]+c2[1]) + n2[2]*(c4[2]+c2[2]);
    
    ux = c4[0] - c3[0];
    uy = c4[1] - c3[1];
    uz = c4[2] - c3[2];
    //if(fabs(sqrt(ux*ux + uy*uy + uz*uz))<1e-12) printf("garg 3 %e\n",fabs(sqrt(ux*ux + uy*uy + uz*uz)));
    dd = 1.0 / sqrt(ux*ux + uy*uy + uz*uz);
    n3[0] = ux*dd;
    n3[1] = uy*dd;
    n3[2] = uz*dd;
    pl3 = n3[0]*(c4[0]+c3[0]) \
    + n3[1]*(c4[1]+c3[1]) + n3[2]*(c4[2]+c3[2]);
    
    /* center = intersection of 3 planes */
    ux = n2[1]*n3[2] - n2[2]*n3[1];
    uy = n1[2]*n3[1] - n1[1]*n3[2];
    uz = n1[1]*n2[2] - n1[2]*n2[1];
    
    dd = n1[0]*ux + n2[0]*uy + n3[0]*uz;
    if(fabs((dd))<1e-12)  return(0);
    dd = 0.5 / dd;
    
    cc1 = ux*pl1 + uy*pl2 + uz*pl3;
    cc2 = pl1 * (n2[2]*n3[0] - n2[0]*n3[2]) \
    + pl2 * (n1[0]*n3[2] - n3[0]*n1[2]) \
    + pl3 * (n2[0]*n1[2] - n2[2]*n1[0]);
    cc3 = pl1 * (n2[0]*n3[1] - n2[1]*n3[0]) \
    + pl2 * (n3[0]*n1[1] - n3[1]*n1[0]) \
    + pl3 * (n1[0]*n2[1] - n2[0]*n1[1]);
    
    c[0] = dd * cc1;
    c[1] = dd * cc2;
    c[2] = dd * cc3;
    
    /* radius (squared) */
    *rad = (c[0] - c4[0]) * (c[0] - c4[0]) \
    + (c[1] - c4[1]) * (c[1] - c4[1]) \
    + (c[2] - c4[2]) * (c[2] - c4[2]);
    
    return(1);
}

