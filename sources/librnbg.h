/**
 * Logiciel initial: MMG3D Version 4.0
 * Co-auteurs : Cecile Dobrzynski et Pascal Frey.
 * Proprietaires :IPB - UPMC -INRIA.
 *
 * Copyright \copyright 2004-2005-2006-2007-2008-2009-2010-2011,
 * diffuse sous les termes et conditions de la licence publique generale de GNU
 * Version 3 ou toute version ulterieure.
 *
 * Ce fichier est une partie de MMG3D.
 * MMG3D est un logiciel libre ; vous pouvez le redistribuer et/ou le modifier
 * suivant les termes de la licence publique generale de GNU
 * Version 3 ou toute version ulterieure.
 * MMG3D est distribue dans l espoir qu il sera utile, mais SANS
 * AUCUNE GARANTIE ; sans meme garantie de valeur marchande.
 * Voir la licence publique generale de GNU pour plus de details.
 * MMG3D est diffuse en esperant qu il sera utile,
 * mais SANS AUCUNE GARANTIE, ni explicite ni implicite,
 * y compris les garanties de commercialisation ou
 * d adaptation dans un but specifique.
 * Reportez-vous a la licence publique generale de GNU pour plus de details.
 * Vous devez avoir re\c{c}u une copie de la licence publique generale de GNU
 * en meme temps que ce document.
 * Si ce n est pas le cas, aller voir <http://www.gnu.org/licenses/>.**/
/**
 * Initial software: MMG3D Version 4.0
 * Co-authors: Cecile Dobrzynski et Pascal Frey.
 * Owners: IPB - UPMC -INRIA.
 *
 * Copyright \copyright 2004-2005-2006-2007-2008-2009-2010-2011,
 * spread under the terms and conditions of the license GNU General Public License
 * as published Version 3, or (at your option) any later version.
 *
 * This file is part of MMG3D
 * MMG3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * MMG3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with MMG3D. If not, see <http://www.gnu.org/licenses/>.
 **/
/** librnbg
 *
 * Written by Cedric Lachat and Algiane Froehly
 **/
#ifdef USE_SCOTCH

#ifndef __RENUM__
#define __RENUM__

#include <scotch.h>

#define HASHPRIME 37

typedef struct MeshGraphHash_ {
  int vertNum;
  int vertEnd;
} MeshGraphHash;

int _SCOTCHintSort2asc1(SCOTCH_Num * sortPartTb, int vertNbr);

#endif /* __RENUM__ */
#endif
