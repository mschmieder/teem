/*
  teem: Gordon Kindlmann's research software
  Copyright (C) 2003, 2002, 2001, 2000, 1999, 1998 University of Utah

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "limn.h"

/*
******** limnHestCamOptAdd()
**
** calls hestOptAdd a bunch of times to set up command-line options
** useful for specifying a limnCam.  The flags used are as follows:
** fr: cam->from
** at: cam->at
** up: cam->up
** rh: cam->rightHanded
** or: cam->ortho
** dn: cam->neer
** di: cam->dist
** df: cam->faar
** ar: cam->atRel
** ur: cam->uRange
** vr: cam->vRange
*/
void
limnHestCamOptAdd(hestOpt **hoptP, limnCam *cam,
		  char *frDef, char *atDef, char *upDef,
		  char *dnDef, char *diDef, char *dfDef,
		  char *urDef, char *vrDef) {
  hestOpt *hopt;
  
  hopt = *hoptP;
  hestOptAdd(&hopt, "fr", "eye pos", airTypeDouble, 3, 3, cam->from,
	     frDef, "camera eye point");
  hestOptAdd(&hopt, "at", "at pos", airTypeDouble, 3, 3, cam->at,
	     atDef, "camera look-at point");
  hestOptAdd(&hopt, "up", "up dir", airTypeDouble, 3, 3, cam->up,
	     upDef, "camera pseudo-up vector");
  hestOptAdd(&hopt, "rh", NULL, airTypeInt, 0, 0, &(cam->rightHanded), NULL,
	     "use a right-handed UVN frame (V points down)");
  hestOptAdd(&hopt, "or", NULL, airTypeInt, 0, 0, &(cam->ortho), NULL,
	     "orthogonal (not perspective) projection");
  hestOptAdd(&hopt, "dn", "near", airTypeDouble, 1, 1, &(cam->neer),
	     dnDef, "distance to near clipping plane");
  hestOptAdd(&hopt, "di", "image", airTypeDouble, 1, 1, &(cam->dist),
	     diDef, "distance to image plane");
  hestOptAdd(&hopt, "df", "far", airTypeDouble, 1, 1, &(cam->faar),
	     dfDef, "distance to far clipping plane");
  hestOptAdd(&hopt, "ar", NULL, airTypeInt, 0, 0, &(cam->atRel), NULL,
	     "near, image, and far plane distances are relative to "
	     "the *at* point, instead of the eye point");
  hestOptAdd(&hopt, "ur", "uMin uMax", airTypeDouble, 2, 2, cam->uRange,
	     urDef, "range in U direction of image plane");
  hestOptAdd(&hopt, "vr", "vMin vMax", airTypeDouble, 2, 2, cam->vRange,
	     vrDef, "range in V direction of image plane");
  *hoptP = hopt;
  return;
}
