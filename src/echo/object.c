/*
  teem: Gordon Kindlmann's research software
  Copyright (C) 2002, 2001, 2000, 1999, 1998 University of Utah

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

#include "echo.h"

#define NEW_TMPL(TYPE, BODY)                                     \
EchoObject##TYPE *                                               \
_echoObject##TYPE##_new(void) {                                  \
  EchoObject##TYPE *obj;                                         \
                                                                 \
  obj = (EchoObject##TYPE *)calloc(1, sizeof(EchoObject##TYPE)); \
  obj->type = echoObject##TYPE;                                  \
  do { BODY } while (0);                                         \
  return obj;                                                    \
}

/*
  echoObjectUnknown,
  echoObjectSphere,     /*  1 */
  echoObjectCube,       /*  2 */
  echoObjectUnitCube,   /*  3 */
  echoObjectTriangle,   /*  4 */
  echoObjectRectangle,  /*  5 */
  echoObjectTriMesh,    /*  6: only triangles in the mesh */
  echoObjectIsosurface, /*  7 */
  echoObjectAABBox,     /*  8 */
  echoObjectList,       /*  9 */
  echoObjectInstance,   /* 10 */
  echoObjectLast
*/

NEW_TMPL(Sphere,);
NEW_TMPL(Cube,);
NEW_TMPL(UnitCube,);
NEW_TMPL(Triangle,);
NEW_TMPL(Rectangle,);
NEW_TMPL(TriMesh, /* ??? */ );
NEW_TMPL(Isosurface, /* ??? */);
NEW_TMPL(AABBox,
	 int i;
	 for (i=0; i<ECHO_AABBOX_OBJECT_MAX; i++) {
	   obj->obj[i] = NULL;
	 }
	 obj->len = 0;
	 );
NEW_TMPL(List,
	 obj->obj = NULL;
	 obj->objArr = airArrayNew((void**)&(obj->obj), NULL,
				   sizeof(EchoObject *),
				   ECHO_LIST_OBJECT_INCR);
	 airArrayPointerCB(obj->objArr,
			   airNull,
			   (void *(*)(void *))echoObjectNix);
	 );
NEW_TMPL(Instance,
	 obj->own = AIR_FALSE;
	 obj->obj = NULL;
	 );

EchoObject *(*
_echoObjectNew[ECHO_OBJECT_MAX+1])(void) = {
  NULL,
  (EchoObject *(*)(void))_echoObjectSphere_new,
  (EchoObject *(*)(void))_echoObjectCube_new,
  (EchoObject *(*)(void))_echoObjectUnitCube_new,
  (EchoObject *(*)(void))_echoObjectTriangle_new,
  (EchoObject *(*)(void))_echoObjectRectangle_new,
  (EchoObject *(*)(void))_echoObjectTriMesh_new,
  (EchoObject *(*)(void))_echoObjectIsosurface_new,
  (EchoObject *(*)(void))_echoObjectAABBox_new,
  (EchoObject *(*)(void))_echoObjectList_new,
  (EchoObject *(*)(void))_echoObjectInstance_new
};

EchoObject *
echoObjectNew(int type) {
  
  return _echoObjectNew[type]();
}

/* ---------------------------------------------------- */

#define NIX_TMPL(TYPE, BODY)                                     \
EchoObject##TYPE *                                               \
_echoObject##TYPE##_nix(EchoObject##TYPE *obj) {                 \
                                                                 \
  do { BODY } while (0);                                         \
  free(obj);                                                     \
  return NULL;                                                   \
}

NIX_TMPL(TriMesh, /* ??? */);
NIX_TMPL(Isosurface, /* ??? */);
NIX_TMPL(AABBox,
	 int i;
	 for (i=0; i<obj->len; i++) {
	   echoObjectNix(obj->obj[i]);
	 }
	 );
NIX_TMPL(List,
	 /* due to airArray callbacks, this will nuke all kids */
	 airArrayNuke(obj->objArr);
	 );
NIX_TMPL(Instance,
	 if (obj->own) {
	   echoObjectNix(obj->obj);
	 }
	 );

EchoObject *(*
_echoObjectNix[ECHO_OBJECT_MAX+1])(EchoObject *) = {
  NULL,
  (EchoObject *(*)(EchoObject *))airFree,
  (EchoObject *(*)(EchoObject *))airFree,
  (EchoObject *(*)(EchoObject *))airFree,
  (EchoObject *(*)(EchoObject *))airFree,
  (EchoObject *(*)(EchoObject *))airFree,
  (EchoObject *(*)(EchoObject *))_echoObjectTriMesh_nix,
  (EchoObject *(*)(EchoObject *))_echoObjectIsosurface_nix,
  (EchoObject *(*)(EchoObject *))_echoObjectAABBox_nix,
  (EchoObject *(*)(EchoObject *))_echoObjectList_nix,
  (EchoObject *(*)(EchoObject *))_echoObjectInstance_nix
};

EchoObject *
echoObjectNix(EchoObject *obj) {

  return _echoObjectNix[obj->type](obj);
}

int
echoObjectIsContainer(EchoObject *obj) {
  
  if (!obj)
    return AIR_FALSE;

  if (echoObjectAABBox == obj->type ||
      echoObjectInstance == obj->type ||
      echoObjectList == obj->type) {
    return AIR_TRUE;
  }
  else {
    return AIR_FALSE;
  }
}

void
echoObjectListAdd(EchoObject *parent, EchoObject *child) {
  EchoObjectList *list;
  int idx;
  
  if (!(parent && child &&
	echoObjectList == parent->type))
    return;

  list = (EchoObjectList *)parent;
  idx = airArrayIncrLen(list->objArr, 1);
  list->obj[idx] = child;

  return;
}

void
echoObjectSphereSet(EchoObject *_sphere,
		    echoPos_t x, echoPos_t y, echoPos_t z, echoPos_t rad) {
  EchoObjectSphere *sphere;

  if (_sphere && echoObjectSphere == _sphere->type) {
    sphere = (EchoObjectSphere *)_sphere;
    ELL_3V_SET(sphere->pos, x, y, z);
    sphere->rad = rad;
  }
  return;
}

void
echoObjectRectangleSet(EchoObject *_rect,
		       echoPos_t ogx, echoPos_t ogy, echoPos_t ogz,
		       echoPos_t e0x, echoPos_t e0y, echoPos_t e0z,
		       echoPos_t e1x, echoPos_t e1y, echoPos_t e1z) {
  EchoObjectRectangle *rect;

  if (_rect && echoObjectRectangle == _rect->type) {
    rect = (EchoObjectRectangle *)_rect;
    ELL_3V_SET(rect->origin, ogx, ogy, ogz);
    ELL_3V_SET(rect->edge0, e0x, e0y, e0z);
    ELL_3V_SET(rect->edge1, e1x, e1y, e1z);
    rect->area = ELL_3V_LEN(rect->edge0) * ELL_3V_LEN(rect->edge1);
  }
  return;
}
		       
void
echoObjectTriangleSet(EchoObject *_tri,
		      echoPos_t x0, echoPos_t y0, echoPos_t z0, 
		      echoPos_t x1, echoPos_t y1, echoPos_t z1, 
		      echoPos_t x2, echoPos_t y2, echoPos_t z2) {
  EchoObjectTriangle *tri;

  if (_tri && echoObjectTriangle == _tri->type) {
    tri = (EchoObjectTriangle *)_tri;
    ELL_3V_SET(tri->vert[0], x0, y0, z0);
    ELL_3V_SET(tri->vert[1], x1, y1, z1);
    ELL_3V_SET(tri->vert[2], x2, y2, z2);
  }
  return;
}
		       
