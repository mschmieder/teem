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

#include "nrrd.h"
#include "privateNrrd.h"

/* int _VV = 0;*/

/*
** learned: even when using doubles, because of limited floating point
** precision, you can get different results between quantizing 
** unrescaled (value directly from nrrd, map domain set to nrrd range,
** as with early behavior of unu rmap) and rescaled (value from nrrd
** scaled to fit in existing map domain, as with unu imap -r) value,
** to the exact same index space.
*/

/*
** I won't try to support mapping individual values through a
** colormap, as with a function evaluation on a single passed value.
** That will be handled in an upcoming library...
*/

/* 
** this identifies the different kinds of 1D maps, useful for the
** functions in this file only
*/
enum {
  kindLut=0,
  kindRmap=1,
  kindImap=2
};

/*
******** nrrd1DIrregMapCheck()
**
** return zero only for the valid forms of 1D irregular map.
** imap must be 2D, both sizes >= 2, non-block-type, no non-existant
** values in range.  If the first point's position is non-existant,
** than the first three points positions must be -inf, NaN, and +inf,
** and none of the other points locations can be non-existant, and
** they must increase monotonically.  There must be at least two
** points with existant positions.
*/
int
nrrd1DIrregMapCheck(Nrrd *nmap) {
  char me[]="nrrd1DIrregMapCheck", err[AIR_STRLEN_MED];
  double (*mapLup)(void *v, size_t I);
  int i, entLen, mapLen, baseI, min[2], max[2];
  Nrrd *nrange;

  if (!nmap) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdCheck(nmap)) {
    sprintf(err, "%s: ", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdTypeBlock == nmap->type) {
    sprintf(err, "%s: map is %s type, need scalar", 
	    me, airEnumStr(nrrdType, nrrdTypeBlock));
    biffAdd(NRRD, err); return 1;
  }
  if (2 != nmap->dim) {
    sprintf(err, "%s: map needs to have dimension 2, not %d", me, nmap->dim);
    biffAdd(NRRD, err); return 1;
  }
  entLen = nmap->axis[0].size;
  mapLen = nmap->axis[1].size;
  if (!( entLen >= 2 && mapLen >= 2 )) {
    sprintf(err, "%s: both map's axes sizes should be >= 2 (not %d,%d)",
	    me, entLen, mapLen);
    biffAdd(NRRD, err); return 1;
  }
  min[0] = 1; max[0] = nmap->axis[0].size-1;
  min[1] = 0; max[1] = nmap->axis[1].size-1;
  if (nrrdCrop(nrange=nrrdNew(), nmap, min, max)) {
    sprintf(err, "%s: couldn't crop to isolate range of map", me);
    biffAdd(NRRD, err); nrrdNuke(nrange); return 1;
  }
  if (nrrdHasNonExistSet(nrange)) {
    sprintf(err, "%s: map has non-existent values in its range", me);
    biffAdd(NRRD, err); nrrdNuke(nrange); return 1;
  }
  nrrdNuke(nrange);
  mapLup = nrrdDLookup[nmap->type];
  if (AIR_EXISTS(mapLup(nmap->data, 0))) {
    baseI = 0;
  } else {
    baseI = 3;
    if (!( mapLen >= 5 )) {
      sprintf(err, "%s: length of map w/ non-existant locations must "
	      "be >= 5 (not %d)", me, mapLen);
      biffAdd(NRRD, err); return 1;
    }
    if (!( airFP_NEG_INF == airFPClass_f(mapLup(nmap->data, 0*entLen)) &&
	   airFP_QNAN    == airFPClass_f(mapLup(nmap->data, 1*entLen)) &&
	   airFP_POS_INF == airFPClass_f(mapLup(nmap->data, 2*entLen)) )) {
      sprintf(err, "%s: 1st entry's position non-existant, but position "
	      "of 1st three entries not -inf, NaN, and +inf", me);
      biffAdd(NRRD, err); return 1;
    }
  }
  for (i=baseI; i<mapLen; i++) {
    if (!AIR_EXISTS(mapLup(nmap->data, i*entLen))) {
      sprintf(err, "%s: entry %d has non-existant position", me, i);
      biffAdd(NRRD, err); return 1;
    }
  }
  for (i=baseI; i<mapLen-1; i++) {
    if (!( mapLup(nmap->data, i*entLen) < mapLup(nmap->data, (i+1)*entLen) )) {
      sprintf(err, "%s: map entry %d pos (%g) not < entry %d pos (%g)",
	      me, i, mapLup(nmap->data, i*entLen),
	      i+1, mapLup(nmap->data, (i+1)*entLen));
      biffAdd(NRRD, err); return 1;
    }
  }
  return 0;
}

/*
** _nrrdApply1DSetUp()
**
** some error checking and initializing needed for 1D LUTS, regular,
** and irregular maps.  The intent is that if this succeeds, then
** there is no need for any further error checking.
**
** The only thing this function DOES is allocate the output nrrd,
** and set peripheral information.  The rest is just error checking.
*/
int
_nrrdApply1DSetUp(Nrrd *nout, Nrrd *nin, Nrrd *nmap,
		  int kind, int typeOut, int rescale) {
  char me[]="_nrrdApply1DSetUp", err[AIR_STRLEN_MED], *mapcnt;
  char nounStr[][AIR_STRLEN_SMALL]={"lut",
				    "regular map",
				    "irregular map"};
  char verbStr[][AIR_STRLEN_SMALL]={"lut",
				    "rmap",
				    "imap"};
  int mapAxis, size[NRRD_DIM_MAX], axisMap[NRRD_DIM_MAX], d, colLen;

  if (nout == nin) {
    sprintf(err, "%s: due to laziness, nout==nin always disallowed", me);
    biffAdd(NRRD, err); return 1;
  }
  if (airEnumValCheck(nrrdType, typeOut)) {
    sprintf(err, "%s: invalid requested output type %d", me, typeOut);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdTypeBlock == nin->type || nrrdTypeBlock == typeOut) {
    sprintf(err, "%s: input or requested output type is %s, need scalar",
	    me, airEnumStr(nrrdType, nrrdTypeBlock));
    biffAdd(NRRD, err); return 1;
  }
  if (rescale && !(AIR_EXISTS(nin->min) && AIR_EXISTS(nin->max))) {
    sprintf(err, "%s: want rescaling but not both nin->{min,max} exist", me);
    biffAdd(NRRD, err); return 1;
  }
  if (kindLut == kind || kindRmap == kind) {
    mapAxis = nmap->dim - 1;
    if (!(0 == mapAxis || 1 == mapAxis)) {
      sprintf(err, "%s: dimension of %s should be 1 or 2, not %d", 
	      me, nounStr[kind], nmap->dim);
      biffAdd(NRRD, err); return 1;
    }
    if (!(AIR_EXISTS(nmap->axis[mapAxis].min) &&
	  AIR_EXISTS(nmap->axis[mapAxis].max))) {
      sprintf(err, "%s: axis min and max not both set on axis %d of %s",
	      me, mapAxis, nounStr[kind]);
      biffAdd(NRRD, err); return 1;
    }
    if (!(nmap->axis[mapAxis].min < nmap->axis[mapAxis].max)) {
      sprintf(err, "%s: axis min (%g) not less than max (%g)", me,
	      nmap->axis[mapAxis].min, nmap->axis[mapAxis].max);
      biffAdd(NRRD, err); return 1;
    }
    if (nrrdHasNonExistSet(nmap)) {
      sprintf(err, "%s: %s nrrd has non-existent values",
	      me, nounStr[kind]);
      biffAdd(NRRD, err); return 1;
    }
    colLen = mapAxis ? nmap->axis[0].size : 1;
  } else {
    /* its an irregular map */
    if (nrrd1DIrregMapCheck(nmap)) {
      sprintf(err, "%s: problem with irregular map", me);
      biffAdd(NRRD, err); return 1;
    }
    /* mapAxis has no meaning for irregular maps, but we'll pretend ... */
    mapAxis = nmap->axis[0].size == 2 ? 0 : 1;
    colLen = nmap->axis[0].size-1;
  }
  if (mapAxis + nin->dim > NRRD_DIM_MAX) {
    sprintf(err, "%s: input nrrd dim %d through non-scalar %s exceeds "
	    "NRRD_DIM_MAX %d",
	    me, nin->dim, nounStr[kind], NRRD_DIM_MAX);
    biffAdd(NRRD, err); return 1;
  }
  nrrdAxesGet_nva(nin, nrrdAxesInfoSize, size+mapAxis);
  if (mapAxis) {
    size[0] = colLen;
    axisMap[0] = -1;
  }
  for (d=0; d<nin->dim; d++) {
    axisMap[d+mapAxis] = d;
  }
  /*
  fprintf(stderr, "##%s: pre maybe alloc: nout->data = %p\n", me, nout->data);
  for (d=0; d<mapAxis + nin->dim; d++) {
    fprintf(stderr, "    size[%d] = %d\n", d, size[d]);
  }
  fprintf(stderr, "   nout->dim = %d; nout->type = %d = %s; sizes = %d,%d\n", 
	  nout->dim, nout->type,
	  airEnumStr(nrrdType, nout->type));
  fprintf(stderr, "   typeOut = %d = %s\n", typeOut,
	  airEnumStr(nrrdType, typeOut));
  */
  if (nrrdMaybeAlloc_nva(nout, typeOut, mapAxis + nin->dim, size)) {
    sprintf(err, "%s: couldn't allocate output nrrd", me);
    biffAdd(NRRD, err); return 1;
  }
  /*
  fprintf(stderr, "   nout->dim = %d; nout->type = %d = %s\n",
	  nout->dim, nout->type,
	  airEnumStr(nrrdType, nout->type),
	  nout->axis[0].size, nout->axis[1].size);
  for (d=0; d<nout->dim; d++) {
    fprintf(stderr, "    size[%d] = %d\n", d, nout->axis[d].size);
  }
  fprintf(stderr, "##%s: post maybe alloc: nout->data = %p\n", me, nout->data);
  */
  if (nrrdAxesCopy(nout, nin, axisMap, NRRD_AXESINFO_NONE)) {
    sprintf(err, "%s: trouble copying axes", me);
    biffAdd(NRRD, err); return 1;
  }
  mapcnt = _nrrdContentGet(nmap);
  if (nrrdContentSet(nout, verbStr[kind], nin, "%s", mapcnt)) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); free(mapcnt); return 1;
  }
  free(mapcnt); 
  nrrdPeripheralInit(nout);
  nout->hasNonExist = nrrdNonExistFalse;   /* really! */
  return 0;
}

/*
** _nrrdApply1DLutOrRegMap()
**
** the guts of nrrdApply1DLut and nrrdApply1DRegMap
**
** yikes, does NOT use biff, since we're only supposed to be called
** after copious error checking.  
**
** FOR INSTANCE, this allows nout == nin, which could be a big
** problem if mapAxis == 1.
**
** we don't need a typeOut arg because nout has already been allocated
** as some specific type; we'll look at that.
*/
int
_nrrdApply1DLutOrRegMap(Nrrd *nout, Nrrd *nin, Nrrd *nmap,
			int interp, int rescale) {
  /* char me[]="nrrdApply1DLutOrRegMap" , err[AIR_STRLEN_MED] ;  */
  char *inData, *outData, *mapData, *entData0, *entData1;
  size_t N, I;
  double (*inLoad)(void *v), (*mapLup)(void *v, size_t I),
    (*outInsert)(void *v, size_t I, double d),
    val, mapIdxFrac, mapMin, mapMax;
  int i, mapAxis, mapLen, mapIdx, entSize, entLen, inSize, outSize;

  mapAxis = nmap->dim - 1;             /* axis of nmap containing entries */
  mapData = nmap->data;                /* map data, as char* */
  mapMin = nmap->axis[mapAxis].min;    /* low end of map domain */
  mapMax = nmap->axis[mapAxis].max;    /* high end of map domain */
  mapLen = nmap->axis[mapAxis].size;   /* number of entries in map */
  mapLup = nrrdDLookup[nmap->type];    /* how to get doubles out of map */
  inData = nin->data;                  /* input data, as char* */
  inLoad = nrrdDLoad[nin->type];       /* how to get doubles out of nin */
  inSize = nrrdElementSize(nin);       /* size of one input value */
  outData = nout->data;                /* output data, as char* */
  outInsert = nrrdDInsert[nout->type]; /* putting doubles into output */
  entLen = (mapAxis                    /* number of elements in one entry */
	    ? nmap->axis[0].size
	    : 1);
  outSize = entLen*nrrdElementSize(nout); /* size of entry in output */
  entSize = entLen*nrrdElementSize(nmap); /* size of entry in map */

  /* In below, we do quantize-then-clamp, because clamp-then-quantize
     won't correctly handle non-existant values */
  N = nrrdElementNumber(nin);       /* the number of values to be mapped */
  /* _VV = 0; */
  if (interp) {
    /* regular map */
    for (I=0; I<N; I++) {
      /* _VV = 0*(I > 73600); */
      /* if (_VV && !(I % 100)) fprintf(stderr, "I = %d\n", (int)I); */
      val = inLoad(inData);
      /* if (_VV) fprintf(stderr, "##%s: val = \na% 31.15f --> ", me, val); */
      if (rescale) {
	val = AIR_AFFINE(nin->min, val, nin->max, mapMin, mapMax);
	/* if (_VV) fprintf(stderr, "\nb% 31.15f --> ", val); */
      }
      /* if (_VV) fprintf(stderr, "\nc% 31.15f --> ", val); */
      val = AIR_CLAMP(mapMin, val, mapMax);
      mapIdxFrac = AIR_AFFINE(mapMin, val, mapMax, 0, mapLen-1);
      /* if (_VV) fprintf(stderr, "mapIdxFrac = \nd% 31.15f --> ",
	 mapIdxFrac); */
      mapIdx = mapIdxFrac;
      mapIdx -= mapIdx == mapLen-1;
      mapIdxFrac -= mapIdx;
      /* if (_VV) fprintf(stderr, "(%d,\ne% 31.15f) --> ", mapIdx,
	 mapIdxFrac); */
      entData0 = mapData + mapIdx*entSize;
      entData1 = mapData + (mapIdx+1)*entSize;
      for (i=0; i<entLen; i++) {
	val = ((1-mapIdxFrac)*mapLup(entData0, i) + 
	       mapIdxFrac*mapLup(entData1, i));
	outInsert(outData, i, val);
	/* if (_VV) fprintf(stderr, "\nf% 31.15f\n", val); */
      }
      inData += inSize;
      outData += outSize;
    }    
  } else {
    /* lookup table */
    for (I=0; I<N; I++) {
      val = inLoad(inData);
      if (rescale) {
	val = AIR_AFFINE(nin->min, val, nin->max, mapMin, mapMax);
      }
      AIR_INDEX(mapMin, val, mapMax, mapLen, mapIdx);
      mapIdx = AIR_CLAMP(0, mapIdx, mapLen-1);
      entData0 = mapData + mapIdx*entSize;
      for (i=0; i<entLen; i++) {
	outInsert(outData, i, mapLup(entData0, i));
      }
      inData += inSize;
      outData += outSize;
    }
  }

  return 0;
}

/*
******** nrrdApply1DLut
**
** A "lut" is a simple lookup table: the data points are evenly spaced,
** with cell-centering assumed, and there is no interpolation except
** nearest neighbor.  The axis min and max are used to determine the
** range of values that can be mapped with the lut.
**
** Of the three kinds of 1-D maps, only luts can have output type block.
**
** If the lut nrrd is 1-D, then the output is a scalar nrrd with the
** same dimension as the input.  If the lut nrrd is 2-D, then each
** value in the input is mapped to a vector of values from the lut,
** which is always a scanline along axis 0.  
**
** This allows lut length to be simply 1.
*/
int
nrrdApply1DLut(Nrrd *nout, Nrrd *nin, Nrrd *nlut,
	       int typeOut, int rescale) {
  char me[]="nrrdApply1DLut", err[AIR_STRLEN_MED];
  
  if (!(nout && nlut && nin)) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (_nrrdApply1DSetUp(nout, nin, nlut, kindLut, typeOut, rescale)) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); return 1;
  }
  if (_nrrdApply1DLutOrRegMap(nout, nin, nlut, AIR_FALSE, rescale)) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); return 1;
  }
  return 0;
}

/*
******** nrrdApply1DRegMap
**
** A regular map has data points evenly spaced, with node-centering and
** and linear interpolation assumed.  As with luts, the axis min and
** max determine the range of values that are mapped.  Previously, the
** intent was to allow interpolation with an arbitrary kernel, but
** that would be complicated to generalize to N-D.  Also, for the 
** foreseeable future, this function is to be used by nrrdHistoEq() and
** nothing else in teem, so the need for higher-order filtering is not
** accute.  Besides, it probably belongs in another library anyway.
**
** On second thought, eliminating arbitrary interpolation makes the
** existance of this function a little less motivated, given that it
** is a special case of nrrdApply1DIrregMap().
**
** If the lut nrrd is 1-D, then the output is a scalar nrrd with the
** same dimension as the input.  If the lut nrrd is 2-D, then each
** value in the input is mapped to a linear weighting of vectors 
** from the map; the vectors are the scanlines along axis 0.
**
** NB: this function makes NO provisions for non-existant input values.
** There won't be any memory errors, but the results are undefined.
**
** This allows map length to be simply 1.
*/
int
nrrdApply1DRegMap(Nrrd *nout, Nrrd *nin, Nrrd *nmap,
		  int typeOut, int rescale) {
  char me[]="nrrdApply1DRegMap", err[AIR_STRLEN_MED];

  if (!(nout && nmap && nin)) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (_nrrdApply1DSetUp(nout, nin, nmap, kindRmap, typeOut, rescale)) {
    sprintf(err, "%s: ", me);
    biffAdd(NRRD, err); return 1;
  }
  if (_nrrdApply1DLutOrRegMap(nout, nin, nmap, AIR_TRUE, rescale)) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); return 1;
  }
  return 0;
}

/*
******** nrrd1DIrregAclCheck()
**
** returns zero only on valid accelerators for 1D irregular mappings
*/
int
nrrd1DIrregAclCheck(Nrrd *nacl) {
  char me[]="nrrd1DIrregAclCheck", err[AIR_STRLEN_MED];

  if (!nacl) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdCheck(nacl)) {
    sprintf(err, "%s: ", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdTypeUShort != nacl->type) {
    sprintf(err, "%s: type should be %s, not %s", me,
	    airEnumStr(nrrdType, nrrdTypeUShort),
	    airEnumStr(nrrdType, nacl->type));
    biffAdd(NRRD, err); return 1;
  }
  if (2 != nacl->dim) {
    sprintf(err, "%s: dimension should be 2, not %d", me, nacl->dim);
    biffAdd(NRRD, err); return 1;
  }
  if (!( nacl->axis[0].size == 2 && nacl->axis[1].size >= 2 )) {
    sprintf(err, "%s: sizes (%d,%d) not (2,>=2)", me,
	    nacl->axis[0].size, nacl->axis[1].size);
    biffAdd(NRRD, err); return 1;
  }

  return 0;
}

/*
** _nrrd1DIrregMapDomain()
**
** ALLOCATES an array of doubles storing the existant control point
** locations, and sets its length in *poslenP.  If there are the three
** points with non-existant locations, these are ignored.
**
** Assumes that nrrd1DIrregMapCheck has been called on "nmap".
*/
double *
_nrrd1DIrregMapDomain(int *posLenP, int *baseIP, Nrrd *nmap) {
  char me[]="_nrrd1DIrregMapDomain", err[AIR_STRLEN_MED];
  int i, entLen, baseI, posLen;
  double *pos, (*mapLup)(void *v, size_t I);

  mapLup = nrrdDLookup[nmap->type];
  baseI = AIR_EXISTS(mapLup(nmap->data, 0)) ? 0 : 3;
  if (baseIP) {
    *baseIP = baseI;
  }
  entLen = nmap->axis[0].size;
  posLen = nmap->axis[1].size - baseI;
  if (posLenP) {
    *posLenP = posLen;
  }
  pos = (double*)malloc(posLen * sizeof(double));
  if (!pos) {
    sprintf(err, "%s: couldn't allocate %d doubles\n", me, posLen);
    biffAdd(NRRD, err); return NULL;
  }
  for (i=0; i<posLen; i++) {
    pos[i] = mapLup(nmap->data, (baseI+i)*entLen);
  }
  return pos;
}

/*
** _nrrd1DIrregFindInterval()
**
** The hard part of doing 1D irregular mapping: given an array of
** control point locations, and a value, find which interval the value
** lies in.  The lowest and highest possible indices are given in
** "loI" and "hiI".  Results are undefined if these do not in fact
** bound the location of correct interval, or if loI > hiI, or if the
** query positon "p" is not in the domain vector "pos".  Intervals are
** identified by the integral index of the LOWER of the two control
** points spanning the interval.
**
** This imposes the same structure of half-open intervals that
** is done by AIR_INDEX.  That is, a value p is in interval i
** if pos[i] <= p < pos[i+1] for all but the last interval, and
** pos[i] <= p <= pos[i+1] for the last interval (in which case
** i == hiI)
*/
int
_nrrd1DIrregFindInterval(double *pos, double p, int loI, int hiI) {
  int midI;

  /*
  fprintf(stderr, "##%s(%g): hi: %d/%g-%g | %d/%g-%g\n",
	  "_nrrd1DIrregFindInterval", p,
	  loI, pos[loI], pos[loI+1], 
	  hiI, pos[hiI], pos[hiI+1]);
  */
  while (loI < hiI) {
    midI = (loI + hiI)/2;
    if ( pos[midI] <= p && ((midI <  hiI && p <  pos[midI+1]) ||
			    (midI == hiI && p <= pos[midI+1])) ) {
      /* p is between (pos[midI],pos[midI+1]): we're done */
      loI = hiI = midI;
    } else if (pos[midI] > p) {
      /* p is below interval midI: midI-1 is valid upper bound */
      hiI = midI-1;
    } else {
      /* p is above interval midI: midI+1 is valid lower bound */
      loI = midI+1;
    }
    /*
    fprintf(stderr, "##%s(%g): %d/%g-%g | %d/%g-%g | %d/%g-%g\n",
	    "_nrrd1DIrregFindInterval", p,
	    loI, pos[loI], pos[loI+1], 
	    midI, pos[midI], pos[midI+1], 
	    hiI, pos[hiI], pos[hiI+1]);
    */
  }
  return loI;
}

/*
******** nrrd1DIrregAclGenerate()
**
** Generates the "acl" that is used to speed up the action of
** nrrdApply1DIrregMap().  Basically, the domain of the map
** is quantized into "acllen" bins, and for each bin, the
** lowest and highest possible map interval is stored. This
** either obviates or speeds up the task of finding which
** interval contains a given value.
**
** Assumes that nrrd1DIrregMapCheck has been called on "nmap".
*/
int
nrrd1DIrregAclGenerate(Nrrd *nacl, Nrrd *nmap, int aclLen) {
  char me[]="nrrd1DIrregAclGenerate", err[AIR_STRLEN_MED];
  int i, posLen;
  unsigned short *acl;
  double lo, hi, min, max, *pos;

  if (!(nacl && nmap)) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (!(aclLen >= 2)) {
    sprintf(err, "%s: given acl length (%d) is too small", me, aclLen);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdMaybeAlloc(nacl, nrrdTypeUShort, 2, 2, aclLen)) {
    sprintf(err, "%s: ", me);
    biffAdd(NRRD, err); return 1;
  }
  acl = nacl->data;
  pos = _nrrd1DIrregMapDomain(&posLen, NULL, nmap);
  if (!pos) {
    sprintf(err, "%s: couldn't determine domain", me); 
    biffAdd(NRRD, err); return 1;
  }
  nacl->axis[1].min = min = pos[0];
  nacl->axis[1].max = max = pos[posLen-1];
  for (i=0; i<=aclLen-1; i++) {
    lo = AIR_AFFINE(0, i, aclLen, min, max);
    hi = AIR_AFFINE(0, i+1, aclLen, min, max);
    acl[0 + 2*i] = _nrrd1DIrregFindInterval(pos, lo, 0, posLen-2);
    acl[1 + 2*i] = _nrrd1DIrregFindInterval(pos, hi, 0, posLen-2);
  }
  free(pos);

  return 0;
}

/*
******** nrrdApply1DIrregMap()
**
** Linear interpolation between irregularly spaced control points.
** Obviously, the location of the control point has to be given
** explicitly.  The map nrrd must have dimension 2, and each 
** control point is represented by a scanline along axis 0.  The
** first value is the position of the control point, and the remaining
** value(s) are linearly weighted according to the position of the
** input value among the control point locations.
**
** To allow "coloring" of non-existant values -inf, NaN, and +inf, if
** the very first value of the map (the location of the first control
** point) is non-existant, then the first three control point locations
** must be -inf, NaN, and +inf, in that order, and the information
** about these points will be used for corresponding input values.
** Doing this makes everything slower, however, because airFPClass_f()
** is called on every single value.
**
** This assumes that nrrd1DIrregMapCheck has been called on "nmap",
** and that nrrd1DIrregAclCheck has been called on "nacl" (if it is
** non-NULL).
*/
int
nrrdApply1DIrregMap(Nrrd *nout, Nrrd *nin, Nrrd *nmap, Nrrd *nacl,
		    int typeOut, int rescale) {
  char me[]="nrrdApply1DIrregMap", err[AIR_STRLEN_MED];
  size_t N, I;
  int i, *acl, entLen, posLen, aclLen, mapIdx, aclIdx,
    entSize, colSize, inSize, lo, hi, baseI;
  double val, *pos, mapMin, mapMax, mapIdxFrac,
    (*mapLup)(void *v, size_t I),
    (*inLoad)(void *v), (*outInsert)(void *v, size_t I, double d);
  char *inData, *outData, *entData0, *entData1;

  if (!(nout && nmap && nin)) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (_nrrdApply1DSetUp(nout, nin, nmap, kindImap, typeOut, rescale)) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nacl && nrrd1DIrregAclCheck(nacl)) {
    sprintf(err, "%s: given acl isn't valid", me);
    biffAdd(NRRD, err); return 1;
  }
  
  if (nacl) {
    acl = nacl->data;
    aclLen = nacl->axis[1].size;
  } else {
    acl = NULL;
    aclLen = 0;
  }
  pos = _nrrd1DIrregMapDomain(&posLen, &baseI, nmap);
  if (!pos) {
    sprintf(err, "%s: couldn't determine domain", me); 
    biffAdd(NRRD, err); return 1;
  }
  mapLup = nrrdDLookup[nmap->type];
  
  inData = nin->data;
  inLoad = nrrdDLoad[nin->type];
  inSize = nrrdElementSize(nin);
  mapLup = nrrdDLookup[nmap->type];
  entLen = nmap->axis[0].size;    /* entLen is really 1 + entry length */
  entSize = entLen*nrrdElementSize(nmap);
  colSize = (entLen-1)*nrrdTypeSize[typeOut];
  outData = nout->data;
  outInsert = nrrdDInsert[nout->type];
  mapMin = pos[0];
  mapMax = pos[posLen-1];
  
  N = nrrdElementNumber(nin);
  for (I=0;
       I<N;
       I++, inData += inSize, outData += colSize) {
    val = inLoad(inData);
    /* _VV = ( (AIR_EXISTS(val) && (21 == (int)(-val))) 
       || 22400 < I ); */
    /* if (_VV)
       fprintf(stderr, "##%s: (%d) val = % 31.15f\n", me, (int)I, val); */
    if (!AIR_EXISTS(val)) {
      /* got a non-existant value */
      if (baseI) {
	/* and we know how to deal with them */
	switch (airFPClass_f(val)) {
	case airFP_NEG_INF:
	  mapIdx = 0;
	  break;
	case airFP_SNAN:
	case airFP_QNAN:
	  mapIdx = 1;
	  break;
	case airFP_POS_INF:
	  mapIdx = 2;
	  break;
	default:
	  mapIdx = 0;
	  fprintf(stderr, "%s: PANIC: non-existant value/class %g/%d "
		  "not handled\n",
		  me, val, airFPClass_f(val));
	  exit(1);
	}
	entData0 = (char*)(nmap->data) + mapIdx*entSize;
	for (i=1; i<entLen; i++) {
	  outInsert(outData, i-1, mapLup(entData0, i));
	}
	continue;  /* we're done! (with this value) */
      } else {
	/* we don't know how to properly deal with this non-existant value:
	   we use the first entry, and then fall through to code below */
	mapIdx = 0;
	mapIdxFrac = 0.0;
      }
    }
    else {
      /* we have an existant value */
      if (rescale) {
	val = AIR_AFFINE(nin->min, val, nin->max, mapMin, mapMax);
	/* if (_VV) fprintf(stderr, "   rescaled --> % 31.15f\n", val); */
      }
      val = AIR_CLAMP(mapMin, val, mapMax);
      if (acl) {
	AIR_INDEX(mapMin, val, mapMax, aclLen, aclIdx);
	lo = acl[0 + 2*aclIdx];
	hi = acl[1 + 2*aclIdx];
      } else {
	lo = 0;
	hi = posLen-2;
      }
      if (lo < hi) {
	mapIdx = _nrrd1DIrregFindInterval(pos, val, lo, hi);
      } else {
	/* acl did its job ==> lo == hi */
	mapIdx = lo;
      }
    }
    mapIdxFrac = AIR_AFFINE(pos[mapIdx], val, pos[mapIdx+1], 0.0, 1.0);
    /* if (_VV) fprintf(stderr, "##%s: val=\n% 31.15f --> "
		     "mapIdx,frac = %d,\n% 31.15f\n",
		     me, val, mapIdx, mapIdxFrac); */
    entData0 = (char*)(nmap->data) + (baseI+mapIdx)*entSize;
    entData1 = (char*)(nmap->data) + (baseI+mapIdx+1)*entSize;
    /* if (_VV) fprintf(stderr, "##%s: 2; %d/\n% 31.15f --> entLen=%d "
		     "baseI=%d -->\n",
		     me, mapIdx, mapIdxFrac, entLen, baseI); */
    for (i=1; i<entLen; i++) {
      val = ((1-mapIdxFrac)*mapLup(entData0, i) +
	     mapIdxFrac*mapLup(entData1, i));
      /* if (_VV) fprintf(stderr, "% 31.15f\n", val); */
      outInsert(outData, i-1, val);
    }
    /* if (_VV) fprintf(stderr, "##%s: 3\n", me); */
  }
  free(pos);
  return 0;
}
