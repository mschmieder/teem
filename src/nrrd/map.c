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
#include "private.h"

/*
******** nrrdMinMaxSet()
**
** Sets nrrd->min and nrrd->max to the extremal (existant) values in
** the given nrrd, by calling the appropriate member of nrrdFindMinMax[]
**
** calling this function will result in nrrd->hasNonExist being set
** (because of the nrrdFindMinMax[] functions)
**
** decided NOT to use biff, so that this is a more distinct alternative to 
** nrrdMinMaxCleverSet().
*/
void
nrrdMinMaxSet(Nrrd *nrrd) {
  NRRD_TYPE_BIGGEST _min, _max;

  if (nrrd) {
    if (airEnumValidVal(nrrdType, nrrd->type)
	&& nrrdTypeBlock != nrrd->type) {
      nrrdFindMinMax[nrrd->type](&_min, &_max, nrrd);
      nrrd->min = nrrdDLoad[nrrd->type](&_min);
      nrrd->max = nrrdDLoad[nrrd->type](&_max);
    } else {
      nrrd->min = nrrd->max = AIR_NAN;
    }
  }
  return;
}

/*
** nrrdMinMaxCleverSet()
**
** basically a wrapper around nrrdMinMaxSet(), with bells + whistles:
** 1) will call nrrdMinMaxSet only when one of nrrd->min and nrrd->max
**    are non-existent, with the end result that only the non-existent
**    values are over-written
** 2) obeys the nrrdStateClever8BitMinMax global state to short-cut
**    finding min and max for 8-bit data.  Values for nrrd->min or 
**    nrrd->max which were existant to start with are untouched.
** 3) reports error if there are no existent values in nrrd (AIR_EXISTS()
**    fails on every value)
**
** Like nrrdMinMaxSet(), this will always set nrrd->hasNonExist.
**
** Uses biff.
*/
int
nrrdMinMaxCleverSet(Nrrd *nrrd) {
  char me[]="nrrdMinMaxCleverSet", err[AIR_STRLEN_MED];
  double min, max;

  if (!nrrd) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (!airEnumValidVal(nrrdType, nrrd->type)) {
    sprintf(err, "%s: input nrrd has invalid type (%d)", me, nrrd->type);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdTypeBlock == nrrd->type) {
    sprintf(err, "%s: can't find min/max of type %s", me,
	    airEnumStr(nrrdType, nrrdTypeBlock));
  }
  if (AIR_EXISTS(nrrd->min) && AIR_EXISTS(nrrd->max)) {
    /* both of min and max already set, so we won't look for those, but
       we have to comply with stated behavior of always setting hasNonExist */
    nrrdHasNonExistSet(nrrd);
    return 0;
  }
  if (nrrdStateClever8BitMinMax
      && (nrrdTypeChar == nrrd->type || nrrdTypeUChar == nrrd->type)) {
    if (nrrdTypeChar == nrrd->type) {
      if (!AIR_EXISTS(nrrd->min))
	nrrd->min = SCHAR_MIN;
      if (!AIR_EXISTS(nrrd->max))
	nrrd->max = SCHAR_MAX;
    } else {
      if (!AIR_EXISTS(nrrd->min))
	nrrd->min = 0;
      if (!AIR_EXISTS(nrrd->max))
	nrrd->max = UCHAR_MAX;
    }
    nrrdHasNonExistSet(nrrd);
    return 0;
  }

  /* at this point we need to find either min and/or max (at least
     one of them was non-existent on the way in) */

  /* save incoming values in case they exist */
  min = nrrd->min;
  max = nrrd->max;
  /* this will set nrrd->min, nrrd->max, and hasNonExist */
  nrrdMinMaxSet(nrrd);
  if (!( AIR_EXISTS(nrrd->min) && AIR_EXISTS(nrrd->max) )) {
    sprintf(err, "%s: no existent values!", me);
    biffAdd(NRRD, err); return 1;
  }
  /* re-enstate the existent incoming min and/or max values */
  if (AIR_EXISTS(min))
    nrrd->min = min;
  if (AIR_EXISTS(max))
    nrrd->max = max;

  return 0;
}

/*
******** nrrdConvert()
**
** copies values from one type of nrrd to another, without any
** transformation, except what you get with a cast.  The point is to
** make available on Nrrds the exact same behavior as you have in C
** with casts and assignments.
*/
int
nrrdConvert(Nrrd *nout, Nrrd *nin, int type) {
  char me[] = "nrrdConvert", typeS[AIR_STRLEN_SMALL], err[AIR_STRLEN_MED];
  int size[NRRD_DIM_MAX];
  size_t num;

  if (!( nin && nout 
	 && airEnumValidVal(nrrdType, nin->type)
	 && airEnumValidVal(nrrdType, type) )) {
    sprintf(err, "%s: invalid args", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nin->type == nrrdTypeBlock || type == nrrdTypeBlock) {
    sprintf(err, "%s: can't convert to or from nrrd type %s", me,
	    airEnumStr(nrrdType, nrrdTypeBlock));
    biffAdd(NRRD, err); return 1;
  }
  if (nout == nin && nrrdTypeSize[type] != nrrdTypeSize[nin->type]) {
    sprintf(err, "%s: nout==nin but input,output type sizes unequal", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdStateDisallowFixedPointNonExist 
      && !nrrdTypeFixed[nin->type] && nrrdTypeFixed[type]) {
    /* there's a risk of non-existant values getting converted to
       non-sensical fixed point values */
    if (nrrdHasNonExistSet(nin)) {
      sprintf(err, "%s: can't convert to fixed point with "
	      "non-existant values in input nrrd", me);
      biffAdd(NRRD, err); return 1;
    }
  }

  /* if we're actually converting to the same type, just do a copy */
  if (type == nin->type) {
    if (nout == nin) {
      /* nout == nin is allowed if the input and output type are
	 of the same size, which will certainly be the case if the
	 input and output types are identical, so there's actually
	 no work to do */
    } else {
      if (nrrdCopy(nout, nin)) {
	sprintf(err, "%s: couldn't copy input to output", me);
	biffAdd(NRRD, err); return 1;
      }
    }
    return 0;
  }

  /* allocate space if necessary */
  nrrdAxesGet_nva(nin, nrrdAxesInfoSize, size);
  /* MUST be nrrdMaybeAlloc_nva (not nrrdAlloc_nva) because we allow
     nout==nin if type sizes match */
  if (nrrdMaybeAlloc_nva(nout, type, nin->dim, size)) {
    sprintf(err, "%s: failed to allocate output", me);
    biffAdd(NRRD, err); return 1;
  }

  /* call the appropriate converter */
  num = nrrdElementNumber(nin);
  _nrrdConv[nout->type][nin->type](nout->data, nin->data, num);

  /* copy peripheral information */
  nrrdAxesCopy(nout, nin, NULL, NRRD_AXESINFO_NONE);
  sprintf(typeS, "(%s)", airEnumStr(nrrdType, nout->type));
  if (nrrdContentSet(nout, typeS, nin, "")) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); return 1;
  }
  /* the min and max have probably changed if there was a conversion
     to fixed point, or to a lower precision representation */
  nrrdPeripheralInit(nout);
  nout->blockSize = 0;
  nout->hasNonExist = (nrrdTypeFixed[nout->type]
		       ? nrrdNonExistFalse
		       : nin->hasNonExist);

  /* bye */
  return 0;
}

/*
******** nrrdQuantize()
**
** convert values to 8, 16, or 32 bit unsigned quantities
** by mapping the value range delimited by the nrrd's min
** and max to the representable range 
**
** NOTE: for the time being, this uses a "double" as the intermediate
** value holder, which may mean needless loss of precision
*/
int
nrrdQuantize(Nrrd *nout, Nrrd *nin, int bits) {
  char me[] = "nrrdQuantize", func[]="quantize", err[AIR_STRLEN_MED];
  double valIn, min, max, eps;
  int valOut, type=nrrdTypeUnknown, size[NRRD_DIM_MAX];
  unsigned long long int valOutll;
  size_t I, num;
  unsigned char *outUC;
  unsigned short *outUS;
  unsigned int *outUI;

  if (!(nin && nout)) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (!(8 == bits || 16 == bits || 32 == bits)) {
    sprintf(err, "%s: bits has to be 8, 16, or 32 (not %d)", me, bits);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdTypeBlock == nin->type) {
    sprintf(err, "%s: can't quantize type %s", me,
	    airEnumStr(nrrdType, nrrdTypeBlock));
  }

  /* determine nrrd type from number of bits */
  switch (bits) {
  case 8:  type = nrrdTypeUChar;  break;
  case 16: type = nrrdTypeUShort; break;
  case 32: type = nrrdTypeUInt;   break;
  }
  if (nout == nin && nrrdTypeSize[type] != nrrdTypeSize[nin->type]) {
    sprintf(err, "%s: nout==nin but input,output type sizes unequal", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdMinMaxCleverSet(nin)) {
    sprintf(err, "%s: trouble setting min, max", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdStateDisallowFixedPointNonExist && nin->hasNonExist) {
    sprintf(err, "%s: can't quantize non-existent values (NaN, +/-inf)", me);
    biffAdd(NRRD, err); return 1;
  }
  
  /* allocate space if necessary */
  nrrdAxesGet_nva(nin, nrrdAxesInfoSize, size);
  /* MUST be nrrdMaybeAlloc_nva (not nrrdAlloc_nva) because we allow
     nout==nin if type sizes match */
  if (nrrdMaybeAlloc_nva(nout, type, nin->dim, size)) {
    sprintf(err, "%s: failed to create output", me);
    biffAdd(NRRD, err); return 1;
  }

  /* the skinny */
  num = nrrdElementNumber(nin);
  min = nin->min; 
  max = nin->max;
  eps = (min == max ? 1.0 : 0.0);
  outUC = (unsigned char*)nout->data;
  outUS = (unsigned short*)nout->data;
  outUI = (unsigned int*)nout->data;
  switch(bits) {
  case 8:
    for (I=0; I<num; I++) {
      valIn = nrrdDLookup[nin->type](nin->data, I);
      valIn = AIR_CLAMP(min, valIn, max);
      AIR_INDEX(min, valIn, max+eps, 1 << 8, valOut);
      outUC[I] = valOut;
    }
    break;
  case 16:
    for (I=0; I<num; I++) {
      valIn = nrrdDLookup[nin->type](nin->data, I);
      valIn = AIR_CLAMP(min, valIn, max);
      AIR_INDEX(min, valIn, max+eps, 1 << 16, valOut);
      outUS[I] = valOut;
    }
    break;
  case 32:
    for (I=0; I<num; I++) {
      valIn = nrrdDLookup[nin->type](nin->data, I);
      valIn = AIR_CLAMP(min, valIn, max);
      AIR_INDEX(min, valIn, max+eps, 1LLU << 32, valOutll);
      outUI[I] = valOutll;
    }
    break;
  }

  /* set information in new volume */
  /* nrrdPeripheralInit(nout); nout==nin qualms */
  if (nout != nin) {
    nrrdAxesCopy(nout, nin, NULL, NRRD_AXESINFO_NONE);
  }
  nout->oldMin = min;
  nout->oldMax = max;
  if (nrrdContentSet(nout, func, nin, "%d", bits)) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); return 1;
  }
  nout->min = nout->max = AIR_NAN;
  nout->blockSize = 0;
  nout->hasNonExist = nrrdNonExistFalse;

  return 0;
}

/*
** _nrrdHistoEqCompare()
**
** used by nrrdHistoEq in smart mode to sort the "steady" array
** in _descending_ order
*/
int 
_nrrdHistoEqCompare(const void *a, const void *b) {

  return *((unsigned int*)b) - *((unsigned int*)a);
}

/*
******** nrrdHistoEq()
**
** performs histogram equalization on given nrrd, treating it as a
** big one-dimensional array.  The procedure is as follows: 
** - create a histogram of nrrd (using "bins" bins)
** - integrate the histogram, and normalize and shift this so it is 
**   a monotonically increasing function from min to max, where
**   (min,max) is the range of values in the nrrd
** - map the values in the nrrd through the adjusted histogram integral
** 
** If the histogram of the given nrrd is already as flat as can be,
** the histogram integral will increase linearly, and the adjusted
** histogram integral should be close to the identity function, so
** the values shouldn't change much.
**
** If the nhistP arg is non-NULL, then it is set to point to
** the histogram that was used for calculation. Otherwise this
** histogram is deleted on return.
**
** This is all that is done normally, when "smart" is <= 0.  In
** "smart" mode (activated by setting "smart" to something greater
** than 0), the histogram is analyzed during its creation to detect if
** there are a few bins which keep getting hit with the same value
** over and over.  It may be desirable to ignore these bins in the
** histogram integral because they may not contain any useful
** information, and so they should not effect how values are
** re-mapped.  The value of "smart" is the number of bins that will be
** ignored.  For instance, use the value 1 if the problem with naive
** histogram equalization is a large amount of background (which is
** exactly one fixed value).  
*/
int
nrrdHistoEq(Nrrd *nout, Nrrd *nin, Nrrd **nmapP, int bins, int smart) {
  char me[]="nrrdHistoEq", func[]="heq", err[AIR_STRLEN_MED];
  Nrrd *nhist, *nmap;
  float *ycoord = NULL;
  double val, min, max, *last = NULL;
  int i, idx, *respect = NULL, lort, hirt;
  unsigned int *hist, *steady = NULL;
  size_t I, num;
  airArray *mop;

  if (!(nout && nin)) {
    sprintf(err, "%s: got NULL pointer", me);
    biffAdd(NRRD, err); return 1;
  }
  if (nrrdTypeBlock == nin->type) {
    sprintf(err, "%s: can't histogram equalize type %s", me,
	    airEnumStr(nrrdType, nrrdTypeBlock));
    biffAdd(NRRD, err); return 1;
  }
  if (!(bins > 4)) {
    sprintf(err, "%s: need # bins > 4 (not %d)", me, bins);
    biffAdd(NRRD, err); return 1;
  }
  /* we start by simply copying, because the only thing we're 
     changing is the values themselves, and all peripheral
     information is unchanged by this value remapping, even
     min and max. */
  if (nout != nin) {
    if (nrrdCopy(nout, nin)) {
      sprintf(err, "%s:", me); 
      biffAdd(NRRD, err); return 1;
    }
  }
  mop = airMopInit();
  if (nmapP) {
    airMopAdd(mop, *nmapP, (airMopper)airSetNull, airMopOnError);
  }
  num = nrrdElementNumber(nin);
  if (smart <= 0) {
    nhist = nrrdNew();
    if (nrrdHisto(nhist, nin, bins, nrrdTypeUInt)) {
      sprintf(err, "%s: failed to create histogram", me);
      biffAdd(NRRD, err); airMopError(mop); return 1;
    }
    airMopAdd(mop, nhist, (airMopper)nrrdNuke, airMopAlways);
    hist = nhist->data;
    min = nhist->axis[0].min;
    max = nhist->axis[0].max;
  } else {
    /* for "smart" mode, we have to some extra work while creating the
       histogram to look for bins incessantly hit with the exact same
       value */
    if (nrrdAlloc(nhist=nrrdNew(), nrrdTypeUInt, 1, bins)) {
      sprintf(err, "%s: failed to allocate histogram", me);
      biffAdd(NRRD, err); airMopError(mop); return 1;
    }
    airMopAdd(mop, nhist, (airMopper)nrrdNuke, airMopAlways);
    hist = nhist->data;
    nhist->axis[0].size = bins;
    /* allocate the respect, steady, and last arrays */
    respect = calloc(bins, sizeof(int));
    steady = calloc(2*bins, sizeof(unsigned int));
    last = calloc(bins, sizeof(double));
    airMopMem(mop, &respect, airMopAlways);
    airMopMem(mop, &steady, airMopAlways);
    airMopMem(mop, &last, airMopAlways);
    if (!(respect && steady && last)) {
      sprintf(err, "%s: couldn't allocate smart arrays", me);
      biffAdd(NRRD, err); airMopError(mop); return 1;
    }
    /* steady[0 + 2*i] == how many times has bin i seen the same value
       steady[1 + 2*i] == i (steady will be rearranged by qsort()) */
    for (i=0; i<bins; i++) {
      last[i] = AIR_NAN;
      respect[i] = 1;
      steady[1 + 2*i] = i;
    }
    /* now create the histogram */
    nin->min = nin->max = AIR_NAN;
    if (nrrdMinMaxCleverSet(nin)) {
      sprintf(err, "%s: couldn't find value range in nrrd", me);
      biffAdd(NRRD, err); airMopError(mop); return 1;
    }
    min = nin->min;
    max = nin->max;
    for (I=0; I<num; I++) {
      val = nrrdDLookup[nin->type](nin->data, I);
      if (AIR_EXISTS(val)) {
	AIR_INDEX(min, val, max, bins, idx);
	++hist[idx];
	if (AIR_EXISTS(last[idx])) {
	  steady[0 + 2*idx] = (last[idx] == val
			       ? 1 + steady[0 + 2*idx]
			       : 0);
	}
	last[idx] = val;
      }
    }
    /*
    for (i=0; i<bins; i++) {
      printf("steady(%d) = %d\n", i, steady[0 + 2*i]);
    }
    */
    /* now sort the steady array */
    qsort(steady, bins, 2*sizeof(unsigned int), _nrrdHistoEqCompare);
    /*
    for (i=0; i<=20; i++) {
      printf("sorted steady(%d/%d) = %d\n", i, steady[1+2*i], steady[0+2*i]);
    }
    */
    /* we ignore some of the bins according to "smart" arg */
    for (i=0; i<smart; i++) {
      respect[steady[1+2*i]] = 0;
      /* printf("%s: disrespecting bin %d\n", me, steady[1+2*i]); */
    }
  }
  if (nrrdAlloc(nmap=nrrdNew(), nrrdTypeFloat, 1, bins+1)) {
    sprintf(err, "%s: failed to create map nrrd", me);
    biffAdd(NRRD, err); airMopError(mop); return 1;
  }
  airMopAdd(mop, nmap, (airMopper)nrrdNuke,
	    nmapP ? airMopOnError : airMopAlways);
  nmap->axis[0].min = min;
  nmap->axis[0].max = max;
  ycoord = nmap->data;

  /* integrate the histogram then normalize it */
  for (i=0; i<=bins; i++) {
    if (i == 0) {
      ycoord[i] = 0;
    } else {
      ycoord[i] = ycoord[i-1] + hist[i-1]*(smart 
					   ? respect[i-1] 
					   : 1);
    }
  }
  /* if we've done smart, the integral will have little flat spots
     where we ignored hits in the histogram.  That means the mapping
     will actually be lower at that value than it should be.  In
     truth, we should be using an irregular mapping for this, and the
     control points at the ignored bins should just be missing.  So we
     have to do this silliness to raise those control points in the
     regular map. */
  if (smart) {
    /* there are bins+1 control points, with indices 0 to bins.
       We'll fix control points 1 to bins-1.  ycoord[i] is too low
       if hist[i-1] was not respected (!respect[i-1]) */
    for (i=1; i<=bins-1; i++) {
      if (!respect[i-1]) {
	/* lort and hirt will bracket the index of the bad control point
	   with points corresponding either to respected bins or the
	   endpoints of the histogram */
	for (lort=i; lort>=1 && !respect[lort-1]; lort--);
	for (hirt=i; hirt<=bins-1 && !respect[hirt-1]; hirt++);
	ycoord[i] = AIR_AFFINE(lort, i, hirt, ycoord[lort], ycoord[hirt]);
      }
    }
    /* the very last control point has to be handled differently */
    if (!respect[bins-1]) {
      ycoord[bins] += ycoord[bins-1] - ycoord[bins-2];
    }
  }
  for (i=0; i<=bins; i++) {
    ycoord[i] = AIR_AFFINE(0, ycoord[i], ycoord[bins], min, max);
  }

  /* map the nrrd values through the normalized histogram integral */
  if (nrrdApply1DRegMap(nout, nin, nmap, nin->type, AIR_FALSE)) {
    sprintf(err, "%s: problem remapping", me);
    biffAdd(NRRD, err); airMopError(mop); return 1;
  }
  /*
  for (I=0; I<num; I++) {
    val = nrrdDLookup[nin->type](nin->data, I);
    if (AIR_EXISTS(val)) {
      AIR_INDEX(min, val, max, bins, idx);
      val = AIR_AFFINE(xcoord[idx], val, xcoord[idx+1], 
		       ycoord[idx], ycoord[idx+1]);
    }
    nrrdDInsert[nout->type](nout->data, I, val);
  }
  */

  /* if user is interested, set pointer to map nrrd,
     otherwise it will be nixed by airMop */
  if (nmapP) {
    *nmapP = nmap;
  }

  /* fiddling with content is the only thing we'll do */
  if (nrrdContentSet(nout, func, nin, "%d,%d", bins, smart)) {
    sprintf(err, "%s:", me);
    biffAdd(NRRD, err); airMopError(mop); return 1;
  }
  nout->min = nin->min;
  nout->max = nin->max;

  airMopOkay(mop);
  return 0;
}
