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


#include <stdio.h>
#include <biff.h>
#include <hest.h>
#include <nrrd.h>
#include <gage.h>
#include <ten.h>

#define SPACING(spc) (AIR_EXISTS(spc) ? spc: nrrdDefSpacing)

/* copied this from ten.h; I don't want gage to depend on ten */
#define PROBE_MAT2LIST(l, m) ( \
   (l)[1] = (m)[0],          \
   (l)[2] = (m)[3],          \
   (l)[3] = (m)[6],          \
   (l)[4] = (m)[4],          \
   (l)[5] = (m)[7],          \
   (l)[6] = (m)[8] )

int
probeParseKind(void *ptr, char *str, char err[AIR_STRLEN_HUGE]) {
  char me[] = "probeParseKind";
  gageKind **kindP;
  
  if (!(ptr && str)) {
    sprintf(err, "%s: got NULL pointer", me);
    return 1;
  }
  kindP = ptr;
  airToLower(str);
  if (!strcmp("scalar", str)) {
    *kindP = gageKindScl;
  } else if (!strcmp("vector", str)) {
    *kindP = gageKindVec;
  } else if (!strcmp("tensor", str)) {
    *kindP = tenGageKind;
  } else {
    sprintf(err, "%s: not \"scalar\", \"vector\", or \"tensor\"", me);
    return 1;
  }

  return 0;
}

hestCB probeKindHestCB = {
  sizeof(gageKind *),
  "kind",
  probeParseKind,
  NULL
}; 

char *probeInfo = ("Shows off the functionality of the gage library. "
		   "Uses gageProbe() to query scalar or vector volumes "
		   "to learn various measured or derived quantities. ");

int
main(int argc, char *argv[]) {
  gageKind *kind;
  char *me, *outS, *whatS, *err;
  hestParm *hparm;
  hestOpt *hopt = NULL;
  NrrdKernelSpec *k00, *k11, *k22;
  float x, y, z, scale[3];
  int what, a, idx, ansLen, E=0, xi, yi, zi, otype,
    six, siy, siz, sox, soy, soz, iBaseDim, oBaseDim, renorm;
  gage_t *answer;
  Nrrd *nin, *nout;
  gageContext *ctx;
  gagePerVolume *pvl;
  double t0, t1, gmc;
  airArray *mop;

  mop = airMopNew();
  me = argv[0];
  hparm = hestParmNew();
  airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
  hparm->elideSingleOtherType = AIR_TRUE;
  hestOptAdd(&hopt, "i", "nin", airTypeOther, 1, 1, &nin, NULL,
	     "input volume", NULL, NULL, nrrdHestNrrd);
  hestOptAdd(&hopt, "k", "kind", airTypeOther, 1, 1, &kind, NULL,
	     "\"kind\" of volume (\"scalar\", \"vector\", or \"tensor\")",
	     NULL, NULL, &probeKindHestCB);
  hestOptAdd(&hopt, "q", "query", airTypeString, 1, 1, &whatS, NULL,
	     "the quantity (scalar, vector, or matrix) to learn by probing");
  hestOptAdd(&hopt, "s", "sclX sclY sxlZ", airTypeFloat, 3, 3, scale,
	     "1.0 1.0 1.0",
	     "scaling factor for resampling on each axis "
	     "(>1.0 : supersampling)");
  hestOptAdd(&hopt, "k00", "kern00", airTypeOther, 1, 1, &k00,
	     "tent", "kernel for gageKernel00",
	     NULL, NULL, nrrdHestKernelSpec);
  hestOptAdd(&hopt, "k11", "kern11", airTypeOther, 1, 1, &k11,
	     "cubicd:1,0", "kernel for gageKernel11",
	     NULL, NULL, nrrdHestKernelSpec);
  hestOptAdd(&hopt, "k22", "kern22", airTypeOther, 1, 1, &k22,
	     "cubicdd:1,0", "kernel for gageKernel22",
	     NULL, NULL, nrrdHestKernelSpec);
  hestOptAdd(&hopt, "rn", NULL, airTypeInt, 0, 0, &renorm, NULL,
	     "renormalize kernel weights at each new sample location. "
	     "\"Accurate\" kernels don't need this; doing it always "
	     "makes things go slower");
  hestOptAdd(&hopt, "gmc", "min gradmag", airTypeDouble, 1, 1, &gmc,
	     "0.0", "For curvature-based queries, use zero when gradient "
	     "magnitude is below this");
  hestOptAdd(&hopt, "t", "type", airTypeEnum, 1, 1, &otype, "float",
	     "type of output volume", NULL, nrrdType);
  hestOptAdd(&hopt, "o", "nout", airTypeString, 1, 1, &outS, NULL,
	     "output volume");
  hestParseOrDie(hopt, argc-1, argv+1, hparm,
		 me, probeInfo, AIR_TRUE, AIR_TRUE, AIR_TRUE);
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  what = airEnumVal(kind->enm, whatS);
  if (-1 == what) {
    /* -1 indeed always means "unknown" for any gageKind */
    fprintf(stderr, "%s: couldn't parse \"%s\" as measure of \"%s\" volume\n",
	    me, whatS, kind->name);
    hestUsage(stderr, hopt, me, hparm);
    hestGlossary(stderr, hopt, hparm);
    airMopError(mop);
    return 1;
  }

  ansLen = kind->ansLength[what];
  iBaseDim = gageKindScl == kind ? 0 : 1;
  oBaseDim = 1 == ansLen ? 0 : 1;
  six = nin->axis[0+iBaseDim].size;
  siy = nin->axis[1+iBaseDim].size;
  siz = nin->axis[2+iBaseDim].size;
  sox = scale[0]*six;
  soy = scale[1]*siy;
  soz = scale[2]*siz;
  nin->axis[0+iBaseDim].spacing = SPACING(nin->axis[0+iBaseDim].spacing);
  nin->axis[1+iBaseDim].spacing = SPACING(nin->axis[1+iBaseDim].spacing);
  nin->axis[2+iBaseDim].spacing = SPACING(nin->axis[2+iBaseDim].spacing);

  /***
  **** Except for the gageProbe() call in the inner loop below,
  **** and the gageContextNix() call at the very end, all the gage
  **** calls which set up the context and state are here.
  ***/
  ctx = gageContextNew();
  airMopAdd(mop, ctx, (airMopper)gageContextNix, airMopAlways);
  gageSet(ctx, gageParmGradMagCurvMin, gmc);
  gageSet(ctx, gageParmVerbose, 1);
  gageSet(ctx, gageParmRenormalize, renorm ? AIR_TRUE : AIR_FALSE);
  gageSet(ctx, gageParmCheckIntegrals, AIR_TRUE);
  E = 0;
  if (!E) E |= !(pvl = gagePerVolumeNew(ctx, nin, kind));
  if (!E) E |= gagePerVolumeAttach(ctx, pvl);
  if (!E) E |= gageKernelSet(ctx, gageKernel00, k00->kernel, k00->parm);
  if (!E) E |= gageKernelSet(ctx, gageKernel11, k11->kernel, k11->parm); 
  if (!E) E |= gageKernelSet(ctx, gageKernel22, k22->kernel, k22->parm);
  if (!E) E |= gageQuerySet(ctx, pvl, 1 << what);
  if (!E) E |= gageUpdate(ctx);
  if (E) {
    airMopAdd(mop, err = biffGetDone(GAGE), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble:\n%s\n", me, err);
    airMopError(mop);
    return 1;
  }
  answer = gageAnswerPointer(ctx, pvl, what);
  gageSet(ctx, gageParmVerbose, 0);
  /***
  **** end gage setup.
  ***/

  fprintf(stderr, "%s: kernel support = %d^3 samples\n", me,
	  2*(ctx->havePad + 1));
  fprintf(stderr, "%s: effective scaling is %g %g %g\n", me,
	  (float)sox/six, (float)soy/siy, (float)soz/siz);
  if (ansLen > 1) {
    fprintf(stderr, "%s: creating %d x %d x %d x %d output\n", 
	   me, ansLen, sox, soy, soz);
    if (!E) E |= nrrdMaybeAlloc(nout=nrrdNew(), otype, 4,
				ansLen, sox, soy, soz);
  } else {
    fprintf(stderr, "%s: creating %d x %d x %d output\n", me, sox, soy, soz);
    if (!E) E |= nrrdMaybeAlloc(nout=nrrdNew(), otype, 3, sox, soy, soz);
  }
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (E) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble:\n%s\n", me, err);
    airMopError(mop);
    return 1;
  }
  t0 = airTime();
  for (zi=0; zi<=soz-1; zi++) {
    fprintf(stderr, "%d/%d ", zi, soz-1); fflush(stderr);
    z = AIR_AFFINE(0, zi, soz-1, 0, siz-1);
    for (yi=0; yi<=soy-1; yi++) {
      y = AIR_AFFINE(0, yi, soy-1, 0, siy-1);
      for (xi=0; xi<=sox-1; xi++) {
	x = AIR_AFFINE(0, xi, sox-1, 0, six-1);
	idx = xi + sox*(yi + soy*zi);
	ctx->verbose = 0*( (!xi && !yi && !zi) ||
			   /* ((100 == xi) && (8 == yi) && (8 == zi)) */
			   ((61 == xi) && (51 == yi) && (46 == zi))
			   /* ((40==xi) && (30==yi) && (62==zi)) || */
			   /* ((40==xi) && (30==yi) && (63==zi)) */ ); 
	if (gageProbe(ctx, x, y, z)) {
	  fprintf(stderr, 
		  "%s: trouble at i=(%d,%d,%d) -> f=(%g,%g,%g):\n%s\n(%d)\n",
		  me, xi, yi, zi, x, y, z, gageErrStr, gageErrNum);
	  airMopError(mop);
	  return 1;
	}
	if (1 == ansLen) {
	  nrrdFInsert[nout->type](nout->data, idx,
				  nrrdFClamp[nout->type](*answer));
	} else {
	  for (a=0; a<=ansLen-1; a++) {
	    nrrdFInsert[nout->type](nout->data, a + ansLen*idx, 
				    nrrdFClamp[nout->type](answer[a]));
	  }
	}
      }
    }
  }

  /* HEY: this isn't actually correct in general, but is true
     for gageKindScl and gageKindVec */
  nrrdContentSet(nout, "probe", nin, "%s", airEnumStr(kind->enm, what));
  nout->axis[0+oBaseDim].spacing = 
    ((double)six/sox)*SPACING(nin->axis[0+iBaseDim].spacing);
  nout->axis[0+oBaseDim].label = airStrdup(nin->axis[0+iBaseDim].label);
  nout->axis[1+oBaseDim].spacing = 
    ((double)six/sox)*SPACING(nin->axis[1+iBaseDim].spacing);
  nout->axis[1+oBaseDim].label = airStrdup(nin->axis[1+iBaseDim].label);
  nout->axis[2+oBaseDim].spacing = 
    ((double)six/sox)*SPACING(nin->axis[2+iBaseDim].spacing);
  nout->axis[2+oBaseDim].label = airStrdup(nin->axis[2+iBaseDim].label);

  fprintf(stderr, "\n");
  t1 = airTime();
  fprintf(stderr, "probe rate = %g KHz\n", sox*soy*soz/(1000.0*(t1-t0)));
  if (nrrdSave(outS, nout, NULL)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble saving output:\n%s\n", me, err);
    airMopError(mop);
    return 1;
  }

  airMopOkay(mop);
  return 0;
}
