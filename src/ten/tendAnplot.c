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

#include "ten.h"
#include "tenPrivate.h"

#define INFO "Graph anisotropy metric in barycentric coords"
char *_tend_anplotInfoL =
  (INFO
   ".  The metrics all vary from 0.0 to 1.0, and will be sampled "
   "in the lower right half of the image.");

int
tend_anplotMain(int argc, char **argv, char *me, hestParm *hparm) {
  int pret;
  hestOpt *hopt = NULL;
  char *perr, *err;
  airArray *mop;

  int res, aniso;
  Nrrd *nout;
  char *outS;

  hestOptAdd(&hopt, "r", "res", airTypeInt, 1, 1, &res, "256",
	     "resolution of anisotropy plot");
  hestOptAdd(&hopt, "a", "aniso", airTypeEnum, 1, 1, &aniso, NULL,
	     "Which anisotropy metric to plot.  " TEN_ANISO_DESC,
	     NULL, tenAniso);
  hestOptAdd(&hopt, "o", "nout", airTypeString, 1, 1, &outS, NULL,
	     "output image (floating point)");

  mop = airMopInit();
  airMopAdd(mop, hopt, (airMopper)hestOptFree, airMopAlways);
  USAGE(_tend_anplotInfoL);
  PARSE();
  airMopAdd(mop, hopt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);
  if (tenAnisoPlot(nout, aniso, res)) {
    airMopAdd(mop, err=biffGetDone(TEN), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble making plot:\n%s\n", me, err);
    airMopError(mop); exit(1);
  }
  
  if (nrrdSave(outS, nout, NULL)) {
    airMopAdd(mop, err=biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: trouble writing:\n%s\n", me, err);
    airMopError(mop); exit(1);
  }

  airMopOkay(mop);
  return 0;
}
/* TEND_CMD(anplot, INFO); */
unrrduCmd tend_anplotCmd = { "anplot", INFO, tend_anplotMain };

