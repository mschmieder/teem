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

#include "private.h"

char *histaxName = "histax";
#define INFO "Replace each scanline along an axis with its histogram"
char *histaxInfo = INFO;
char *histaxInfoL = (INFO
		   ". ");

int
histaxMain(int argc, char **argv, char *me) {
  hestOpt *opt = NULL;
  char *out, *err;
  Nrrd *nin, *nout;
  int axis, type, bins;
  airArray *mop;

  OPT_ADD_NIN(nin, "input nrrd");
  OPT_ADD_AXIS(axis, "axis to histogram along");
  hestOptAdd(&opt, "b", "bins", airTypeInt, 1, 1, &bins, NULL,
	     "# of bins in histogram");
  OPT_ADD_TYPE(type, "output type");
  OPT_ADD_NOUT(out, "output nrrd");

  mop = airMopInit();
  airMopAdd(mop, opt, (airMopper)hestOptFree, airMopAlways);

  USAGE(histaxInfoL);
  PARSE();
  airMopAdd(mop, opt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);

  if (nrrdHistoAxis(nout, nin, axis, bins, type)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error doing axis histogramming:\n%s", me, err);
    airMopError(mop);
    return 1;
  }
  

  SAVE(nout, NULL);

  airMopOkay(mop);
  return 0;
}
