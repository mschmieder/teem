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

/* NB: not the same as char sliceName[] = "slice"; Read your C FAQs */
char *sliceName = "slice";
#define INFO "Slice at a position along an axis"
char *sliceInfo = INFO;
char *sliceInfoL = (INFO
		    ". Output nrrd dimension is one less than input nrrd "
		    "dimension.  Per-axis information is preserved.");

int
sliceMain(int argc, char **argv, char *me) {
  hestOpt *opt = NULL;
  char *out, *err;
  Nrrd *nin, *nout;
  int axis, _pos[2], pos;
  airArray *mop;

  OPT_ADD_NIN(nin, "input nrrd");
  OPT_ADD_AXIS(axis, "axis to slice along");
  hestOptAdd(&opt, "p", "pos", airTypeOther, 1, 1, _pos, NULL,
	     "position to slice at:\n "
	     "\b\bo <int> gives 0-based index\n "
	     "\b\bo M-<int> give index relative "
	     "to the last sample on the axis (M == #samples-1).",
	     NULL, NULL, &unuPosHestCB);
  OPT_ADD_NOUT(out, "output nrrd");

  mop = airMopInit();
  airMopAdd(mop, opt, (airMopper)hestOptFree, airMopAlways);

  USAGE(sliceInfoL);
  PARSE();
  airMopAdd(mop, opt, (airMopper)hestParseFree, airMopAlways);
  if (!( AIR_INSIDE(0, axis, nin->dim-1) )) {
    fprintf(stderr, "%s: axis %d not in range [0,%d]\n", me, axis, nin->dim-1);
    return 1;
  }
  pos = _pos[0]*(nin->axis[axis].size-1) + _pos[1];

  nout = nrrdNew();
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);

  if (nrrdSlice(nout, nin, axis, pos)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error slicing nrrd:\n%s", me, err);
    airMopError(mop);
    return 1;
  }

  SAVE(out, nout, NULL);

  airMopOkay(mop);
  return 0;
}
