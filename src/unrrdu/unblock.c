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

#include "privateUnrrdu.h"

char *unblockName = "unblock";
#define INFO "Expand \"blocks\" into scanlines on axis 0"
char *unblockInfo = INFO;
char *unblockInfoL = (INFO
		      ". Based on the requested output type, the number of "
		      "samples along axis 0 will be determined automatically. "
		      "Axis N information will be bumped up to axis N+1. "
		      "Underlying data is unchanged.");

int
unblockMain(int argc, char **argv, char *me) {
  hestOpt *opt = NULL;
  char *out, *err;
  Nrrd *nin, *nout;
  int type, blockSize;
  airArray *mop;

  OPT_ADD_NIN(nin, "input nrrd");
  OPT_ADD_TYPE(type, "type to unblock to", NULL);
  hestOptAdd(&opt, "bs", "blocksize", airTypeInt, 1, 1, &blockSize, "0",
	     "Useful only if *output* type is also block: the size of "
	     "blocks in output nrrd");
  OPT_ADD_NOUT(out, "output nrrd");

  mop = airMopInit();
  airMopAdd(mop, opt, (airMopper)hestOptFree, airMopAlways);

  USAGE(unblockInfo);
  PARSE();
  airMopAdd(mop, opt, (airMopper)hestParseFree, airMopAlways);

  nout = nrrdNew();
  nout->blockSize = blockSize;
  airMopAdd(mop, nout, (airMopper)nrrdNuke, airMopAlways);

  if (nrrdUnblock(nout, nin, type)) {
    airMopAdd(mop, err = biffGetDone(NRRD), airFree, airMopAlways);
    fprintf(stderr, "%s: error unblocking nrrd:\n%s", me, err);
    airMopError(mop);
    return 1;
  }

  SAVE(out, nout, NULL);

  airMopOkay(mop);
  return 0;
}
