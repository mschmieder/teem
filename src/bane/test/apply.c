/*
  The contents of this file are subject to the University of Utah Public
  License (the "License"); you may not use this file except in
  compliance with the License.
  
  Software distributed under the License is distributed on an "AS IS"
  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied.  See
  the License for the specific language governing rights and limitations
  under the License.

  The Original Source Code is "teem", released March 23, 2001.
  
  The Original Source Code was developed by the University of Utah.
  Portions created by UNIVERSITY are Copyright (C) 2001, 1998 University
  of Utah. All Rights Reserved.
*/


#include "../bane.h"

char *me;

void
usage() {
  /*                      0     1     2       3   (4) */
  fprintf(stderr, "usage: %s <nin> <measr> <nout>\n", me);
  exit(1);
}

int
main(int argc, char *argv[]) {
  int measr;
  char *iStr, *mStr, *oStr;
  Nrrd *nin, *nout;
  
  me = argv[0];
  if (4 != argc)
    usage();

  iStr = argv[1];
  mStr = argv[2];
  oStr = argv[3];
  
  if (!(nin = nrrdNewLoad(iStr))) {
    fprintf(stderr, "%s: trouble reading input nrrd:\n%s\n", me, 
	    biffGet(NRRD));
    usage();
  }

  if (1 != sscanf(mStr, "%d", &measr)) {
    fprintf(stderr, "%s: couldn't parse %s as int\n", me, mStr);
    usage();
  }
  
  if (baneApplyMeasr(nout = nrrdNew(), nin, measr)) {
    fprintf(stderr, "%s: trouble:\n%s\n", me, biffGet(BANE));
    exit(1);
  }

  if (nrrdSave(oStr, nout)) {
    fprintf(stderr, "%s: trouble writing output nrrd:\n%s\n", me,
	    biffGet(NRRD));
    usage();
  }
  exit(0);
}
