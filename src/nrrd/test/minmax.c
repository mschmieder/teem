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


#include "../nrrd.h"

void
usage(char *me) { 
  /*                      0    1  (2) */
  fprintf(stderr, "usage: %s <nin>\n", me);
  exit(1);
}

int
main(int argc, char **argv) {
  char *me, *err;
  Nrrd *nrrd;

  me = argv[0];
  if (2 != argc)
    usage(me);

  nrrdStateVerboseIO = 10;
  
  if (nrrdLoad(nrrd=nrrdNew(), argv[1], NULL)) {
    fprintf(stderr, "%s: trouble loading \"%s\":\n%s", 
	    me, argv[1], err = biffGet(NRRD));
    free(err);
    exit(1);
  }

  if (nrrdMinMaxCleverSet(nrrd)) {
    fprintf(stderr, "%s: trouble finding min/max \"%s\":\n%s", 
	    me, argv[1], err = biffGet(NRRD));
    free(err);
    exit(1);
  }

  printf("%s: min = %g; max = %g, nrrd->hasNonExist = %d\n", 
	 me, nrrd->min, nrrd->max, nrrd->hasNonExist);

  nrrdNuke(nrrd);

  exit(0);
}
