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
#include "../ell.h"

char *me;

void
usage() {
  /*                      0   1   2   3   (4) */
  fprintf(stderr, "usage: %s <A> <B> <C>\n", me);
  fprintf(stderr, "for cubic x^3 + Ax^2 + Bx + C == 0\n");
  exit(1);
}

int
main(int argc, char **argv) {
  char buf[512];
  double ans0, ans1, ans2, A, B, C;
  int ret;
  double r[3];

  me = argv[0];
  if (argc != 4) {
    usage();
  }

  sprintf(buf, "%s %s %s", argv[1], argv[2], argv[3]);
  if (3 != sscanf(buf, "%lf %lf %lf", &A, &B, &C)) {
    fprintf(stderr, "%s: couldn't parse 3 floats from command line\n", me);
    exit(1);
  }

  ellDebug = AIR_TRUE;
  ret = ellCubic(r, A, B, C, AIR_TRUE);
  ans0 = C + r[0]*(B + r[0]*(A + r[0]));
  switch(ret) {
  case ellCubicRootSingle:
    printf("1 single root: %g -> %f\n", r[0], ans0);
    break;
  case ellCubicRootTriple:
    printf("1 triple root: %g -> %f\n", r[0], ans0);
    break;
  case ellCubicRootSingleDouble:
    ans1 = C + r[1]*(B + r[1]*(A + r[1]));
    printf("1 single root %g -> %f, 1 double root %g -> %f\n", 
	   r[0], ans0, r[1], ans1);
    break;
  case ellCubicRootThree:
    ans1 = C + r[1]*(B + r[1]*(A + r[1]));
    ans2 = C + r[2]*(B + r[2]*(A + r[2]));
    printf("3 distinct roots:\n %g -> %f\n %g -> %f\n %g -> %f\n",
	   r[0], ans0, r[1], ans1, r[2], ans2);
    break;
  default:
    printf("%s: something fatally wacky happened\n", me);
    exit(1);
  }
  exit(0);
}
