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


#include "../nrrd.h"

int
main(int argc, char *argv[]) {
  int tt;
  /*
  char str1[AIR_STRLEN_MED], str2[AIR_STRLEN_MED], 
    str3[AIR_STRLEN_MED];
  */

  for (tt=nrrdTypeChar; tt<=nrrdTypeDouble; tt++) {
    printf(" ----- %s -----\n", airEnumStr(nrrdType, tt));
    printf("nrrdTypeSize: %d\n", nrrdTypeSize[tt]);
    printf("nrrdTypeUnsigned: %d\n", nrrdTypeUnsigned[tt]);
    printf("nrrdTypeFixed: %d\n", nrrdTypeFixed[tt]);
    printf("nrrdTypeMin: % 31.15f\n", nrrdTypeMin[tt]);
    printf("nrrdTypeMax: % 31.15f\n", nrrdTypeMax[tt]);
    printf("nrrdTypeNumberValues: % 31.15f\n", nrrdTypeNumberValues[tt]);
  }
  printf("sizeof(size_t) = %d\n", (int)sizeof(size_t));

  /*
  c = -10;
  uc = 10;
  s = -10;
  us = 10;
  i = -10;
  ui = 10;
  lli = -10;
  ulli = 10;
  f = 3.14159324234098320948172304987123;
  d = 3.14159324234098320948172304987123;
  printf("c: %d\n", c);
  printf("uc: %u\n", uc);
  printf("s: %hd\n", s);
  printf("us: %hu\n", us);
  printf("i: %d\n", i);
  printf("ui: %u\n", ui);
  printf("lli: %lld\n", lli);
  printf("ulli: %llu\n", ulli);
  printf("f: %f\n", f);
  printf("d: %lf\n", d);

  sprintf(str1, "-10");
  sprintf(str2, "10");
  sprintf(str3, "3.14159324234098320948172304987123");

  sscanf(str1, "%d", &c);
  sscanf(str2, "%u", &uc);
  sscanf(str1, "%hd", &s);
  sscanf(str2, "%hu", &us);
  sscanf(str1, "%d", &i);
  sscanf(str2, "%u", &ui);
  sscanf(str1, "%lld", &lli);
  sscanf(str2, "%llu", &ulli);
  sscanf(str3, "%f", &f);
  sscanf(str3, "%lf", &d);

  printf("\n");
  printf("c: %d\n", c);
  printf("uc: %u\n", uc);
  printf("s: %hd\n", s);
  printf("us: %hu\n", us);
  printf("i: %d\n", i);
  printf("ui: %u\n", ui);
  printf("lli: %lld\n", lli);
  printf("ulli: %llu\n", ulli);
  printf("f: %f\n", f);
  printf("d: %lf\n", d);
  */

  exit(0);
}
