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


#include <stdio.h>
#include <air.h>
#include "../hest.h"

int
main(int argc, char **argv) {
  int res[2], v, numIn;
  char **in, *out, *blah[3], *option = NULL;
  int n, *ints, numN;
  hestOpt *opt = NULL;
  hestParm *parm;
  char *err = NULL, info[] = 
    "This program does nothing in particular, though it does attempt "
    "to pose as some sort of command-line image processing program. "
    "As usual, any implied functionality is purely coincidental, "
    "especially since this is the output of a unicyclist.";

  parm = hestParmNew();
  parm->respFileEnable = AIR_TRUE;
  parm->verbosity = 3;

  opt = NULL;
  hestOptAdd(&opt, "v,verbose",     "level", airTypeInt,    0,  1,  &v,   "0",
	     "verbosity level");
  hestOptAdd(&opt, "out",   "file",  airTypeString, 1,  1,  &out, "output.ppm",
	     "PPM image output");
  hestOptAdd(&opt, "blah",  "input", airTypeString, 3,  3,  blah,  "a b c",
	     "input image file(s)");
  hestOptAdd(&opt, "option","opt", airTypeString, 0, 1, &option, "default",
	     "this is just a test");
  /*
  hestOptAdd(&opt, NULL,    "input", airTypeString, 1, -1,  &in,  NULL,
	     "input image file(s)", &numIn);
  */
  hestOptAdd(&opt, "ints",  "N",     airTypeInt,    1,  -1, &ints,  "10 20 30",
	     "a list of integers", &numN);
  hestOptAdd(&opt, "res",   "sx sy", airTypeInt,    2,  2,  res,  NULL,
	     "image resolution");
  
  printf("what 0\n");
  if (1 == argc) {
    /* didn't get anything at all on command line */
    /* print program information ... */
    hestInfo(stderr, argv[0], info, parm);
    /* ... and usage information ... */
    hestUsage(stderr, opt, argv[0], parm);
    hestGlossary(stderr, opt, parm);
    /* ... and avoid memory leaks */
    opt = hestOptFree(opt);
    parm = hestParmFree(parm);
    exit(1);
  }

  printf("what 1\n");

  /* else we got something, see if we can parse it */
  if (hestParse(opt, argc-1, argv+1, &err, parm)) {
    fprintf(stderr, "ERROR: %s\n", err); free(err);
    /* print usage information ... */
    hestUsage(stderr, opt, argv[0], parm);
    hestGlossary(stderr, opt, parm);
    /* ... and then avoid memory leaks */
    opt = hestOptFree(opt);
    parm = hestParmFree(parm);
    printf(" ---- in = %lx\n", (unsigned long)in);
    printf(" ---- blah[0] = %lx\n", (unsigned long)(blah[0]));
    printf(" ---- option = %lx\n", (unsigned long)option);
    exit(1);
  }
  printf("what 2\n");

  printf("(err = %s)\n", err ? err : "(null)");
  printf("  v = %d\n", v);
  printf("out = \"%s\"\n", out ? out : "(null)");
  printf("blah = \"%s\" \"%s\" \"%s\"\n", blah[0], blah[1], blah[2]);
  printf("option = \"%s\"\n", option ? option : "(null)");
  printf("res = %d %d\n", res[0], res[1]);
  /*
  printf(" ---- in = %lx\n", (unsigned long)in);
  printf(" in = %d files:", numIn);
  for (n=0; n<=numIn-1; n++) {
    printf(" \"%s\"", in[n] ? in[n] : "(null)");
  }
  printf("\n");
  */
  printf(" ints = %d ints:", numN);
  for (n=0; n<=numN-1; n++) {
    printf(" %d", ints[n]);
  }
  printf("\n");
  printf("what 3\n");

  /* free the memory allocated by parsing ... */
  hestParseFree(opt);
  /* ... and the other stuff */
  opt = hestOptFree(opt);
  parm = hestParmFree(parm);
  exit(0);
}
