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

#include "unrrdu.h"
#include "privateUnrrdu.h"

/* bad bad bad Gordon */
extern int _nrrdOneLine(int *lenP, NrrdIO *io, FILE *file);

#define INFO "Print header of one or more nrrd files"
char *_unrrdu_headInfoL = 
(INFO  ".  The value of this is simply to print the contents of nrrd "
 "headers.  This avoids the use of \"head -N\", where N has to be "
 "determined manually, which always risks printing raw binary data "
 "(following the header) to screen, which tends to clobber terminal "
 "settings, as well as be annoying.");

int
unrrdu_headDoit(char *me, NrrdIO *io, char *inS, FILE *fout) {
  char err[AIR_STRLEN_MED];
  airArray *mop;
  int len, magic;
  FILE *fin;

  mop = airMopNew();
  if (!( fin = airFopen(inS, stdin, "rb") )) {
    sprintf(err, "%s: couldn't fopen(\"%s\",\"rb\"): %s\n", 
	    me, inS, strerror(errno));
    biffAdd(me, err); airMopError(mop); return 1;
  }
  airMopAdd(mop, fin, (airMopper)airFclose, airMopAlways);

  if (_nrrdOneLine(&len, io, fin)) {
    sprintf(err, "%s: error getting first line of file \"%s\"", me, inS);
    biffAdd(me, err); airMopError(mop); return 1;
  }
  if (!len) {
    sprintf(err, "%s: immediately hit EOF\n", me);
    biffAdd(me, err); airMopError(mop); return 1;
  }
  magic = airEnumVal(nrrdMagic, io->line);
  if (!( nrrdMagicOldNRRD == magic || nrrdMagicNRRD0001 == magic )) {
    sprintf(err, "%s: first line (\"%s\") isn't a nrrd magic\n", 
	    me, io->line);
    biffAdd(me, err); airMopError(mop); return 1;
  }
  while (len > 1) {
    fprintf(fout, "%s\n", io->line);
    _nrrdOneLine(&len, io, fin);
  };
  
#ifdef _WIN32
  /* seems that only on windows does the writing process's fwrite() to
     stdout fail if we exit without consuming everything from stdin */
  c = fgetc(fin);
  while (EOF != c) {
    c = fgetc(fin);
  }
#endif

  airMopOkay(mop);
  return 0;
}

int
unrrdu_headMain(int argc, char **argv, char *me, hestParm *hparm) {
  hestOpt *opt = NULL;
  char *err, **inS, *outS="-";
  NrrdIO *io;
  airArray *mop;
  int pret, ni, ninLen;
  FILE *fout;
#ifdef _WIN32
  int c;
#endif

  mop = airMopNew();
  hestOptAdd(&opt, NULL, "nin1", airTypeString, 1, -1, &inS, NULL,
	     "input nrrd(s)", &ninLen);
  airMopAdd(mop, opt, (airMopper)hestOptFree, airMopAlways);

  USAGE(_unrrdu_headInfoL);
  PARSE();
  airMopAdd(mop, opt, (airMopper)hestParseFree, airMopAlways);

  io = nrrdIONew();
  airMopAdd(mop, io, (airMopper)nrrdIONix, airMopAlways);

  if (!( fout = airFopen(outS, stdout, "wb") )) {
    sprintf(err, "%s: couldn't fopen(\"%s\",\"wb\"): %s\n", 
	    me, outS, strerror(errno));
    biffAdd(me, err); airMopError(mop); return 1;
  }
  airMopAdd(mop, fout, (airMopper)airFclose, airMopAlways);

  for (ni=0; ni<ninLen; ni++) {
    if (ninLen > 1) {
      fprintf(fout, "==> %s <==\n", inS[ni]);
    }
    /* HEY: would be better if we continued reading from other
       files after one is bad, but oh well */
    if (unrrdu_headDoit(me, io, inS[ni], fout)) {
      airMopAdd(mop, err = biffGetDone(me), airFree, airMopAlways);
      fprintf(stderr, "%s: trouble reading from \"%s\":\n%s",
	      me, inS[ni], err);
      airMopError(mop); return 1;
    }
    if (ninLen > 1 && ni < ninLen-1) {
      fprintf(fout, "\n");
    }
  }

  airMopOkay(mop);
  return 0;
}

UNRRDU_CMD(head, INFO);
