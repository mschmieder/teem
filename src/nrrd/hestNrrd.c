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

#include "nrrd.h"
#include "privateNrrd.h"

/* ---------------------------- Nrrd ----------------------------- */

int
_nrrdHestNrrdParse(void *ptr, char *str, char err[AIR_STRLEN_HUGE]) {
  char me[] = "_nrrdHestNrrdParse", *nerr;
  Nrrd **nrrdP;
  airArray *mop;
  
  if (!(ptr && str)) {
    sprintf(err, "%s: got NULL pointer", me);
    return 1;
  }
  nrrdP = (Nrrd **)ptr;
  mop = airMopInit();
  *nrrdP = nrrdNew();
  airMopAdd(mop, *nrrdP, (airMopper)nrrdNuke, airMopOnError);
  if (nrrdLoad(*nrrdP, str)) {
    airMopAdd(mop, nerr = biffGetDone(NRRD), airFree, airMopOnError);
    strncpy(err, nerr, AIR_STRLEN_HUGE-1);
    airMopError(mop);
    return (strstr(err, _NRRD_IMM_EOF) ? 2 : 1);
  }
  airMopOkay(mop);
  return 0;
}

hestCB
_nrrdHestNrrd = {
  sizeof(Nrrd *),
  "nrrd",
  _nrrdHestNrrdParse,
  (airMopper)nrrdNuke
}; 

hestCB *
nrrdHestNrrd = &_nrrdHestNrrd;

/* ------------------------ NrrdKernelSpec -------------------------- */

int
_nrrdHestKernelSpecParse(void *ptr, char *str, char err[AIR_STRLEN_HUGE]) {
  NrrdKernelSpec **ksP;
  char me[]="_nrrdHestKernelSpecParse", *nerr;

  if (!(ptr && str)) {
    sprintf(err, "%s: got NULL pointer", me);
    return 1;
  }
  ksP = (NrrdKernelSpec **)ptr;
  *ksP = nrrdKernelSpecNew();
  if (nrrdKernelParse(&((*ksP)->kernel), (*ksP)->parm, str)) {
    nerr = biffGetDone(NRRD);
    strncpy(err, nerr, AIR_STRLEN_HUGE-1);
    free(nerr);
    return 1;
  }
  return 0;
}

hestCB
_nrrdHestKernelSpec = {
  sizeof(NrrdKernelSpec*),
  "kernel specification",
  _nrrdHestKernelSpecParse,
  (airMopper)nrrdKernelSpecNix
};

hestCB *
nrrdHestKernelSpec = &_nrrdHestKernelSpec;

/* --------------------------- NrrdIter ----------------------------- */

int
_nrrdLooksLikeANumber(char *str) {
  char c;
  /* 0: -+                (no restriction, but that's a little daft)
     1: 0123456789        n > 0
     2: .                 0 <= n <= 1
     3: eE                0 <= n <= 1
     4: everything else   0 == n
  */
  int count[5];
  
  count[0] = count[1] = count[2] = count[3] = count[4] = 0;
  while (*str) {
    c = tolower(*str);
    switch (c) {
    case '-': case '+':
      count[0]++;
      break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      count[1]++;
      break;
    case '.':
      count[2]++;
      break;
    case 'e':
      count[3]++;
      break;
    default:
      count[4]++;
      break;
    }
    str++;
  }
  if (count[1] > 0 &&
      AIR_INSIDE(0, count[2], 1) &&
      AIR_INSIDE(0, count[3], 1) &&
      count[4] == 0) {
    return AIR_TRUE;
  } else {
    return AIR_FALSE;
  }
}

int
_nrrdHestIterParse(void *ptr, char *str, char err[AIR_STRLEN_HUGE]) {
  char me[]="_nrrdHestIterParse", *nerr;
  Nrrd *nrrd;
  NrrdIter **iterP;
  airArray *mop;
  double val;
  int ret;
  
  if (!(ptr && str)) {
    sprintf(err, "%s: got NULL pointer", me);
    return 1;
  }
  iterP = (NrrdIter **)ptr;
  mop = airMopInit();
  *iterP = nrrdIterNew();
  airMopAdd(mop, *iterP, (airMopper)nrrdIterNuke, airMopOnError);

  /* the challenge here is determining if a given string represents a
     filename or a number.  Obviously there are cases where it could
     be both, so we'll assume its a filename first.  Because: there
     are different ways of writing the same number, such as "3" -->
     "+3", "3.1" --> "3.10", so someone accidently using the file when
     they mean to use the number has easy ways of changing the number
     representation.  On the other hand, to change the file name
     representation, they could prefix it with "./".  Another problem
     is that one really wants a general robust test to see if a given
     string is a valid number representation AND NOTHING BUT THAT, and
     sscanf() is not that test.  In any case, if there are to be
     improved smarts about this matter, they need to be implemented
     below and nowhere else. */
  
  nrrd = nrrdNew();
  ret = nrrdLoad(nrrd, str);
  if (!ret) {
    /* first attempt at nrrdLoad() was SUCCESSFUL */
    nrrdIterSetNrrd(*iterP, nrrd);
  } else {
    /* so it didn't load as a nrrd- if its because fopen() failed,
       then we'll try it as a number.  If its for another reason,
       then we complain */
    nrrdNuke(nrrd);
    if (2 != ret) {
      /* it failed because of something besides the fopen(), so complain */
      nerr = biffGetDone(NRRD);
      strncpy(err, nerr, AIR_STRLEN_HUGE-1);
      return 1;
    } else {
      /* fopen() failed, so it probably wasn't meant to be a filename */
      free(biffGetDone(NRRD));
      if (_nrrdLooksLikeANumber(str)) {
	/* printf("|%s| looks like a number\n", str); */
	if (1 == airSingleSscanf(str, "%lf", &val)) {
	  nrrdIterSetValue(*iterP, val);
	} else {
	  /* oh, this is bad. */
	  fprintf(stderr, "%s: PANIC, is it a number or not?", me);
	  exit(1);
	}
      } else {
	/* it doesn't look like a number, but the fopen failed, so
	   we'll let it fail again and pass back the error messages */
	if (nrrdLoad(nrrd = nrrdNew(), str)) {
	  nerr = biffGetDone(NRRD);
	  strncpy(err, nerr, AIR_STRLEN_HUGE-1);
	  return 1;
	} else {
	  /* what the hell? */
	  fprintf(stderr, "%s: PANIC, is it a nrrd or not?", me);
	  exit(1);
	}
      }
    }
  }
  airMopAdd(mop, iterP, (airMopper)airSetNull, airMopOnError);
  airMopOkay(mop);
  return 0;
}

hestCB
_nrrdHestIter = {
  sizeof(NrrdIter *),
  "nrrd/value",
  _nrrdHestIterParse,
  (airMopper)nrrdIterNuke
}; 

hestCB *
nrrdHestIter = &_nrrdHestIter;
