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


#include "biff.h"

/*
** _biffEntry struct
**
** hold information and messages associated with one key
*/
typedef struct {
  char key[BIFF_MAXKEYLEN+1]; /* the key */
  char **err;                 /* array of error strings; the err array itself
				 is NOT null-terminated */
  int num;                    /* length of "err" == # strings stored */
  airArray *AA;               /* air array for err and num */
} _biffEntry;

_biffEntry **_biffErr ;       /* master array of _biffEntry pointers */
int _biffNum,                 /* length of _biffErr == # keys maintained */
  _biffIdx;                   /* hack: index of latest key found */
airArray *_biffAA;            /* air array of _biffErr and _biffNum */

#define _BIFF_INCR 2

/*
** _biffInit()
**
** allocates data structers needed by biff.  Panics and exit(1)s if 
** anything goes wrong.  Can be harmlessly called multiple times.
*/
void
_biffInit() {
  char me[]="_biffInit";

  if (!_biffAA) {
    _biffAA = airArrayNew((void**)&_biffErr, &_biffNum, 
			  sizeof(_biffEntry*), _BIFF_INCR);
    if (!_biffAA) {
      fprintf(stderr, "%s: PANIC: couldn't allocate internal data\n", me);
      exit(1);
    }
  }
}

/*
** _biffCheckKey()
**
** makes sure given key is kosher.  Panics and exit(1)s if given a NULL key
*/
void
_biffCheckKey(char *key) {
  char me[] = "_biffCheckKey";

  if (!key) {
    fprintf(stderr, "%s: PANIC: given NULL key\n", me);
    exit(1);
  }
  if (strlen(key) > BIFF_MAXKEYLEN) {
    fprintf(stderr, "%s: WARNING: truncating key \"%s\" to %d chars\n",
	    me, key, BIFF_MAXKEYLEN);
    key[BIFF_MAXKEYLEN] = 0;
  }
}

/*
** _biffFindKey()
**
** returns a pointer to the entry which contains the given key, or
** NULL if it was not found
*/
_biffEntry *
_biffFindKey(char *key) {
  int i = -1;
  _biffEntry *e;

  if (_biffNum) {
    for (i=0; i<=_biffNum-1; i++) {
      /* printf("HEY: comparing key[%d]=\"%s\" to \"%s\"\n", 
	 i, _biffErr[i]->key, key); */
      if (!strcmp(_biffErr[i]->key, key))
	break;
    }
    if (i == _biffNum) {
      i = -1;
    }
  }
  /* printf("HEY: index(\"%s\") = %d\n", key, i); */
  if (-1 == i) {
    e = NULL;
    _biffIdx = -1;
  }
  else {
    e = _biffErr[i];
    _biffIdx = i;
  }
  return(e);
}

/*
** _biffNewEntry()
**
** creates and initializes one new _biffEntry, returning a pointer to it
** panics and exit(1)s if there is a problem.
*/
_biffEntry *
_biffNewEntry(char *key) {
  char me[]="_biffInitEntry";
  _biffEntry *e;

  e = calloc(1, sizeof(_biffEntry));
  if (!e) {
    fprintf(stderr, "%s: couldn't make entry for new key \"%s\"\n", me, key);
    exit(1);
  }
  strcpy(e->key, key);
  e->AA = airArrayNew((void**)&(e->err), &(e->num), sizeof(char*), _BIFF_INCR);
  if (!e->AA) {
    fprintf(stderr, "%s: couldn't make array for new key \"%s\"\n", me, key);
    exit(1);
  }
  airArrayPointerCB(e->AA, NULL, free);
  return(e);
}

/*
** _biffNukeEntry()
**
** deletes given entry, and all info contained therein
*/
void
_biffNukeEntry(_biffEntry *e) {

  if (e) {
    airArraySetLen(e->AA, 0);
    airArrayNuke(e->AA);
    free(e);
  }
}

/*
** _biffAddKey()
**
** adds a key to _biffErr, and returns a pointer to the new entry
** assumes that given key does NOT appear in current list.
** panics and exit(1)s if there is a problem
*/
_biffEntry *
_biffAddKey(char *key) {
  char me[]="_biffAddKey";
  int i, newIdx;
  _biffEntry *e;

  /* find index of new key */
  for (i=0; i<=_biffNum-1; i++) {
    if (strcmp(key, _biffErr[i]->key) < 0) {
      /* we've hit the one which comes after the new key */
      break;
    }
  }
  /* if the for loop was never broken, _biffNum is the correct new index */
  newIdx = i;
  /* printf("HEY: index(new key \"%s\") = %d\n", key, i); */
  
  if (airArrayIncrLen(_biffAA, 1)) {
    fprintf(stderr, "%s: PANIC: couldn't accomodate one more key\n", me);
    exit(1);
  }

  /* _biffNum is now one bigger */
  for (i=_biffNum-2; i>=newIdx; i--) {
    _biffErr[i+1] = _biffErr[i];
  }
  e = _biffErr[newIdx] = _biffNewEntry(key);

  return(e);
}

/*
** _biffAddErr()
**
** adds a given message to the given entry.  The message is processed to
** convert all whitespace into ' ', and to eliminate whitespace at the
** end of the message.
** panics and exit(1)s if there is a problem
*/
void
_biffAddErr(_biffEntry *e, char *err) {
  char *buf, me[]="_biffAddErr";
  int i, len;

  /* printf("%s: HEY(before): err[%s]->num = %d\n", me, e->key, e->num); */
  if (airArrayIncrLen(e->AA, 1)) {
    fprintf(stderr, "%s: PANIC: couldn't add message for key %s\n",
	    me, e->key);
    exit(1);
  }
  /* printf("%s: HEY(after): err[%s]->num = %d\n", me, e->key, e->num); */
  buf = airStrdup(err);
  len = strlen(buf);
  for (i=0; i<=len-1; i++) {
    if (isspace(buf[i]))
      buf[i] = ' ';
  }
  i = len-1;
  while (isspace(buf[i])) {
    buf[i--] = 0;
  }
  /* printf("%s: HEY(after): err[%s]->num = %d\n", me, e->key, e->num); */
  /* printf("%s: HEY: err[%s][%d] now \"%s\"\n", me, e->key, e->num-1, buf); */
  e->err[e->num-1] = buf;
}

/***********************************************************************/
/***********************************************************************/

/*
******** biffSet()
**
** Sets given message "err" to be only message at "key".  "key" can be
** a new or existing key, but if it is an existing key, then existing
** messages at that key are lost
*/
void
biffSet(char *key, char *err) {
  _biffEntry *e;

  _biffInit();
  _biffCheckKey(key);

  e = _biffFindKey(key);
  if (!e) {
    /* not a key we remember seeing */
    e = _biffAddKey(key);
  }

  /* delete any existing messages at this index */
  airArraySetLen(e->AA, 0);

  /* add the new message */
  _biffAddErr(e, err);
}

/*
******** biffAdd()
**
** just like biffSet(), but doesn't delete existing messages
*/
void
biffAdd(char *key, char *err) {
  _biffEntry *e;

  _biffInit();
  _biffCheckKey(key);
  
  e = _biffFindKey(key);
  if (!e) {
    e = _biffAddKey(key);
  }

  /* add the new message */
  _biffAddErr(e, err);
}

/*
******** biffGet()
**
** creates a string which records all the errors at given key and
** returns it.  Returns NULL in case of error.  This function should
** be considered a glorified strdup(): it is the callers responsibility
** to free this string later
*/
char *
biffGet(char *key) {
  int i, max, len, sum;
  char me[] = "biffGet", *ret = NULL, *buf;
  _biffEntry *e;

  _biffInit();
  _biffCheckKey(key);

  /* find the index */
  e = _biffFindKey(key);
  if (!e) {
    /* error: not a key we remember seeing */
    fprintf(stderr, "%s: WARNING: no information for key \"%s\"\n", me, key);
    return(NULL);
  }
  if (!e->num) {
    /* there's a key, but no error messages.  Odd */
    return(airStrdup(""));
  }

  max = sum = 0;
  for (i=0; i<=e->num-1; i++) {
    len = strlen(e->err[i]) + strlen(e->key) + 4;
    sum += len;
    max = AIR_MAX(max, len);
  }
  buf = (char*)calloc(max+1, sizeof(char));
  ret = (char*)calloc(sum+1, sizeof(char));
  if (!(buf && ret)) {
    fprintf(stderr, "%s: PANIC: unable to allocate buffers\n", me);
    exit(1);
  }
  for (i=e->num-1; i>=0; i--) {
    sprintf(buf, "[%s] %s\n", key, e->err[i]);
    strcat(ret, buf);
  }
  free(buf);

  return(ret);
}

/*
******** biffDone()
**
** frees everything associated with given key, and shrinks list of keys
*/
void
biffDone(char *key) {
  char me[]="biffDone";
  int i, idx;
  _biffEntry *e;

  _biffInit();
  _biffCheckKey(key);

  e = _biffFindKey(key);
  if (!e) {
    fprintf(stderr, "%s: WARNING: no information for key \"%s\"\n", me, key);
    return;
  }
  idx = _biffIdx;

  _biffNukeEntry(e);
  for (i=idx; i<=_biffNum-2; i++) {
    _biffErr[i] = _biffErr[i+1];
  }
  airArrayIncrLen(_biffAA, -1);
}

void
biffMove(char *destKey, char *err, char *srcKey) {
  int i, len, max;
  char me[] = "biffMove", *buf;
  _biffEntry *dest, *src;

  _biffInit();
  _biffCheckKey(destKey);
  _biffCheckKey(srcKey);

  /* if srcKey and destKey are the same, this degenerates to biffAdd() */
  if (!strcmp(destKey, srcKey)) {
    biffAdd(srcKey, err);
    return;
  }

  dest = _biffFindKey(destKey);
  if (!dest) {
    dest = _biffAddKey(destKey);
  }
  src = _biffFindKey(srcKey);
  if (!src) {
    fprintf(stderr, "%s: WARNING: key \"%s\" unknown\n", me, srcKey);
    return;
  }

  max = 0;
  for (i=0; i<=src->num-1; i++) {
    len = strlen(src->err[i]) + strlen(src->key) + 4;
    max = AIR_MAX(max, len);
  }
  buf = (char*)calloc(max+1, sizeof(char));
  if (!buf) {
    fprintf(stderr, "%s: PANIC: can't allocate buffer\n", me);
    exit(1);
  }

  for (i=0; i<=src->num-1; i++) {
    sprintf(buf, "[%s] %s", srcKey, src->err[i]);
    /* printf("%s: HEY: moving \"%s\" to %s\n", me, buf, destKey); */
    _biffAddErr(dest, buf);
  }
  if (err) {
    _biffAddErr(dest, err);
  }
  biffDone(srcKey);
  
  free(buf);
}

char *
biffGetDone(char *key) {
  char *ret;

  _biffInit();
  _biffCheckKey(key);

  ret = biffGet(key);
  biffDone(key);

  return(ret);
}

