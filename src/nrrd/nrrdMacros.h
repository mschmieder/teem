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

#ifdef __cplusplus
extern "C" {
#endif


#ifndef NRRD_MACROS_HAS_BEEN_INCLUDED
#define NRRD_MACROS_HAS_BEEN_INCLUDED

/*
******** NRRD_AXIS_POS, NRRD_AXIS_IDX
**
** the guts of nrrdAxisPos() and nrrdAxisIdx(), for converting
** between "index space" location and "position" or "world space" location,
** given the centering, min and max "position", and number of samples.
**
** Unlike nrrdAxisPos() and nrrdAxisIdx(), this assumes that center
** is either nrrdCenterCell or nrrdCenterNode, but not nrrdCenterUnknown.
*/
#define NRRD_AXIS_POS(center, min, max, size, idx)    \
  (nrrdCenterCell == (center)                         \
   ? AIR_AFFINE(0, (idx) + 0.5, (size), (min), (max)) \
   : AIR_AFFINE(0, (idx), (size)-1, (min), (max)))

#define NRRD_AXIS_IDX(center, min, max, size, pos)    \
  (nrrdCenterCell == (center)                         \
   ? AIR_AFFINE((min), (pos), (max), 0, (size)) - 0.5 \
   : AIR_AFFINE((min), (pos), (max), 0, (size)-1))

/*
******** NRRD_AXIS_SPACING
**
** the guts of nrrdAxisSpacing(), determines the inter-sample
** spacing, given centering, min and max "position", and number of samples
**
** Unlike nrrdAxisSpacing, this assumes that center is either
** nrrdCenterCell or nrrdCenterNode, but not nrrdCenterUnknown.
*/
#define NRRD_AXIS_SPACING(center, min, max, size)  \
  (nrrdCenterCell == center                        \
   ? ((max) - (min))/(size)                        \
   : ((max) - (min))/((size) - 1))                 \

/*
******** NRRD_COORD_UPDATE
**
** This is for doing the "carrying" associated with gradually
** incrementing an array of coordinates.  Assuming that the given
** coordinate array "coord" has been incrementing by adding 1 to THE
** FIRST, THE ZERO-ETH, ELEMENT (this is a strong assumption), then,
** this macro is good for propagating the change up to higher axes
** (which really only happens when the position has stepped over the
** limit on a lower axis.)  Relies on the array of axes sizes "size",
** as as the length "dim" of "coord" and "size".
**
** This may be turned into something more general purpose soon. 
*/
#define NRRD_COORD_UPDATE(coord, size, dim)    \
do {                                           \
  int d;                                       \
  for (d=0;                                    \
       d < (dim)-1 && (coord)[d] == (size)[d]; \
       d++) {                                  \
    (coord)[d] = 0;                            \
    (coord)[d+1]++;                            \
  }                                            \
} while (0)

/*
******** NRRD_COORD_INCR
**
** same as NRRD_COORD_UPDATE, but starts by incrementing coord[idx]
*/
#define NRRD_COORD_INCR(coord, size, dim, idx) \
do {                                           \
  int d;                                       \
  for (d=idx, (coord)[d]++;                    \
       d < (dim)-1 && (coord)[d] == (size)[d]; \
       d++) {                                  \
    (coord)[d] = 0;                            \
    (coord)[d+1]++;                            \
  }                                            \
} while (0)

/*
******** NRRD_COORD_INDEX
**
** Given a coordinate array "coord", as well as the array sizes "size"
** and dimension "dim", calculates the linear index, and stores it in
** I.
*/
#define NRRD_COORD_INDEX(I, coord, size, dim) \
do {                                          \
  int d;                                      \
  for (d=(dim)-1, (I)=(coord)[d--];           \
       d >= 0;                                \
       d--) {                                 \
    (I) = (coord)[d] + (size)[d]*(I);         \
  }                                           \
} while (0)

#define NRRD_COORD_GEN(coord, size, dim, I)   \
do {                                          \
  int d;                                      \
  for (d=0; d<=(dim)-1; d++) {                \
    (coord)[d] = I % (size)[d];               \
    I /= (size)[d];                           \
  }                                           \
} while (0)

#endif /* NRRD_MACROS_HAS_BEEN_INCLUDED */

#ifdef __cplusplus
}
#endif
