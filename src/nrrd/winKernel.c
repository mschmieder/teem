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

#define _SINC(x) (sin(M_PI*x)/(M_PI*x))

double
_nrrdWindSincInt(double *parm) {

  /* This isn't true, but there aren't good accurate, closed-form
     approximations for these integrals ... */
  return 1.0;
}

double
_nrrdDWindSincInt(double *parm) {

  /* ... or their derivatives */
  return 0.0;
}

double
_nrrdWindSincSup(double *parm) {
  double S;

  S = parm[0];
  return parm[1]*S;
}

#define POW1(S) (S)
#define POW2(S) (S*S)
#define POW3(S) (S*S*S)

#define WS_1_F(name, mac, spow)              \
float                                        \
_nrrd##name##_1_f(float x, double *parm) {   \
  float R, S;                                \
                                             \
  S = parm[0]; R = parm[1];                  \
  x /= S;                                    \
  return mac(x, R)/spow(S);                  \
}

#define WS_N_F(name, mac, spow)                                     \
void                                                                \
_nrrd##name##_N_f(float *f, float *x, size_t len, double *parm) {   \
  float S, R, t;                                                    \
  size_t i;                                                         \
                                                                    \
  S = parm[0]; R = parm[1];                                         \
  for (i=0; i<len; i++) {                                           \
    t = x[i]/S;                                                     \
    f[i] = mac(t, R)/spow(S);                                       \
  }                                                                 \
}
#define WS_1_D(name, mac, spow)             \
double                                      \
_nrrd##name##_1_d(double x, double *parm) { \
  double R, S;                              \
                                            \
  S = parm[0]; R = parm[1];                 \
  x /= S;                                   \
  return mac(x, R)/spow(S);                 \
}

#define WS_N_D(name, mac, spow)                                     \
void                                                                \
_nrrd##name##_N_d(double *f, double *x, size_t len, double *parm) { \
  double S, R, t;                                                   \
  size_t i;                                                         \
                                                                    \
  S = parm[0]; R = parm[1];                                         \
  for (i=0; i<len; i++) {                                           \
    t = x[i]/S;                                                     \
    f[i] = mac(t, R)/spow(S);                                       \
  }                                                                 \
}

/* ------------------------------------------------------------ */

#define _HANN(x, R) \
   (x > R ? 0 : (x < -R ? 0 : (x == 0 ? 1.0 : \
    (1 + cos(M_PI*x/R))*_SINC(x)/2)))

WS_1_D(Hann, _HANN, POW1)
WS_1_F(Hann, _HANN, POW1)
WS_N_F(Hann, _HANN, POW1)
WS_N_D(Hann, _HANN, POW1)

NrrdKernel
_nrrdKernelHann = {
  "hann",
  2, _nrrdWindSincSup,  _nrrdWindSincInt,   
  _nrrdHann_1_f, _nrrdHann_N_f, _nrrdHann_1_d, _nrrdHann_N_d
};
NrrdKernel *
nrrdKernelHann = &_nrrdKernelHann;

/* ------------------------------------------------------------ */

#define _DHANN(x, R)                                            \
   (x > R ? 0.0 : (x < -R ? 0.0 : (x == 0 ? 0.0 :              \
    (R*(1 + cos(M_PI*x/R))*(M_PI*x*cos(M_PI*x) - sin(M_PI*x))  \
       - M_PI*x*sin(M_PI*x)*sin(M_PI*x/R))/(2*R*M_PI*x*x))))

WS_1_D(DHann, _DHANN, POW2)
WS_1_F(DHann, _DHANN, POW2)
WS_N_F(DHann, _DHANN, POW2)
WS_N_D(DHann, _DHANN, POW2)

NrrdKernel
_nrrdKernelDHann = {
  "hannD",
  2, _nrrdWindSincSup, _nrrdDWindSincInt,  
  _nrrdDHann_1_f,  _nrrdDHann_N_f,  _nrrdDHann_1_d,  _nrrdDHann_N_d
};
NrrdKernel *
nrrdKernelHannD = &_nrrdKernelDHann;

/* ------------------------------------------------------------ */

#define _DDHANN_A(x, R) \
  (2*M_PI*R*cos(M_PI*x)*(R + R*cos(M_PI*x/R) + M_PI*x*sin(M_PI*x/R)))
#define _DDHANN_B(x, R)                                      \
  (cos(M_PI*x/R)*(M_PI*M_PI*x*x + R*R*(M_PI*M_PI*x*x - 2)) + \
   R*(R*(M_PI*M_PI*x*x - 2) - 2*M_PI*x*sin(M_PI*x/R)))
#define _DDHANN(x, R)                                                     \
   (x > R ? 0 : (x < -R ? 0 : (x == 0 ? -M_PI*M_PI*(3 + 2*R*R)/(6*R*R) :  \
    -(_DDHANN_A(x,R) + sin(M_PI*x)*_DDHANN_B(x,R)/x)/(2*M_PI*R*R*x*x)       \
    )))

WS_1_D(DDHann, _DDHANN, POW3)
WS_1_F(DDHann, _DDHANN, POW3)
WS_N_F(DDHann, _DDHANN, POW3)
WS_N_D(DDHann, _DDHANN, POW3)

NrrdKernel
_nrrdKernelDDHann = {
  "hannDD",
  2, _nrrdWindSincSup, _nrrdDWindSincInt,  
  _nrrdDDHann_1_f, _nrrdDDHann_N_f, _nrrdDDHann_1_d, _nrrdDDHann_N_d
};
NrrdKernel *
nrrdKernelHannDD = &_nrrdKernelDDHann;

/* ------------------------------------------------------------ */

#define _BLACK(x, R)                                            \
   (x > R ? 0 : (x < -R ? 0 : (x == 0 ? 1.0 :                   \
    (0.42 + cos(M_PI*x/R)/2 + 0.08*cos(2*M_PI*x/R))*_SINC(x))))

WS_1_D(Black, _BLACK, POW1)
WS_1_F(Black, _BLACK, POW1)
WS_N_F(Black, _BLACK, POW1)
WS_N_D(Black, _BLACK, POW1)

NrrdKernel
_nrrdKernelBlackman = {
  "blackman",
  2, _nrrdWindSincSup,  _nrrdWindSincInt,   
  _nrrdBlack_1_f, _nrrdBlack_N_f, _nrrdBlack_1_d, _nrrdBlack_N_d
};
NrrdKernel *
nrrdKernelBlackman = &_nrrdKernelBlackman;

/* ------------------------------------------------------------ */

#define _DBLACK_A(x, R)                                   \
  R*x*cos(M_PI*x)*(2.638937829015426 + M_PI*cos(M_PI*x/R) \
                   + 0.5026548245743669*cos(2*M_PI*x/R))
#define _DBLACK_B(x, R)                                                     \
  sin(M_PI*x)*(-0.84*R - R*cos(M_PI*x/R) - 0.16*R*cos(2*M_PI*x/R) -         \
               M_PI*x*sin(M_PI*x/R) - 1.0053096491487339*x*sin(2*M_PI*x/R))
#define _DBLACK(x, R) (_DBLACK_A(x,R) + _DBLACK_B(x,R))/(2*M_PI*R*x*x)

WS_1_D(DBlack, _DBLACK, POW2)
WS_1_F(DBlack, _DBLACK, POW2)
WS_N_F(DBlack, _DBLACK, POW2)
WS_N_D(DBlack, _DBLACK, POW2)

NrrdKernel
_nrrdKernelDBlack = {
  "blackmanD",
  2, _nrrdWindSincSup, _nrrdDWindSincInt,  
  _nrrdDBlack_1_f,  _nrrdDBlack_N_f,  _nrrdDBlack_1_d,  _nrrdDBlack_N_d
};
NrrdKernel *
nrrdKernelBlackmanD = &_nrrdKernelDBlack;

/* ------------------------------------------------------------ */

#define _DDBLACK(x, R)                                                              \
  ((R*x*cos(M_PI*x)*(-2.638937829015426*R - M_PI*R*cos((M_PI*x)/R)                  \
		    - 0.5026548245743669*R*cos((2*M_PI*x)/R)                        \
		    - M_PI*M_PI*x*sin((M_PI*x)/R)                                   \
		    - 3.158273408348595*x*sin((2*M_PI*x)/R))                        \
   + sin(M_PI*x)*((-4.934802200544679*x*x                                           \
		   + R*R*(1 - 4.934802200544679*x*x))*cos((M_PI*x)/R)               \
		  + (-3.158273408348595*x*x                                         \
		     + R*R*(0.16 - 0.7895683520871487*x*x))*cos((2*M_PI*x)/R)       \
		  + R*(0.84*R - 4.14523384845753*R*x*x                              \
		       + M_PI*x*sin((M_PI*x)/R)                                     \
		       + 1.0053096491487339*x*sin((2*M_PI*x)/R))))/(M_PI*R*R*x*x*x))

WS_1_D(DDBlack, _DDBLACK, POW3)
WS_1_F(DDBlack, _DDBLACK, POW3)
WS_N_F(DDBlack, _DDBLACK, POW3)
WS_N_D(DDBlack, _DDBLACK, POW3)

NrrdKernel
_nrrdKernelDDBlack = {
  "blackDD",
  2, _nrrdWindSincSup, _nrrdDWindSincInt,  
  _nrrdDDBlack_1_f, _nrrdDDBlack_N_f, _nrrdDDBlack_1_d, _nrrdDDBlack_N_d
};
NrrdKernel *
nrrdKernelBlackmanDD = &_nrrdKernelDDBlack;
