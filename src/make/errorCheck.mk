#
# teem: Gordon Kindlmann's research software
# Copyright (C) 2002, 2001, 2000, 1999, 1998 University of Utah
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#

####
#### errorCheck.mk: checks on validity of the variables that make uses
#### use

# all the architectures currently supported
KNOWN_ARCH = irix6.n32 irix6.64 linux cygwin solaris darwin aix

# there is no default architecture
checkArchSet = $(if $(TEEM_ARCH),,\
$(warning *)\
$(warning *)\
$(warning * Environment variable TEEM_ARCH not set.)\
$(warning * Possible settings currently supported:)\
$(warning * $(KNOWN_ARCH))\
$(warning *)\
$(warning *)\
$(error Make quitting))

# verify that we can recognize the architecture setting
checkArchValid = $(if $(strip $(findstring $(TEEM_ARCH),$(KNOWN_ARCH))),,\
$(warning *)\
$(warning *)\
$(warning * Environment variable TEEM_ARCH="$(TEEM_ARCH)" unknown)\
$(warning * Possible settings currently supported:)\
$(warning * $(KNOWN_ARCH))\
$(warning *)\
$(warning *)\
$(error Make quitting))

## Other error checking ...
##
checkShext = $(if $(filter undefined,$(origin TEEM_LINK_SHARED)),,\
$(if $(TEEM_SHEXT),,\
$(warning *)\
$(warning *)\
$(warning * Cannot do shared library linking with TEEM_SHEXT unset)\
$(warning * Set it in teem/make/$(ARCH).mk or as environment var)\
$(warning *)\
$(warning *)\
$(error Make quitting)))

checkPurify = $(if $(filter undefined,$(origin TEEM_PURIFY)),,\
$(if $(TEEM_PURIFY_CMD),,\
$(warning *)\
$(warning *)\
$(warning * Purification requested, but TEEM_PURIFY_CMD not set)\
$(warning * Set it in teem/make/$(ARCH).mk or as environment var)\
$(warning *)\
$(warning *)\
$(error Make quitting)))

