# This variable will help provide a master list of all the sources.
# Add new source files here.
set(MEET_SOURCES
  enumall.c
  meetNrrd.c
  meetGage.c
  meetPull.c
  meet.h
  )

set(MEET_Dependencies hest
                      biff
                      nrrd
                      ell
                      unrrdu
                      gage
                      dye
                      limn
                      echo
                      hoover
                      seek
                      ten
                      pull
                      mite)

if(TEEM_BUILD_EXPERIMENTAL_LIBS)
  list(APPEND MEET_Dependencies alan
                                tijk
                                bane
                                elf
                                coil
                                push)
endif()

ADD_TEEM_LIBRARY(meet ${MEET_SOURCES}
                 DEPENDENCIES ${MEET_Dependencies})
