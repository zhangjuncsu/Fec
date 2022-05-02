ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := Fec
SOURCES  := main.cpp \
	dw.cpp \
	FEC_AlnGraphBoost.C \
	fec_correction.cpp \
	options.cpp \
	reads_correction_aux.cpp \
	reads_correction_can.cpp \
	overlaps_check.cpp \
	nanopore/Fec_extension.cpp \
	nanopore/edlib_wrapper.c \
	nanopore/edlib_ex.c \
	nanopore/edlib_ex_aux.c \
	nanopore/Fec_aligner.c \
	nanopore/align.c \
	nanopore/Fec_daligner.c \
	nanopore/edlib.cpp \
	nanopore/kalloc.c \
	nanopore/kstring.c \
	nanopore/util.c 

SRC_INCDIRS  := . libboost

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lfec
TGT_PREREQS := libfec.a

SUBMAKEFILES :=
