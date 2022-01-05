ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ./$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ./$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libfec.a

SOURCES      := common/alignment.cpp \
		common/buffer_line_iterator.cpp \
		common/defs.cpp \
		common/diff_gapalign.cpp \
		common/fasta_reader.cpp \
		common/gapalign.cpp \
		common/lookup_table.cpp \
		common/packed_db.cpp \
		common/sequence.cpp \
		common/split_database.cpp \
		common/xdrop_gapalign.cpp

SRC_INCDIRS  := common \

SUBMAKEFILES := src/fec.mk
