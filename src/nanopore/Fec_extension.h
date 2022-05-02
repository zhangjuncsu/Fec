#pragma once

#include "Fec_aligner.h"
#include "Fec_daligner.h"
#include "edlib_wrapper.h"
#include "common/alignment.h"

BOOL
cns_extension(OcAlignData* align_data,
			  OcDalignData* dalign_data,
			  FullEdlibAlignData* falign_data,
			  const char* query,
              const int query_start, 
              const int query_size, 
			  const char* target,
              const int target_start, 
              const int target_size, 
			  int min_align_size,
			  M5Record& m5, 
			  bool recuse_long_indel = false);