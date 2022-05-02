#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include "../common/defs.h"

struct ConsensusOptions
{
    int         input_type;
    const char* m4;
    const char* reads;
    const char* corrected_reads;
    int         num_threads;
    index_t     batch_size;
    double      min_mapping_ratio;
    int         min_align_size;
    int         min_cov;
    index_t     min_size;
    bool        print_usage_info;
    int         tech;
	int			num_partition_files;
    int    cache;
    int     second_overlapping;
    bool resuce_long_indel;
    bool full;
    const char* minimap2_filter;
    const char* overlap_filter;
    std::string minimap2_command;
};

void
Usage(const char* prog);

int
parse_arguments(int argc, char* argv[], ConsensusOptions& t);

void
print_options(ConsensusOptions& t);

typedef ConsensusOptions ReadsCorrectionOptions;

#endif // OPTIONS_H
