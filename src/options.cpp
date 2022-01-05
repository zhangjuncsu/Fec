#include "options.h"

#include <cstring>
#include <unistd.h>

#include <iostream>

using namespace std;

static int num_threads = 1;
static index_t batch_size = 100000;
static double mapping_ratio = 0.6;
static int align_size = 1000;
static int cov = 4;
static int min_size = 2000;
static bool print_usage = false;
static int num_partition_files = 100;
static char* overlap_filter = "0.6";
static int cache = 1;
static int second_overlapping = 1;

int MINIMAP2_k = 19;
int MINIMAP2_w = 5;
float MINIMAP2_f = 0.0002;

static const char num_threads_n   = 't';
static const char batch_size_n    = 'p';
static const char mapping_ratio_n = 'r';
static const char align_size_n    = 'a';
static const char cov_n           = 'c';
static const char min_size_n      = 'l';
static const char usage_n         = 'h';
static const char num_partition_files_n = 'k';
static const char cache_n		='e';
static const char second_overlapping_n = 's';
static const char minimap_filter_n = 'm';
static const char overlap_filter_n = 'f';
static const char kmer_n = 'K';
static const char window_n = 'w';
static const char hpc_n = 'H';

void print_default_options() {
	cerr << '-' << num_threads_n << ' ' << num_threads 
		 << ' '
		 << '-' << batch_size_n << ' ' << batch_size 
		 << ' '
		 << '-' << mapping_ratio_n << ' ' << mapping_ratio 
		 << ' '
		 << '-' << align_size_n << ' ' << align_size 
		 << ' '
		 << '-' << cov_n << ' ' << cov 
		 << ' '
		 << '-' << min_size_n << ' ' << min_size 
		 << ' '
		 << '-' << num_partition_files_n << ' ' << num_partition_files 
		 << ' '
		 << '-' << cache_n << ' ' << cache 
		 << ' '
		 << '-' << second_overlapping_n << ' ' << second_overlapping 
		 << ' '
		 << '-' << minimap_filter_n << ' ' << MINIMAP2_f 
		 << ' '
		 << '-' << overlap_filter_n << ' ' << overlap_filter 
		 << ' ' 
		 << '-' << kmer_n << ' ' << MINIMAP2_k 
		 << ' ' 
		 << '-' << window_n << ' ' << MINIMAP2_w 
		 << ' ' 
		 << '-' << hpc_n 
		 <<"\n";
}

void Usage(const char* prog) {
	cerr << "USAGE:\n"
		 << prog << ' ' << "[options]" << ' ' << "input" << ' ' << "reads" << ' ' << "output" << "\n";
	cerr << "\n" << "OPTIONS:" << "\n";
	
	cerr << "-" << num_threads_n << " <Integer>\t" << "number of threads (CPUs)" << "\n";
	
	cerr << "-" << batch_size_n << " <Integer>\t" << "batch size that the reads will be partitioned" << "\n";
	
	cerr << "-" << mapping_ratio_n << " <Real>\t" << "minimum mapping ratio" << "\n";
	
	cerr << "-" << align_size_n << " <Integer>\t" << "minimum overlap size" << "\n";
	
	cerr << "-" << cov_n << " <Integer>\t" << "minimum coverage under consideration" << "\n";
	
	cerr << "-" << min_size_n << " <Integer>\t" << "minimum length of corrected sequence" << "\n";
	
	cerr << "-" << num_partition_files_n << " <Integer>\t" 
		 << "number of partition files when partitioning overlap results" 
		 << " (if < 0, then it will be set to system limit value)"
		 << "\n";

	cerr << "-" << cache_n << " <Integer>\t" << "use cache or not: 0 not use, 1 use" << "\n";

	cerr << "-" << second_overlapping_n << " <Integer>\t" << "perform second overlapping: 0 no, 1 yes" << "\n";
	
	cerr << "-" << minimap_filter_n << " <String>\t" << "the second round minimap command" <<"\n";

	cerr << "-" << overlap_filter_n << " <String>\t" << "command to filter overlaps" <<"\n";

	cerr << "-" << kmer_n << " <Integer>\t" << "kmer length used in minimap2" <<"\n";

	cerr << "-" << window_n << " <Integer>\t" << "window used in minimap2" <<"\n";

	cerr << "-" << hpc_n << "\t\thomopolymen compress used in minimap2" <<"\n";

	cerr << "-" << usage_n << "\t\t" << "print usage info." << "\n";
	
	cerr << "\nDefault Options:\n";
	print_default_options();
}

ConsensusOptions init_consensus_options() {
	ConsensusOptions t;
	t.m4                    = NULL;
	t.reads                 = NULL;
	t.corrected_reads       = NULL;
	t.num_threads           = num_threads;
	t.batch_size            = batch_size;
	t.min_mapping_ratio     = mapping_ratio;
	t.min_align_size        = align_size;
	t.min_cov               = cov;
	t.min_size              = min_size;
	t.print_usage_info      = print_usage;
	t.num_partition_files 	= num_partition_files;
	t.cache				= cache;
	t.second_overlapping = second_overlapping;
	t.overlap_filter = overlap_filter;
	t.minimap2_command = "minimap2 ";
    return t;
}

int
parse_arguments(int argc, char* argv[], ConsensusOptions& t)
{
	bool parse_success = true;
	t = init_consensus_options();
	
	int opt_char;
    char err_char;
    opterr = 0;
	while((opt_char = getopt(argc, argv, "i:e:s:t:p:r:a:c:l:k:m:K:w:f:hH")) != -1) {
		switch (opt_char) {
			case cache_n:
				if(optarg[0] == '0')
					t.cache = 0;
				else if(optarg[0] == '1')
					t.cache = 1;
				else{
					fprintf(stderr, "invalid argument to option '%c': %s\n", cache_n, optarg);
					return 1;
				}
				break;
			case second_overlapping_n:
				if(optarg[0] == '0')
					t.second_overlapping = 0;
				else if(optarg[0] == '1')
					t.second_overlapping = 1;
				else{
					fprintf(stderr, "invalid argument to option '%c': %s\n", second_overlapping_n, optarg);
					return 1;
				}
				break;
			case num_threads_n:
				t.num_threads = atoi(optarg);
				break;
			case minimap_filter_n:
				MINIMAP2_f = atof(optarg);
				break;
			case overlap_filter_n:
				t.overlap_filter = optarg;
				break;
			case batch_size_n:
				t.batch_size = atoll(optarg);
				break;
			case mapping_ratio_n:
				t.min_mapping_ratio = atof(optarg);
				break;
			case align_size_n:
				t.min_align_size = atoi(optarg);
				break;
			case cov_n:
				t.min_cov = atoi(optarg);
				break;
			case min_size_n:
				t.min_size = atoll(optarg);
				break;
			case usage_n:
				t.print_usage_info = true;
				break;
			case num_partition_files_n:
				t.num_partition_files = atoi(optarg);
				break;
			case kmer_n:
				MINIMAP2_k = atoi(optarg);
				break;
			case window_n:
				MINIMAP2_w = atoi(optarg);
				break;
			case hpc_n:
				t.minimap2_command += "-H ";
				break;
			case '?':
                err_char = (char)optopt;
				fprintf(stderr, "unrecognised option '%c'\n", err_char);
                return 1;
                break;
            case ':':
                err_char = (char)optopt;
				fprintf(stderr, "argument to option '%c' is missing.\n", err_char);
                return 1;
                break;
		}
	}
	
	if (t.num_threads <= 0)
	{
		std::cerr << "cpu threads must be greater than 0\n";
		parse_success = false;
	}
	if (t.batch_size <= 0)
	{
		std::cerr << "batch size must be greater than 0\n";
		parse_success = false;
	}
	if (t.min_mapping_ratio < 0.0)
	{
		std::cerr << "mapping ratio must be >= 0.0\n";
		parse_success = false;
	}
	if (t.min_cov < 0)
	{
		std::cerr << "coverage must be >= 0\n";
		parse_success = false;
	}
	if (t.min_cov < 0)
	{
		std::cerr << "sequence size must be >= 0\n";
		parse_success = false;
	}
	
	if (argc < 3) return 1;
	
	t.m4 = argv[argc - 3];
	t.reads = argv[argc - 2];
	t.corrected_reads = argv[argc - 1];
	t.minimap2_command += "-k " + std::to_string(MINIMAP2_k) + " -w " + std::to_string(MINIMAP2_w) + " -t " + std::to_string(t.num_threads);
	t.minimap2_command += " -f " + std::to_string(MINIMAP2_f);
	
	if (parse_success) return 0;
	return 1;
}

void
print_options(ConsensusOptions& t)
{
	if (t.m4) cout << "reads\t" << t.m4 << "\n";
	if (t.reads) cout << "output\t" << t.reads << "\n";
	if (t.corrected_reads) cout << "m4\t" << t.corrected_reads << "\n";
	cout << "number of threads:\t" << t.num_threads << "\n";
	cout << "batch size:\t" << t.batch_size << "\n";
	cout << "mapping ratio:\t" << t.min_mapping_ratio << "\n";
	cout << "align size:\t" << t.min_align_size << "\n";
	cout << "cov:\t" << t.min_cov << "\n";
	cout << "min size:\t" << t.min_size << "\n";
	cout << "partition files:\t" << t.num_partition_files << "\n";
	cout << "use cache:\t" << t.cache << "\n";
	cout << "perform second overlapping:\t" << t.second_overlapping << "\n";
	cout << "filter command:\t" << t.overlap_filter << "\n";
	cout << "minimap2 command:\t" << t.minimap2_command << "\n";
}
