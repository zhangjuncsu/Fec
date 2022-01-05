#ifndef OVERLAPS_CHECK_H
#define OVERLAPS_CHECK_H

#include <array>
#include <vector>
#include <unordered_map>

#include "../common/packed_db.h"

struct PAF
{
	char strand;
	std::string tname, qname;
	int qid, qs, qe, qlen, tid, ts, te, tlen, rm, bl, mq;
};

void 
generate_partition_index_file_name(const char* m4_file_name, std::string& ret);

void
generate_partition_file_name(const char* m4_file_name, const index_t part, std::string& ret);

struct PartitionFileInfo
{
    std::string file_name;
    index_t min_seq_id;
    index_t max_seq_id;
};

void
load_partition_files_info(const char* idx_file_name, std::vector<PartitionFileInfo>& file_info_vec);

void partition_paf(PackedDB &reads,
				   const char *reads_file,
				   const char *paf_file,
				   const double min_cov_ration,
				   const index_t batch_size,
				   const int min_read_size,
				   int num_files,
				   const int min_cov,
				   const int threads,
				   std::vector<std::string> &names, 
				   std::unordered_map<std::string, std::array<int, 2>> &name2id, 
				   bool use_cache, 
				   bool so, 
				   const char* mmc, 
				   const char* fc);

#endif //OVERLAP_CHECK