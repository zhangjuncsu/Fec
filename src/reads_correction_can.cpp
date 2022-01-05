#include "reads_correction_can.h"

#include <cstring>

#include <array>
#include <atomic>
#include "../common/fasta_reader.h"

#include "FEC_AlnGraphBoost.H"
#include "fec_correction.h"
#include "overlaps_check.h"
#include "overlaps_store.h"
#include "overlaps_check.h"
#include <iostream>

#define CACHE_SIZE 2000

using namespace std;

std::atomic<size_t>INDEX {0};
std::unordered_map<idx_t, std::array<idx_t, 2>>OFFSETS;
std::vector<idx_t>GROUP_IDS;
std::vector<size_t>GROUP_LINE;

size_t GET_ID(std::vector<idx_t> &ids){
	auto curr = INDEX.fetch_add(1);
	if(curr < GROUP_LINE.size() - 1){
		assert(GROUP_LINE[curr + 1] - GROUP_LINE[curr] <= ids.size());
		for(size_t i = GROUP_LINE[curr]; i < GROUP_LINE[curr + 1]; ++i){
			ids[i - GROUP_LINE[curr]] = GROUP_IDS[i];
		}
		return GROUP_LINE[curr + 1] - GROUP_LINE[curr];
	}else{
		return (size_t)0;
	}
}

struct CmpExtensionCandidateBySid
{
	bool operator()(const ExtensionCandidate& a, const ExtensionCandidate& b)
	{
		return a.sid < b.sid;
	}
};

void*
reads_correction_func_can_cache(void* arg)
{
    ConsensusThreadData& cns_data = *static_cast<ConsensusThreadData*>(arg);
	ExtensionCandidate* candidates = cns_data.candidates;
	const index_t num_candidates = cns_data.num_candidates;
    index_t i = 0, j;
	std::vector<idx_t>ids (CACHE_SIZE);
	size_t size;
	while((size = GET_ID(ids)) != 0){
		cns_data.cache_.SetIds(ids);
		for(size_t k = 0; k < size; ++k){
			idx_t id = ids[k];
			auto iter_off = OFFSETS.find(id);
			assert(iter_off != OFFSETS.end());
			idx_t b = iter_off->second[0];
			idx_t e = iter_off->second[1];
			ns_meap_cns::consensus_one_read_can(&cns_data, id, b, e);
			if(cns_data.cns_results.size() >= MAX_CNS_RESULTS){
				pthread_mutex_lock(&cns_data.out_lock);
				for(std::vector<CnsResult>::iterator iter = cns_data.cns_results.begin(); iter != cns_data.cns_results.end(); ++iter){
					(*cns_data.out) << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n";
					std::string& seq = iter->seq;
					(*cns_data.out) << seq << "\n";
				}
				cns_data.cns_results.clear();
				pthread_mutex_unlock(&cns_data.out_lock);
			}
		}
		cns_data.cache_.Clear();
	}
	return NULL;
}

void*
reads_correction_func_can(void* arg)
{
    ConsensusThreadData& cns_data = *static_cast<ConsensusThreadData*>(arg);
	ExtensionCandidate* candidates = cns_data.candidates;
	const index_t num_candidates = cns_data.num_candidates;
    index_t i = 0, j;
    while (i < num_candidates)
    {
        const index_t sid = candidates[i].sid;
        j = i + 1;
        while (j < num_candidates && candidates[j].sid == sid) ++j;
        if (j - i < cns_data.rco.min_cov) { i = j; continue; }
        if (candidates[i].ssize < cns_data.rco.min_size * 0.95) { i = j; continue; }
		ns_meap_cns::consensus_one_read_can(&cns_data, sid, i, j);
		if (cns_data.cns_results.size() >= MAX_CNS_RESULTS)
		{
			pthread_mutex_lock(&cns_data.out_lock);
			for (std::vector<CnsResult>::iterator iter = cns_data.cns_results.begin(); iter != cns_data.cns_results.end(); ++iter)
			{
				(*cns_data.out) << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n";
				std::string& seq = iter->seq;
				(*cns_data.out) << seq << "\n";
			}
			cns_data.cns_results.clear();
			pthread_mutex_unlock(&cns_data.out_lock);
		}
        i = j;
    }
    return NULL;
}

void
consensus_one_partition_can_cache(const char* m4_file_name,
						const index_t min_read_id,
						const index_t max_read_id,
						ReadsCorrectionOptions& rco,
						PackedDB& reads,
						std::ostream& out, 
						bool cache)
{
	idx_t num_ec, i = 0, j;
	std::unordered_map<idx_t, bool>done;
	std::unordered_map<idx_t, std::array<idx_t, 2>>offsets;
	std::vector<idx_t>group_ids;
	std::vector<size_t>group_line;
	group_line.push_back(0);
	ExtensionCandidate* ec_list = load_partition_data<ExtensionCandidate>(m4_file_name, num_ec);
    ConsensusThreadData* pctds[rco.num_threads];
	std::sort(ec_list, ec_list + num_ec, CmpExtensionCandidateBySid());
	while(i < num_ec){
		idx_t sid = ec_list[i].sid;
		j = i + 1;
		while (j < num_ec && ec_list[j].sid == sid) ++j;
		if(j - i < rco.min_cov) { i = j; continue; }
		if(ec_list[i].ssize < rco.min_size * 0.95) { i = j; continue; }
		done[sid] = false;
		offsets[sid] = {i, j};
		i = j;
	}
	for(auto id = done.begin(); id != done.end(); ++id){
		if(!id->second){
			group_ids.push_back(id->first);
			id->second = true;
			auto iter_off = offsets.find(id->first);
			assert(iter_off != offsets.end());
			for(i = iter_off->second[0]; i < iter_off->second[1]; ++i){
				if(ec_list[i].sid == iter_off->first){
					auto iter_done = done.find(ec_list[i].qid);
					if(iter_done != done.end() && !iter_done->second){
						group_ids.push_back(ec_list[i].qid);
						iter_done->second = true;
						if(group_ids.size() - group_line.back() >= CACHE_SIZE) break;
					}
				}
			}
			group_line.push_back(group_ids.size());
		}
	}
	if(group_ids.size() > group_line.back()){
		group_line.push_back(group_ids.size());
	}
	for(int tid = 0; tid < rco.num_threads; ++tid){
		pctds[tid] = new ConsensusThreadData(&rco, tid, &reads, ec_list, 0, &out);
		pctds[tid]->useCache = cache;
	}
	OFFSETS = offsets;
	GROUP_LINE = group_line;
	GROUP_IDS = group_ids;
    pthread_t thread_ids[rco.num_threads];
    for (int i = 0; i < rco.num_threads; ++i)
        pthread_create(&thread_ids[i], NULL, reads_correction_func_can_cache, static_cast<void*>(pctds[i]));
    for (int i = 0; i < rco.num_threads; ++i)
        pthread_join(thread_ids[i], NULL);
	for (int i = 0; i < rco.num_threads; ++i)
	{
		std::vector<CnsResult>& cns_results = pctds[i]->cns_results;
		for (std::vector<CnsResult>::iterator iter = cns_results.begin(); iter != cns_results.end(); ++iter)
		{
			out << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n";
			std::string& seq = iter->seq;
			out << seq << "\n";
		}
	}

    delete[] ec_list;
    for (int i = 0; i < rco.num_threads; ++i) delete pctds[i];
}

void
consensus_one_partition_can(const char* m4_file_name,
						const index_t min_read_id,
						const index_t max_read_id,
						ReadsCorrectionOptions& rco,
						PackedDB& reads,
						std::ostream& out)
{
	idx_t num_ec;
	ExtensionCandidate* ec_list = load_partition_data<ExtensionCandidate>(m4_file_name, num_ec);
    ConsensusThreadData* pctds[rco.num_threads];
	build_cns_thrd_data_can(ec_list, num_ec, min_read_id, max_read_id, &rco, &reads, &out, pctds);
    for(int i = 0; i < rco.num_threads; ++i) pctds[i]->useCache = false;
	pthread_t thread_ids[rco.num_threads];
    for (int i = 0; i < rco.num_threads; ++i)
        pthread_create(&thread_ids[i], NULL, reads_correction_func_can, static_cast<void*>(pctds[i]));
    for (int i = 0; i < rco.num_threads; ++i)
        pthread_join(thread_ids[i], NULL);
	for (int i = 0; i < rco.num_threads; ++i)
	{
		std::vector<CnsResult>& cns_results = pctds[i]->cns_results;
		for (std::vector<CnsResult>::iterator iter = cns_results.begin(); iter != cns_results.end(); ++iter)
		{
			out << ">" << iter->id << "_" << iter->range[0] << "_" << iter->range[1] << "_" << iter->seq.size() << "\n";
			std::string& seq = iter->seq;
			out << seq << "\n";
		}
	}

    delete[] ec_list;
    for (int i = 0; i < rco.num_threads; ++i) delete pctds[i];
}

static void
partition_can_cache(PackedDB& reads, const char* reads_file, const char* m4_path, const int batch_size, const int min_size, const int num_partition_files, const int min_cov, const double min_cov_ratio, const int threads, std::vector<std::string> &names, std::unordered_map<std::string, std::array<int, 2>>& name2id, bool use_cache, bool so, const char* mmc, const char* of)
{
	char job_finished[2048];
	sprintf(job_finished, "%s.partition_finished", m4_path);
	if (access(job_finished, F_OK) == 0) {
		fprintf(stderr, "Partition Candidates Has Been Finished, Skipt It.\n");
		return;
	}
	partition_paf(reads, reads_file, m4_path, min_cov_ratio, batch_size, min_size, num_partition_files, min_cov, threads, names, name2id, use_cache, so, mmc, of);
	FILE* out = fopen(job_finished, "w");
	fclose(out);
}

int reads_correction_can(ReadsCorrectionOptions& rco)
{
	std::vector<std::string>names;
	std::unordered_map<std::string, std::array<int, 2>>name2id;
	PackedDB reads;
	reads.load_fasta_db(rco.reads, MAX_SEQ_SIZE, names, name2id);
	partition_can_cache(reads, rco.reads, rco.m4, rco.batch_size, rco.min_size, rco.num_partition_files, rco.min_cov, rco.min_mapping_ratio, rco.num_threads, names ,name2id, rco.cache, rco.second_overlapping, rco.minimap2_command.c_str(), rco.overlap_filter);
	std::string idx_file_name;
	generate_partition_index_file_name(rco.m4, idx_file_name);
	std::vector<PartitionFileInfo> partition_file_vec;
	load_partition_files_info(idx_file_name.c_str(), partition_file_vec);
	std::ofstream out;
	open_fstream(out, rco.corrected_reads, std::ios::out);
	char process_info[2048];
	char job_finished[2048];
	for (std::vector<PartitionFileInfo>::iterator iter = partition_file_vec.begin(); iter != partition_file_vec.end(); ++iter)
	{
		sprintf(job_finished, "%s.correction_finished", iter->file_name.c_str());
		sprintf(process_info, "processing %s", iter->file_name.c_str());
		DynamicTimer dtimer(process_info);
		INDEX.fetch_and(0);
		std::cerr << "cache: " << rco.cache <<"\n";
		if(rco.cache == 1)
			consensus_one_partition_can_cache(iter->file_name.c_str(), iter->min_seq_id, iter->max_seq_id, rco, reads, out, true);
		else
			consensus_one_partition_can(iter->file_name.c_str(), iter->min_seq_id, iter->max_seq_id, rco, reads, out);
		FILE* job_finished_out = fopen(job_finished, "w");
		fclose(job_finished_out);
	}
	
	return 0;
}
