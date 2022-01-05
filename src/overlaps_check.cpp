#include <set>
#include <array>
#include <sstream>
#include <fstream>
#include <algorithm>

#include <unistd.h>
#include <sys/resource.h>

#include "overlaps_store.h"
#include "overlaps_check.h"
#include "../common/alignment.h"

const int MIN_COV = 20;
const int MAX_CANDIDATE = 200;
const std::string TARGET_READS = "target_reads.fasta";

size_t get_peak_RSS(){
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    return (size_t)(rusage.ru_maxrss * 1024);
}

int
fix_file_counts(int num_files) {
	if (num_files < 0) {
		num_files = sysconf(_SC_OPEN_MAX) - 10;
	}
	return num_files;
}

void 
generate_partition_index_file_name(const char* m4_file_name, std::string& ret)
{
    ret = m4_file_name;
    ret += ".partition_files";
}

void
generate_partition_file_name(const char* m4_file_name, const index_t part, std::string& ret)
{
    ret = m4_file_name;
    ret += ".part";
    std::ostringstream os;
    os << part;
    ret += os.str();
}

int get_id_by_name(const std::unordered_map<std::string, std::array<int, 2>> &name2id, const std::string name, int &size){
    auto iter = name2id.find(name);
    if (iter != name2id.end()){
        size = iter->second[1];
        return iter->second[0];
    }else{
        size = 0;
        return -1;
    }
}

std::istream &operator>>(std::istream &in, PAF &paf)
{
    std::string line;
    if (!getline(in, line)) return in;
    std::istringstream ins(line);
    ins >> paf.qname 
        >> paf.qlen 
        >> paf.qs 
        >> paf.qe 
        >> paf.strand 
        >> paf.tname 
        >> paf.tlen 
        >> paf.ts 
        >> paf.te 
        >> paf.rm 
        >> paf.bl 
        >> paf.mq;
    return in;
}

std::ostream &operator<<(std::ostream &out, const PAF &paf)
{
    const char delim = '\t';
    out << paf.qname << delim
        << paf.qid << delim
        << paf.qlen << delim
        << paf.qs << delim
        << paf.qe << delim
        << paf.strand << delim
        << paf.tname << delim
        << paf.tid << delim
        << paf.tlen << delim
        << paf.ts << delim
        << paf.te << delim
        << paf.rm << delim
        << paf.bl << delim
        << paf.mq << std::endl;
    return out;
}

inline void paf_to_candidate(const PAF &paf, ExtensionCandidate &ec)
{
    ec.qid = paf.qid;
    ec.qsize = paf.qlen;
    ec.qoff = paf.qs;
    ec.qend = paf.qe;
    ec.sid = paf.tid;
    ec.ssize = paf.tlen;
    ec.sdir = 0;
    ec.soff = paf.ts;
    ec.send = paf.te;
    ec.score = paf.bl;
    ec.sext = paf.ts;
    if (paf.strand == '+'){
        ec.qdir = 0;
        ec.qext = paf.qs;
    }else{
        ec.qdir = 1;
        ec.qext = paf.qe - 1;
    }
}

inline bool check_pafrecord_mapping_range(const PAF &paf, const double min_cov_ratio)
{
    const int qm = paf.qe - paf.qs;
    const int qs = paf.qlen * min_cov_ratio;
    const int rm = paf.te - paf.ts;
    const int rs = paf.tlen * min_cov_ratio;
    return qm >= qs || rm >= rs;
}

void normalize_pafrecord(PAF &src, PAF &dst, bool is_target)
{
    if (is_target){
        dst = src;
    }else{
        dst.qid = src.tid;
        dst.qname = src.tname;
        dst.qs = src.ts;
        dst.qe = src.te;
        dst.qlen = src.tlen;
        dst.tid = src.qid;
        dst.tname = src.qname;
        dst.ts = src.qs;
        dst.te = src.qe;
        dst.tlen = src.qlen;
        dst.strand = src.strand;
        dst.rm = src.rm;
        dst.mq = src.mq;
        dst.bl = src.bl;
    }
}

void add_overlap(const char *paf_file, const double min_cov_ratio, const int min_read_size, std::unordered_map<std::string, std::array<int, 2>> &name2id, std::unordered_map<int, std::vector<uint64_t>> &paf_ovls)
{
    std::ifstream in;
    PAF paf;
    in.open(paf_file, std::ios::in);
    if(in){
        while (in >> paf){
            if(paf.qname == paf.tname) continue;
            if(paf.qlen < min_read_size || paf.tlen < min_read_size) continue;
            int qsize, tsize;
            int qid = get_id_by_name(name2id, paf.qname, qsize);
            int tid = get_id_by_name(name2id, paf.tname, tsize);
            if (qid == -1 || tid == -1) continue;
            paf.qid = qid;
            paf.tid = tid;
            paf_ovls[tid].push_back((uint64_t)paf.qid << 32 | paf.bl);
        }
        in.close();
    }
}

idx_t get_num_reads_from_paf(const char *paf_file, const int min_read_size, const double min_cov_ratio, idx_t &count, std::unordered_map<std::string, std::array<int, 2>> &name2id, std::unordered_map<int, std::vector<uint64_t>> &paf_ovls)
{
    count = 0;
    std::ifstream in;
    open_fstream(in, paf_file, std::ios::in);
    PAF paf, npaf;
    int max_id = -1;
    while (in >> paf){
        if(paf.qname == paf.tname) continue;
        if (paf.qlen < min_read_size || paf.tlen < min_read_size) continue;
        int qsize, tsize;
        int qid = get_id_by_name(name2id, paf.qname, qsize);
        int tid = get_id_by_name(name2id, paf.tname, tsize);
        if (qid == -1 || tid == -1) continue;
        paf.qid = qid;
        paf.tid = tid;
        count++;
        max_id = std::max(max_id, paf.qid);
        max_id = std::max(max_id, paf.tid);
        normalize_pafrecord(paf, npaf, false);
        paf_ovls[npaf.tid].push_back((uint64_t)npaf.qid << 32 | npaf.bl);
        normalize_pafrecord(paf, npaf, true);
        paf_ovls[npaf.tid].push_back((uint64_t)npaf.qid << 32 | npaf.bl);
    }
    close_fstream(in);
    std::cout << __func__ << " size of paf overlap: " << paf_ovls.size() << std::endl;
    return max_id + 1;
}

void check_overlap(PackedDB &reads,
                   const int min_cov,
                   std::unordered_map<int, std::vector<uint64_t>> &ovls,
                   std::vector<std::string> &names,
                   std::unordered_map<std::string, std::array<int, 2>> name2id)
{
    std::vector<char> read;
    read.reserve(MAX_SEQ_SIZE);
    std::ofstream out;
    open_fstream(out, TARGET_READS, std::ios::out);
    for (size_t i = 0; i < names.size(); ++i){
        if (ovls.find(i) == ovls.end() || ovls[i].size() < MIN_COV){
            int size = name2id[names[i]][1];
            out << ">" << names[i] << "\n";
            reads.get_decode_sequence(i, 0, size, true, read.data());
            std::string r(read.data(), size);
            out << r << "\n";
        }
    }
    close_fstream(out);
}

struct CmpU64{
    bool operator()(const uint64_t &a, const uint64_t &b){
        return ((uint32_t)a) > ((uint32_t)b);
    }
};

void partition_paf(PackedDB &reads,const char *reads_file, const char *paf_file, const double min_cov_ratio, const index_t batch_size, const int min_read_size, int num_files, const int min_cov, const int threads, std::vector<std::string> &names, std::unordered_map<std::string, std::array<int, 2>> &name2id, bool use_cache, bool so, const char* mmc, const char* of)
{
    DynamicTimer dt(__func__);

    num_files = fix_file_counts(num_files);
    idx_t number, count = 0;
    std::unordered_map<int, std::vector<uint64_t>>ovls;
    const idx_t num_reads = get_num_reads_from_paf(paf_file, min_read_size, min_cov_ratio, number, name2id, ovls);
    std::cerr << "first overlap memory usage: " << get_peak_RSS() / 1024 / 1024 << " Mb" << std::endl;
    if(so)
        check_overlap(reads, min_cov, ovls, names, name2id);
    for (auto iter = ovls.begin(); iter != ovls.end(); iter++){
        if (iter->second.size() < MIN_COV){
            iter->second.clear();
        }
        if(iter->second.size() > MAX_CANDIDATE){
            sort(iter->second.begin(), iter->second.end(), CmpU64());
            iter->second.resize(MAX_CANDIDATE);
        }
    }
    std::string paf_file2 = "minimap2output.paf";
    if(so){
        std::ostringstream cmd;
        cmd << mmc << " -K 2g -m 100 -p 0 -g 10000 --max-chain-skip 25 ";
        cmd << TARGET_READS << " " << reads_file;
        cmd << " | awk '{ if($9 - $8 >= " << of << " * $7) print $0}' > " << paf_file2;
        system(cmd.str().c_str());
        add_overlap(paf_file2.c_str(), min_cov_ratio, min_read_size, name2id, ovls);
    }
    idx_t num_batches;
    std::unordered_map<int, bool> done;
    std::vector<int> group_ids;
    std::vector<int> ids;
    std::vector<std::vector<int>> index;
    if(use_cache){
        for (auto iter = ovls.begin(); iter != ovls.end(); iter++){
            if(iter->second.size() < min_cov){
                iter->second.resize(0);
                done[iter->first] = true;
                continue;
            }
            done[iter->first] = false;
            ids.push_back(iter->first);
            if(iter->second.size() > MAX_CANDIDATE){
                sort(iter->second.begin(), iter->second.end(), CmpU64());
                iter->second.resize(MAX_CANDIDATE);
            }
        }
        sort(ids.begin(), ids.end());
        for (auto id : ids){
            if (!done[id]){
                group_ids.push_back(id);
                done[id] = true;
                auto ecids = ovls[id];
                for (auto ecid : ecids){
                    auto iter_done = done.find(ecid >> 32);
                    if (iter_done != done.end() && !iter_done->second){
                        group_ids.push_back(ecid >> 32);
                        iter_done->second = true;
                    }
                }
            }
        }
        index.resize(names.size());
        for(size_t i = 0, size = group_ids.size(); i < size; ++i){
            if(!ovls[group_ids[i]].empty()) index[group_ids[i]].push_back(i);
            for(auto qid: ovls[group_ids[i]]){
                index[group_ids[i]].push_back(qid >> 32);
            }
        }
        std::cerr << "second overlap memory usage: " << get_peak_RSS() / 1024 / 1024 << " Mb" << std::endl;
        num_batches = (group_ids.size() + batch_size - 1) / batch_size;
    }else{
        num_batches = (num_reads + batch_size - 1) / batch_size;
    }
    std::string idx_file_name;
    generate_partition_index_file_name(paf_file, idx_file_name);
    std::ofstream idx_file;
    open_fstream(idx_file, idx_file_name.c_str(), std::ios::out);

    ExtensionCandidate ec;
    PartitionResultsWriter<ExtensionCandidate> prw(num_files);
    for (index_t i = 0; i < num_batches; i += num_files)
    {
        const index_t sfid = i;
        const index_t efid = std::min(sfid + num_files, num_batches);
        const int nf = efid - sfid;
        const index_t L = batch_size * sfid;
        const index_t R = batch_size * efid;
        std::cout << "Lid: " << L << ", Rid: " << R << std::endl;
        prw.OpenFiles(sfid, efid, paf_file, generate_partition_file_name);
        PAF paf, npaf;
        std::ifstream in;
        open_fstream(in, paf_file, std::ios::in);
        while(in >> paf){
            int qsize, tsize;
            int qid = get_id_by_name(name2id, paf.qname, qsize);
            int tid = get_id_by_name(name2id, paf.tname, tsize);
            if(qid == -1 || tid == -1) continue;
            paf.qid = qid; paf.tid = tid;
            if(use_cache){
                if(!index[paf.tid].empty()){
                    if(index[paf.tid][0] >= L && index[paf.tid][0] < R && std::find(index[paf.tid].begin() + 1, index[paf.tid].end(), paf.qid) != index[paf.tid].end()){
                        paf_to_candidate(paf, ec);
                        prw.WriteOneResult((index[paf.tid][0] - L) / batch_size, paf.tid, ec);
                    }
                }
                if(!index[paf.qid].empty()){
                    if(index[paf.qid][0] >= L && index[paf.qid][0] < R && std::find(index[paf.qid].begin() + 1, index[paf.qid].end(), paf.tid) != index[paf.qid].end()){
                        normalize_pafrecord(paf, npaf, false);
                        paf_to_candidate(npaf, ec);
                        prw.WriteOneResult((index[paf.qid][0] - L) / batch_size, paf.qid, ec);
                    }
                }
            }else{
                if(tid >= L && tid < R){
                    paf_to_candidate(paf, ec);
                    prw.WriteOneResult((ec.sid - L) / batch_size, ec.sid, ec);
                }
                if(qid >= L && qid < R){
                    normalize_pafrecord(paf, npaf, false);
                    paf_to_candidate(npaf, ec);
                    prw.WriteOneResult((ec.sid - L) / batch_size, ec.sid, ec);
                }
            }
        }
        close_fstream(in);
        if(so && access(paf_file2.c_str(), F_OK) != -1){
            in.open(paf_file2, std::ios::in);
            if(in){
                while(in >> paf){
                    int qsize, tsize;
                    int qid = get_id_by_name(name2id, paf.qname, qsize);
                    int tid = get_id_by_name(name2id, paf.tname, tsize);
                    if(qid == -1 || tid == -1) continue;
                    paf.qid = qid; paf.tid = tid;
                    if(use_cache){
                        if(!index[paf.tid].empty()){
                            if(index[paf.tid][0] >= L && index[paf.tid][0] < R && std::find(index[paf.tid].begin() + 1, index[paf.tid].end(), paf.qid) != index[paf.tid].end()){
                                paf_to_candidate(paf, ec);
                                prw.WriteOneResult((index[paf.tid][0] - L) / batch_size, paf.tid, ec);
                            }
                        }
                        if(!index[paf.qid].empty()){
                            if(index[paf.qid][0] >= L && index[paf.qid][0] < R && std::find(index[paf.qid].begin() + 1, index[paf.qid].end(), paf.tid) != index[paf.qid].end()){
                                normalize_pafrecord(paf, npaf, false);
                                paf_to_candidate(npaf, ec);
                                prw.WriteOneResult((index[paf.qid][0] - L) / batch_size, paf.qid, ec);
                            }
                        }
                    }else{
                        if(tid >= L && tid < R){
                            paf_to_candidate(paf, ec);
                            prw.WriteOneResult((ec.sid - L) / batch_size, ec.sid, ec);
                        }
                        if(qid >= L && qid < R){
                            normalize_pafrecord(paf, npaf, false);
                            paf_to_candidate(npaf, ec);
                            prw.WriteOneResult((ec.sid - L) / batch_size, ec.sid, ec);
                        }
                    }
                }
                in.close();
            }
        }

        for (int k = 0; k < nf; ++k)
        {
            if (prw.max_seq_ids[k] == std::numeric_limits<index_t>::min())
                continue;
            idx_file << prw.file_names[k] << "\t" << prw.min_seq_ids[k] << "\t" << prw.max_seq_ids[k] << "\n";
            fprintf(stderr, "%s contains reads %d --- %d\n", prw.file_names[k].c_str(), (int)prw.min_seq_ids[k], (int)prw.max_seq_ids[k]);
        }

        prw.CloseFiles();
    }
    close_fstream(idx_file);
    std::cout << "group ids: " << group_ids.size() << "\tdone: " << done.size() << "\tovls: " << ovls.size() << "\tcount: " << count << std::endl;
}

void
load_partition_files_info(const char* idx_file_name, std::vector<PartitionFileInfo>& file_info_vec)
{
    file_info_vec.clear();
    std::ifstream in;
    open_fstream(in, idx_file_name, std::ios::in);
    PartitionFileInfo pfi;
    while (in >> pfi.file_name)
    {
        in >> pfi.min_seq_id >> pfi.max_seq_id;
        file_info_vec.push_back(pfi);
    }
    close_fstream(in);
}