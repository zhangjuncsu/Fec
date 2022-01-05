#ifndef DW_H
#define DW_H

#include <algorithm>
#include <unordered_map>
#include <array>
#include <atomic>
#include <string>
#include <iostream>

#include "../common/alignment.h"
#include "../common/defs.h"
#include "../common/packed_db.h"

namespace ns_banded_sw {

struct SW_Parameters
{
    idx_t segment_size;
    idx_t row_size;
    idx_t column_size;
    idx_t segment_aln_size;
    idx_t max_seq_size;
    idx_t max_aln_size;
    idx_t d_path_size;
    idx_t aln_path_size;
};

SW_Parameters
get_sw_parameters_small();

SW_Parameters
get_sw_parameters_large();

struct Alignment
{
    int aln_str_size;
    int dist;
    int aln_q_s;
    int aln_q_e;
    int aln_t_s;
    int aln_t_e;
    char* q_aln_str;
    char* t_aln_str;

    void init()
    {
        aln_str_size = 0;
        aln_q_s = aln_q_e = 0;
        aln_t_s = aln_t_e = 0;
    }

    Alignment(const idx_t max_aln_size)
    {
        safe_malloc(q_aln_str, char, max_aln_size);
        safe_malloc(t_aln_str, char, max_aln_size);
    }
    ~Alignment()
    {
        safe_free(q_aln_str);
        safe_free(t_aln_str);
    }
};

struct OutputStore
{
    char* left_store1;
    char* left_store2;
    char* right_store1;
    char* right_store2;
    char* out_store1;
    char* out_store2;
    char* out_match_pattern;

    int left_store_size;
    int right_store_size;
    int out_store_size;
    int query_start, query_end;
    int target_start, target_end;
    int mat, mis, ins, del;
	double ident;
    
    OutputStore(const idx_t max_aln_size)
    {
        safe_malloc(left_store1, char, max_aln_size);
        safe_malloc(left_store2, char, max_aln_size);
        safe_malloc(right_store1, char, max_aln_size);
        safe_malloc(right_store2, char, max_aln_size);
        safe_malloc(out_store1, char, max_aln_size);
        safe_malloc(out_store2, char, max_aln_size);
        safe_malloc(out_match_pattern, char, max_aln_size);
    }

    ~OutputStore()
    {
        safe_free(left_store1);
        safe_free(left_store2);
        safe_free(right_store1);
        safe_free(right_store2);
        safe_free(out_store1);
        safe_free(out_store2);
        safe_free(out_match_pattern);
    }

    void init()
    {
        left_store_size = right_store_size = out_store_size = 0;
    }
};

struct DPathData
{
    int pre_k, x1, y1, x2, y2;
};

struct DPathData2
{
    int d, k, pre_k, x1, y1, x2, y2;
};

struct PathPoint
{
    int x, y;
};

struct DiffRunningData
{
    SW_Parameters   swp;
    char*           query;
    char*           target;
    int*            DynQ;
    int*            DynT;
    Alignment*      align;
    OutputStore*    result;
    DPathData2*     d_path;
    PathPoint*      aln_path;
	
	DiffRunningData(const SW_Parameters& swp_in);
	~DiffRunningData();
};

void fill_m4record_from_output_store(const OutputStore& result, 
									 const idx_t qid, 
									 const idx_t sid,
									 const char qstrand,
									 const char sstrand,
									 const idx_t qsize,
									 const idx_t ssize,
									 const idx_t q_off_in_aln, 
									 const idx_t s_off_in_aln, 
									 const idx_t q_ext,
									 const idx_t s_ext,
									 M4Record& m4);

struct CandidateStartPosition
{
	idx_t qoff;
	idx_t toff;
	idx_t tstart;
	idx_t tsize;
	idx_t tid;
	int left_q, left_t;
	int right_q, right_t;
	int num1, num2;
	int score;
	idx_t toff_in_aln;
	char chain;
};

void print_candidate(const CandidateStartPosition& csp);

int Align(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, Alignment* align, 
		  int* V, int* U, DPathData2* d_path, PathPoint* aln_path, const int right_extend);

void dw_in_one_direction(const char* query, const int query_size, const char* target, const int target_size,
						 int* U, int* V, Alignment* align, DPathData2* d_path, PathPoint* aln_path, 
						 SW_Parameters* swp, OutputStore* result, const int right_extend);

int  dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        int* U, int* V, Alignment* align, DPathData2* d_path, 
        PathPoint* aln_path, OutputStore* result, SW_Parameters* swp,
	    double error_rate, const int min_aln_size);

bool GetAlignment(const char* query, const int query_start, const int query_size,
				  const char* target, const int target_start, const int target_size,
				  DiffRunningData* drd, M5Record& m5, double error_rate, 
				  const int min_aln_size);

template<typename T, size_t N>
struct ArrayHash {
    size_t operator()(const std::array<T, N> & r) const {
        size_t h = 17;
        for (size_t i = 0; i < N; ++i) {
            h = h * 31 + r[i];
        }
        return h;
    }
};

template<typename T, size_t N>
struct ArrayCompare {
    size_t operator()(const std::array<T, N> &a, const std::array<T, N> & b) const {
        for (size_t i = 0; i < N; ++i) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }
};

class AlignmentResult{
public:
	void Swap(bool sameDirection);
	static void Rearrange(std::string &alq, std::string &alt);
	bool TrimEnds(size_t checklen = 1000, int stub = 4);

	size_t query_start{0};
	size_t query_end{0};
	size_t query_size{0};
	size_t target_start{0};
	size_t target_end{0};
	size_t target_size{0};
	std::string aligned_query;
	std::string aligned_target;
    std::string match_pattern;
};

class AlignmentCache{
public: 
	AlignmentResult GetAlignment(idx_t qid, idx_t tid, bool sameDirection);
	void SetAlignment(idx_t qid, idx_t tid, const AlignmentResult &result);
	bool HasAlignment(idx_t qid, idx_t tid) const;
        void Clear() {cache_.clear(); ids_.clear();}
	void SetIds(std::vector<idx_t>ids) {ids_ = ids;}
	bool HasIds(idx_t id) {return find(ids_.begin(), ids_.end(), id) != ids_.end();} 
        size_t size() {return cache_.size();}

protected:
	idx_t ToKey(idx_t qid, idx_t sid) const {return (idx_t)(qid << 32 | sid);}
	std::unordered_map<idx_t, AlignmentResult> cache_;
	std::vector<idx_t>ids_;
};

bool GetAlignment(const char* query, const int query_start, const int query_size,
				  const char* target, const int target_start, const int target_size,
				  DiffRunningData* drd, M5Record& m5, double error_rate, 
				  const int min_aln_size, 
                  idx_t tid, int qid, AlignmentCache& cache, bool sameDirection, int& ac, int& cc, bool useCache);


struct TimeCounter {
    TimeCounter(const std::string& n) : name(n) {}
    ~TimeCounter()  { 
        std::cerr << "Record Time(" << name << "):" << count << ", " << ticks / 10  << ", " << value << "\n";
    }

    void Inc(double s) {
        ticks.fetch_add(int(s*10000));
        count.fetch_add(1);
    }

    void Add(long long v) { value.fetch_add(v); }
    struct Mark {
        Mark(TimeCounter& rt) : recTime(rt) { 
            clock_gettime(CLOCK_MONOTONIC, &start);
        }
        ~Mark() { 
            
            clock_gettime(CLOCK_MONOTONIC, &finish);
            
            elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            recTime.Inc(elapsed); 
        }
        TimeCounter& recTime;
        
        struct timespec start, finish;
        double elapsed;

    };

    std::atomic<long long> ticks {0};
    std::atomic<long long> count {0};
    std::atomic<long long> value {0};
    std::string name;
};

} // end of namespace ns_banded_sw

#endif  // DW_H
