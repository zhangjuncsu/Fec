#include "dw.h"

namespace ns_banded_sw {

#define GAP_ALN 4

SW_Parameters
get_sw_parameters_small()
{
    SW_Parameters swp;
    swp.segment_size = 500;
    swp.row_size = 4096;
    swp.column_size = 4096;
    swp.segment_aln_size = 4096;
    swp.max_seq_size = MAX_SEQ_SIZE;
    swp.max_aln_size = MAX_SEQ_SIZE;
    swp.d_path_size = 5000000;
    swp.aln_path_size = 5000000;

    return swp;
}

SW_Parameters
get_sw_parameters_large()
{
    SW_Parameters swp;
    swp.segment_size = 1000;
    swp.row_size = 4096;
    swp.column_size = 4096;
    swp.segment_aln_size = 4096;
    swp.max_seq_size = MAX_SEQ_SIZE;
    swp.max_aln_size = MAX_SEQ_SIZE;
    swp.d_path_size = 5000000;
    swp.aln_path_size = 5000000;

    return swp;
}

DiffRunningData::DiffRunningData(const SW_Parameters& swp_in)
{
	swp = swp_in;
	safe_malloc(query, char, swp.max_seq_size);
	safe_malloc(target, char, swp.max_seq_size);
	safe_malloc(DynQ, int, swp.row_size);
	safe_malloc(DynT, int, swp.column_size);
	align = new Alignment(swp.segment_aln_size);
	result = new OutputStore(swp.max_aln_size);
	safe_malloc(d_path, DPathData2, swp.d_path_size);
	safe_malloc(aln_path, PathPoint, swp.aln_path_size);
}

DiffRunningData::~DiffRunningData()
{
	safe_free(query);
	safe_free(target);
	safe_free(DynQ);
	safe_free(DynT);
	delete align;
	delete result;
	safe_free(d_path);
	safe_free(aln_path);
}

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
									 M4Record& m4)
{
	m4qid(m4) = qid;
	m4sid(m4) = sid;
	m4ident(m4) = result.ident;
	m4vscore(m4) = 100;
	m4qdir(m4) = 1 - (qstrand == 'F');
	m4qoff(m4) = q_off_in_aln + result.query_start;
	m4qend(m4) = q_off_in_aln + result.query_end;
	m4qsize(m4) = qsize;
	m4sdir(m4) = 1 - (sstrand == 'F');
	m4soff(m4) = s_off_in_aln + result.target_start;
	m4send(m4) = s_off_in_aln + result.target_end;
	m4ssize(m4) = ssize;
	m4qext(m4) = q_ext;
	m4sext(m4) = s_ext;
	
	if (m4qdir(m4) == 1)
	{
		m4qoff(m4) = qsize - (q_off_in_aln + result.query_end);
		m4qend(m4) = qsize - (q_off_in_aln + result.query_start);
		m4qext(m4) = qsize - 1 - q_ext;
	}
	if (m4sdir(m4) == 1)
	{
		m4soff(m4) = ssize - (s_off_in_aln + result.target_end);
		m4send(m4) = ssize - (s_off_in_aln + result.target_start);
		m4sext(m4) = ssize - 1 - s_ext;
	}
}

void print_candidate(const CandidateStartPosition& csp)
{
	std::cout << "qoff = " << csp.qoff
			  << ", toff = " << csp.toff
			  << ", tstart = " << csp.tstart
			  << ", tsize = " << csp.tsize
			  << ", tid = " << csp.tid
			  << ", num1 = " << csp.num1
			  << ", num2 = " << csp.num2
			  << ", score = " << csp.score
			  << ",  chain = " << csp.chain
			  << ", l1 = " << csp.left_q
			  << ", r1 = " << csp.right_q
			  << ", l2 = " << csp.left_t
			  << ", r2 = " << csp.right_t
			  << "\n";
}

int CompareDPathData2(const void* a, const void* b)
{
    const DPathData2* d1 = (const DPathData2*)a;
    const DPathData2* d2 = (const DPathData2*)b;
    return (d1->d == d2->d) ? (d1->k - d2->k) : (d1->d - d2->d);
}

struct SCompareDPathData2
{
    bool operator()(const DPathData2& a, const DPathData2& b)
    { return (a.d == b.d) ? (a.k < b.k) : (a.d < b.d); }
};

DPathData2* GetDPathIdx(const int d, const int k, const unsigned int max_idx, DPathData2* base)
{
    DPathData2 target;
    target.d = d;
    target.k = k;
    DPathData2* ret = (DPathData2*)bsearch(&target, base, max_idx, sizeof(DPathData2), CompareDPathData2);
    return ret;
}

int Align(const char* query, const int q_len, const char* target, const int t_len, 
          const int band_tolerance, const int get_aln_str, Alignment* align, 
		  int* V, int* U, DPathData2* d_path, PathPoint* aln_path, 
		  const int right_extend, double error_rate)
{
    int k_offset;
    int  d;
    int  k, k2;
    int best_m;
    int min_k, new_min_k, max_k, new_max_k, pre_k;
    int x, y;
    int ck, cd, cx, cy, nx, ny; 
    int max_d, band_size;
    unsigned long d_path_idx = 0, max_idx = 0;
    int aln_path_idx, aln_pos, i, aligned = 0;
    DPathData2* d_path_aux;
    
    max_d = (int)(2.0 * error_rate * (q_len + t_len));
    k_offset = max_d;
    band_size = band_tolerance * 2;
    align->init();
    best_m = -1;
    min_k = 0;
    max_k = 0;
    d_path_idx = 0;
    max_idx = 0;
    
    for (d = 0; d < max_d; ++d)
    {
        if (max_k - min_k > band_size) break;

        for (k = min_k; k <= max_k; k += 2)
        {
            if( k == min_k || (k != max_k && V[k - 1 + k_offset] < V[k + 1 + k_offset]) )
            { pre_k = k + 1; x = V[k + 1 + k_offset]; }
            else 
            { pre_k = k - 1; x = V[k - 1 + k_offset] + 1; }
            y = x - k;
            d_path[d_path_idx].d = d;
            d_path[d_path_idx].k = k;
            d_path[d_path_idx].x1 = x;
            d_path[d_path_idx].y1 = y;
			
			if (right_extend)
				while( x < q_len && y < t_len && query[x] == target[y]) { ++x; ++y; }
			else
				while( x < q_len && y < t_len && query[-x] == target[-y]) { ++x; ++y; }

            d_path[d_path_idx].x2 = x;
            d_path[d_path_idx].y2 = y;
            d_path[d_path_idx].pre_k = pre_k;
            ++d_path_idx;

            V[k + k_offset] = x;
            U[k + k_offset] = x + y;
            best_m = std::max(best_m, x + y);
            if (x >= q_len || y >= t_len)
            { aligned = 1; max_idx = d_path_idx; break; }
        }

        // for banding
        new_min_k = max_k;
        new_max_k = min_k;
        for (k2 = min_k; k2 <= max_k; k2 += 2)
            if (U[k2 + k_offset] >= best_m - band_tolerance)
            { new_min_k = std::min(new_min_k, k2); new_max_k = std::max(new_max_k, k2); }
        max_k = new_max_k + 1;
        min_k = new_min_k - 1;

        if (aligned)
        {
            align->aln_q_e = x;
            align->aln_t_e = y;
            align->dist = d;
            align->aln_str_size = (x + y + d) / 2;
            align->aln_q_s = 0;
            align->aln_t_s = 0;

            if (get_aln_str)
            {
                cd = d;
                ck = k;
                aln_path_idx = 0;
                while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
                {
                    d_path_aux = GetDPathIdx(cd, ck, max_idx, d_path);
                    aln_path[aln_path_idx].x = d_path_aux->x2;
                    aln_path[aln_path_idx].y = d_path_aux->y2;
                    ++aln_path_idx;
                    aln_path[aln_path_idx].x = d_path_aux->x1;
                    aln_path[aln_path_idx].y = d_path_aux->y1;
                    ++aln_path_idx;
                    ck = d_path_aux->pre_k;
                    cd -= 1;
                }
                --aln_path_idx;
                cx = aln_path[aln_path_idx].x;
                cy = aln_path[aln_path_idx].y;
                align->aln_q_s = cx;
                align->aln_t_s = cy;
                aln_pos = 0;
                while (aln_path_idx > 0)
                {
                    --aln_path_idx;
                    nx = aln_path[aln_path_idx].x;
                    ny = aln_path[aln_path_idx].y;
                    if (cx == nx && cy == ny) continue;
                    if (cx == nx && cy != ny)
                    {
						if (right_extend)
						{
							for (i = 0; i < ny - cy; ++i) align->q_aln_str[aln_pos + i] = GAP_ALN;
							for (i = 0; i < ny - cy; ++i) align->t_aln_str[aln_pos + i] = target[cy + i];
						}
						else
						{
							for (i = 0; i < ny - cy; ++i) align->q_aln_str[aln_pos + i] = GAP_ALN;
							for (i = 0; i < ny - cy; ++i) align->t_aln_str[aln_pos + i] = target[-(cy + i)];
						}
                        aln_pos += ny - cy;
                    }
                    else if (cx != nx && cy == ny)
                    {
						if (right_extend)
						{
							for (i = 0; i < nx - cx; ++i) align->q_aln_str[aln_pos + i] = query[cx + i];
							for (i = 0; i < nx - cx; ++i) align->t_aln_str[aln_pos + i] = GAP_ALN;
						}
						else
						{
							for (i = 0; i < nx - cx; ++i) align->q_aln_str[aln_pos + i] = query[-(cx + i)];
							for (i = 0; i < nx - cx; ++i) align->t_aln_str[aln_pos + i] = GAP_ALN;
						}
                        aln_pos += nx - cx;
                    }
                    else
                    {
						if (right_extend)
						{
							for (i = 0; i < nx - cx; ++i) align->q_aln_str[aln_pos + i] = query[cx + i];
							for (i = 0; i < ny - cy; ++i) align->t_aln_str[aln_pos + i] = target[cy + i];
						}
						else
						{
							for (i = 0; i < nx - cx; ++i) align->q_aln_str[aln_pos + i] = query[-(cx + i)];
							for (i = 0; i < ny - cy; ++i) align->t_aln_str[aln_pos + i] = target[-(cy + i)];
						}
                        aln_pos += ny - cy;
                    }
                    cx = nx;
                    cy = ny;
                }
                align->aln_str_size = aln_pos;
            }
            break;
        }
    }
    if (align->aln_q_e == q_len || align->aln_t_e == t_len) return 1;
    else return 0;
}

void dw_in_one_direction(const char* query, const int query_size, const char* target, const int target_size,
						 int* U, int* V, Alignment* align, DPathData2* d_path, PathPoint* aln_path, 
						 SW_Parameters* swp, OutputStore* result, const int right_extend, double error_rate)
{
	const idx_t ALN_SIZE = swp->segment_size;
	const idx_t U_SIZE = swp->row_size;
	const idx_t V_SIZE = swp->column_size;
    int extend1 = 0, extend2 = 0;
    const char* seq1 = query;
    const char* seq2 = target;
    int extend_size = std::min(query_size, target_size);
    int seg_size;
    int flag_end = 1;
    int align_flag;
    int i, j, k, num_matches;
    while (flag_end)
    {
        if (extend_size > (ALN_SIZE + 100))
        { seg_size = ALN_SIZE; }
        else
        { seg_size = extend_size; flag_end = 0; }
        memset(U, 0, sizeof(int) * U_SIZE);
        memset(V, 0, sizeof(int) * V_SIZE);
        if (right_extend) { seq1 = query + extend1; seq2 = target + extend2; }
        else { seq1 = query - extend1; seq2 = target - extend2; }
        align_flag = Align(seq1, seg_size, seq2, seg_size, 0.3 * seg_size, 400, align, U, V, d_path, aln_path, right_extend, error_rate);
        if (align_flag)
        {
            for (k = align->aln_str_size - 1, i = 0, j = 0, num_matches = 0; k > -1 && num_matches < 4; --k)
            {
                if (align->q_aln_str[k] != GAP_ALN) ++i;
                if (align->t_aln_str[k] != GAP_ALN) ++j;
                if (align->q_aln_str[k] == align->t_aln_str[k]) ++num_matches;
                else num_matches = 0;
            }
            if (flag_end)
            {
                i = ALN_SIZE - align->aln_q_e + i;
                j = ALN_SIZE - align->aln_t_e + j;
                if (i == ALN_SIZE) align_flag = 0;
				extend1 = extend1 + ALN_SIZE - i; extend2 = extend2 + ALN_SIZE - j;
            }
            else
            {
                i = extend_size - align->aln_q_e;
                j = extend_size - align->aln_t_e;
                if (i == extend_size) align_flag = 0;
				extend1 += (extend_size - i); extend2 += (extend_size - j);
                k = align->aln_str_size - 1;
            }
            if (align_flag)
            {
                if (right_extend)
                {
                    memcpy(result->right_store1 + result->right_store_size, align->q_aln_str, k + 1);
                    memcpy(result->right_store2 + result->right_store_size, align->t_aln_str, k + 1);
                    result->right_store_size += (k + 1);
                }
                else
                {
                    memcpy(result->left_store1 + result->left_store_size, align->q_aln_str, k + 1);
                    memcpy(result->left_store2 + result->left_store_size, align->t_aln_str, k + 1);
                    result->left_store_size += (k + 1);
                }
                extend_size = std::min(query_size - extend1, target_size - extend2);
            }
        }
        if (!align_flag) break;
    }
}

int  dw(const char* query, const int query_size, const int query_start,
        const char* target, const int target_size, const int target_start,
        int* U, int* V, Alignment* align, DPathData2* d_path, 
        PathPoint* aln_path, OutputStore* result, SW_Parameters* swp,
	    double error_rate, const int min_aln_size)
{
    result->init();
    align->init();
    // left extend
    dw_in_one_direction(query + query_start - 1, query_start,
						target + target_start - 1, target_start,
						U, V, align, d_path, aln_path, swp, result, 
						0, error_rate);
    align->init();
    // right extend
    dw_in_one_direction(query + query_start, query_size - query_start,
						target + target_start, target_size - target_start,
						U, V, align, d_path, aln_path, swp, result, 
						1, error_rate);

    // merge the results
    int i, j, k, idx = 0;
    const char* encode2char = "ACGT-";
    for (k = result->left_store_size - 1, i = 0, j = 0; k > - 1; --k, ++idx)
    {
		unsigned char ch = result->left_store1[k];
		r_assert(ch >= 0 && ch <= 4);
		ch = encode2char[ch];
		result->out_store1[idx] = ch;
		if (ch != '-') ++i;
		
		ch = result->left_store2[k];
		r_assert(ch >= 0 && ch <= 4);
		ch = encode2char[ch];
		result->out_store2[idx] = ch;
		if (ch != '-')++j;
    }
    result->query_start = query_start - i;
	if (result->query_start < 0)
	{
		std::cerr << "query_start = " << query_start << ", i = " << i << "\n";
	}
	r_assert(result->query_start >= 0);
    result->target_start = target_start - j;
if (result->target_start < 0)
	{
		std::cerr << "target_start = " << target_start << ", j = " << j << "\n";
	}
	r_assert(result->target_start >= 0);
    for (k = 0, i = 0, j = 0; k < result->right_store_size; ++k, ++idx)
    {

		unsigned char ch = result->right_store1[k];
		r_assert(ch >= 0 && ch <= 4);
		ch = encode2char[ch];
		result->out_store1[idx] = ch;
		if (ch != '-') ++i;
		
		ch = result->right_store2[k];
		r_assert(ch >= 0 && ch <= 4);
		ch = encode2char[ch];
		result->out_store2[idx] = ch;
		if (ch != '-') ++j;
    }
    result->out_store_size = idx;
    result->query_end = query_start + i;
    result->target_end = target_start + j;

	if (result->out_store_size >= min_aln_size)
    {
        int mat = 0, mis = 0, ins = 0, del = 0;
        for (j = 0; j < result->out_store_size; ++j)
        {
            if (result->out_store1[j] == result->out_store2[j])
            {
                ++mat;
                result->out_match_pattern[j] = '|';
            }
            else if (result->out_store1[j] == '-')
            {
                ++ins;
                result->out_match_pattern[j] = '*';
            }
            else if (result->out_store2[j] == '-')
            {
                ++del;
                result->out_match_pattern[j] = '*';
            }
            else
            {
                ++mis;
                result->out_match_pattern[j] = '*';
            }
        }
        result->out_store1[result->out_store_size] = '\0';
        result->out_store2[result->out_store_size] = '\0';
        result->out_match_pattern[result->out_store_size] = '\0';
        result->mat = mat;
        result->mis = mis;
        result->ins = ins;
        result->del = del;
        result->ident = 100.0 * mat / result->out_store_size;

        return 1;
    }
    return 0;
}

bool GetAlignment(const char* query, const int query_start, const int query_size,
				  const char* target, const int target_start, const int target_size,
				  DiffRunningData* drd, M5Record& m5, double error_rate,
				  const int min_aln_size)
{
	int flag = dw(query, query_size, query_start,
				  target, target_size, target_start,
				  drd->DynQ, drd->DynT, 
				  drd->align, drd->d_path,
				  drd->aln_path, drd->result,
				  &drd->swp, error_rate, min_aln_size);
	if (!flag) return false;
	
	int qrb = 0, qre = 0;
	int trb = 0, tre = 0;
	int eit = 0, k = 0;
	const int consecutive_match_region_size = 4;
	for (k = 0; k < drd->result->out_store_size && eit < consecutive_match_region_size; ++k)
	{
		const char qc = drd->result->out_store1[k];
		const char tc = drd->result->out_store2[k];
		if (qc != '-') ++qrb;
		if (tc != '-') ++trb;
		if (qc == tc) ++eit;
		else eit = 0;
	}
	if (eit < consecutive_match_region_size) return false;
	k -= consecutive_match_region_size;
	qrb -= consecutive_match_region_size;
	trb -= consecutive_match_region_size;
	const int start_aln_id = k;
	if (start_aln_id < 0)
	{
		std::cout << qrb << "\t" << trb << "\t" << eit << "\t" << k << "\n";
	}
	
	for (k = drd->result->out_store_size - 1, eit = 0; k >= 0 && eit < consecutive_match_region_size; --k)
	{
		const char qc = drd->result->out_store1[k];
		const char tc = drd->result->out_store2[k];
		if (qc != '-') ++qre;
		if (tc != '-') ++tre;
		if (qc == tc) ++eit;
		else eit = 0;
	}
	if (eit < consecutive_match_region_size) return false;
	k += consecutive_match_region_size;
	qre -= consecutive_match_region_size;
	tre -= consecutive_match_region_size;
	const int end_aln_id = k + 1;
	
	m5qsize(m5) = query_size;
	m5qoff(m5) = drd->result->query_start + qrb;
	m5qend(m5) = drd->result->query_end - qre;
	m5qdir(m5) = FWD;

	m5ssize(m5) = target_size;
	m5soff(m5) = drd->result->target_start + trb;
	m5send(m5) = drd->result->target_end - tre;
	m5sdir(m5) = FWD;

	const int aln_size = end_aln_id - start_aln_id;


	memcpy(m5qaln(m5), drd->result->out_store1 + start_aln_id, aln_size);
	memcpy(m5saln(m5), drd->result->out_store2 + start_aln_id, aln_size);
	memcpy(m5pat(m5), drd->result->out_match_pattern + start_aln_id, aln_size);
	m5qaln(m5)[aln_size] = '\0';
	m5saln(m5)[aln_size] = '\0';
	m5pat(m5)[aln_size] = '\0';
	
	return true;
}

// TimeCounter tc("align");
bool GetAlignment(const char* query, const int query_start, const int query_size,
				  const char* target, const int target_start, const int target_size,
				  DiffRunningData* drd, M5Record& m5, double error_rate, 
				  const int min_aln_size, 
                  idx_t tid, int qid, AlignmentCache& cache, bool sameDirection, int& ac, int& cc, bool useCache)
{
    // TimeCounter::Mark m(tc);
    bool has_id = cache.HasIds(qid) && cache.HasIds(tid);
    if(!useCache || !has_id || !cache.HasAlignment(tid, qid)){
        int flag = GetAlignment(query, query_start, query_size, target, target_start, target_size, drd, m5, error_rate, min_aln_size);
        ac++;
        if(!flag) return false;
        AlignmentResult result;
        result.query_start = m5qoff(m5);
        result.query_end = m5qend(m5);
        result.query_size = m5qsize(m5);
        result.target_start = m5soff(m5);
        result.target_end = m5send(m5);
        result.target_size = m5ssize(m5);
        result.aligned_query = m5qaln(m5);
        result.aligned_target = m5saln(m5);
        result.match_pattern = m5pat(m5);
	    if(cache.HasIds(qid) && cache.HasIds(tid)) cache.SetAlignment(qid, tid, result);
        return true;
    }else {
        cc++;
        AlignmentResult result = cache.GetAlignment(tid, qid, sameDirection);
        m5qoff(m5) = result.query_start;
        m5qend(m5) = result.query_end;
        m5qsize(m5) = result.query_size;
        m5qdir(m5) = FWD;
        m5soff(m5) = result.target_start;
        m5send(m5) = result.target_end;
        m5ssize(m5) = result.target_size;
        m5sdir(m5) = FWD;
        strcpy(m5qaln(m5), result.aligned_query.c_str());
        strcpy(m5saln(m5), result.aligned_target.c_str());
        strcpy(m5pat(m5), result.match_pattern.c_str());
        return true;
    }
}

bool GetAlignmentNanopore(const char* query, const int query_start, const int query_size,
				  const char* target, const int target_start, const int target_size,
				  DiffRunningData* drd, M5Record& m5, double error_rate, 
				  const int min_aln_size, OcAlignData* align_data, OcDalignData* dalign_data, FullEdlibAlignData* falign_data, 
                  idx_t tid, int qid, AlignmentCache& cache, bool sameDirection, int& ac, int& cc, bool useCache, bool recuse_long_indel)
{
    // TimeCounter::Mark m(tc);
    bool has_id = cache.HasIds(qid) && cache.HasIds(tid);
    if(!useCache || !has_id || !cache.HasAlignment(tid, qid)){
        int flag = cns_extension(align_data, dalign_data, falign_data, 
                                query, query_start, query_size, 
                                target, target_start, target_size, 
                                min_aln_size, m5, recuse_long_indel);
        ac++;
        if(!flag) return false;
        AlignmentResult result;
        m5qdir(m5) = FWD;
        m5sdir(m5) = FWD;
        m5qsize(m5) = query_size;
        m5ssize(m5) = target_size;
        result.query_start = m5qoff(m5);
        result.query_end = m5qend(m5);
        result.query_size = m5qsize(m5);
        result.target_start = m5soff(m5);
        result.target_end = m5send(m5);
        result.target_size = m5ssize(m5);
        result.aligned_query = m5qaln(m5);
        result.aligned_target = m5saln(m5);
        result.match_pattern = m5pat(m5);
	    if(cache.HasIds(qid) && cache.HasIds(tid)) cache.SetAlignment(qid, tid, result);
        return true;
    }else {
        cc++;
        AlignmentResult result = cache.GetAlignment(tid, qid, sameDirection);
        m5qoff(m5) = result.query_start;
        m5qend(m5) = result.query_end;
        m5qsize(m5) = result.query_size;
        m5qdir(m5) = FWD;
        m5soff(m5) = result.target_start;
        m5send(m5) = result.target_end;
        m5ssize(m5) = result.target_size;
        m5sdir(m5) = FWD;
        strcpy(m5qaln(m5), result.aligned_query.c_str());
        strcpy(m5saln(m5), result.aligned_target.c_str());
        strcpy(m5pat(m5), result.match_pattern.c_str());
        return true;
    }
}

const unsigned char TABLE[256] = {
     0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  '-',  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

unsigned char Complement(char c){
   return TABLE[c];
}

void AlignmentResult::Rearrange(std::string &alq, std::string &alt){
    //  ...CGX...     ---\    ...CGX... 
    //  ...--C...     .../    ...C--...
    //
    auto move_forward = [](std::string &mdf, const std::string& ref) {
        bool modified = false;
        for (size_t i=1; i<mdf.size(); i++) {
            if (mdf[i-1] == '-' && mdf[i] != '-') {
                int cand = -1;
                for (int j = i-1; j>=0 && mdf[j] == '-'; j--) {
                    if (ref[j] == mdf[i]) {
                        cand = j;
                    }
                }

                if (cand != -1 ) {
                    std::swap(mdf[i], mdf[cand]);
                    modified = true;
                }
            }
        }
        return modified;
    };


    //  qlt: ...-GA...     ---\    ...-GA...   ---\    ...GA... 
    //  alt: ...A--...     .../    ...--A...   ---/    ...-A...  
    //
    auto move_backward = [](const std::string &alq, std::string& alt) {
        bool modified = false;
        for (size_t i=0; i+1<alt.size(); i++) {
            if (alt[i] != '-' && alt[i+1] == '-' && alq[i] == '-') {
                int cand = -1;
                for (size_t j = i+1; j < alt.size() && alt[j] == '-'; j++) {
                    if (alq[j] == alt[i]) {
                        cand = j;
                    }
                }

                if (cand != -1 ) {
                    std::swap(alt[i], alt[cand]);
                    modified = true;
                }
            }
        }
        return modified;
    };

//
    // alq: ...A-...
    // alt: ...-G...
    //
    auto swap_indel_dash = [](std::string& alq, std::string& alt) {
        bool modified = false;
        for(size_t i = 1; i < alq.size() - 1; i++) {
            if(alq[i - 1] != '-' && alq[i] == '-') {
                if(alt[i - 1] == '-' && alt[i] != '-') {
                    std::swap(alq[i - 1], alq[i]);
                    std::swap(alt[i - 1], alt[i]);
                    modified = true;
                }
            }
        }
        return modified;
    };

    //
    // alq: ...AAA-...       ---\    ...-AAA...
    // alt: ...---C...       ---/    ...C---...
    //
    auto swap_dash_forward = [](std::string &alq, std::string& alt) {
        bool modified = false;
        for (size_t i=1; i<alt.size(); i++) {
            if (alt[i-1] == '-' && alt[i] != '-' && alq[i] == '-') {
                for (int j = i; j>=1 && alt[j] == '-'; j--) {
                    
                    std::swap(alt[j], alt[j-1]);
                    std::swap(alq[j], alq[j-1]);
                    modified = true;
                }
            }

        }
        return modified;
    };

    //  alq: ...A-C...       ---\    ...AC...
    //  alt: ...A-C...       ---/    ...AC...
    //
    //  alq: ...-A..A-...       ---\    ...A..A...
    //  alt: ...A-..-A...       ---/    ...A..A...
    auto remove_dash = [] (std::string &alq, std::string &alt) {
        size_t oldsize = alq.size();
        for (size_t i=0; i<alt.size(); ++i) {
            if (alq[i] == alt[i] && alt[i] == '-') {
                alt.erase(alt.begin()+i);
                alq.erase(alq.begin()+i);
                i--;
            }
        }

        for (size_t i=1; i< alt.size(); ++i) {
            if (alt[i-1] == '-' && alq[i] == '-'  && alt[i] == alq[i-1]) {
                alt.erase(alt.begin()+i-1);
                alq.erase(alq.begin()+i);
                i--;
                //std::swap(alt[i-1], alt[i]);
            } else if (alt[i] == '-' && alq[i-1] == '-'  && alt[i-1] == alq[i]) {
                alt.erase(alt.begin()+i);
                alq.erase(alq.begin()+i-1);
                i--;
                //std::swap(alq[i-1], alq[i]);
            }
        }
        return oldsize != alq.size();   // modified
    };

    bool finished = false;

    while (!finished) {
        finished = true;
        finished &= !move_forward(alq, alt);
       finished &= !move_backward(alq, alt);
        break;

    }
}

bool AlignmentResult::TrimEnds(size_t checklen, int stub) {
    int sc = 0;
    size_t as = 0, ts = target_start, qs = query_start;
    size_t ias = 0, its = target_start, iqs = query_start;

    for (ias = 0; ias < checklen; ++ias) {
        if (aligned_query[ias] == aligned_target[ias]) {
            if (aligned_query[ias] != '-') {
                if (sc == 0) {
                    as = ias;
                    qs = iqs;
                    ts = its;
                }
                sc ++;
                if (sc >= stub) break;
            } 
        } else  {
            sc = 0;
        } 

        // increase iqs, its
        if (aligned_query[ias] != '-') iqs++;
        if (aligned_target[ias] != '-') its++;
    }

    int ec = 0;
    size_t ae = 0, te = target_end, qe = query_end;
    size_t iae = aligned_target.size(), ite = target_end, iqe = query_end;

    for (; iae > aligned_target.size() - checklen; iae -- ) {
        if (aligned_query[iae-1] == aligned_target[iae-1]) {
            if (aligned_query[iae-1] != '-') {
                if (ec == 0) {
                    ae = iae;
                    qe = iqe;
                    te = ite;
                }
                ec ++;
                if (ec >= stub) break;
            }
        } else {
            ec = 0;
        }
        
        // decrease iqe, ite
        if (aligned_query[iae-1] != '-') iqe--;
        if (aligned_target[iae-1] != '-') ite--; 
    }


    if (sc >= stub && ec >= stub) {
        target_start = ts;
        query_start = qs;

        target_end = te;
        query_end = qe;

        aligned_target = aligned_target.substr(as, ae-as);
        aligned_query  = aligned_query.substr(as, ae-as);
        return true;
    } else {
        return false;
    }

}

void AlignmentResult::Swap(bool sameDirection){
    std::swap(aligned_query, aligned_target);
    std::swap(query_size, target_size);
    std::swap(query_start, target_start);
    std::swap(query_end, target_end);
    
    if(!sameDirection){
        std::swap(target_start, target_end);
        target_start = target_size - target_start;
        target_end = target_size - target_end;
        std::swap(query_start, query_end);
        query_start = query_size - query_start;
        query_end = query_size - query_end;
	    size_t size = aligned_query.size();
	    for(size_t i = 0; i < (size + 1)/ 2; ++i){
		    unsigned char c = TABLE[aligned_query[i]];
		    aligned_query[i] = TABLE[aligned_query[size - i - 1]];
		    aligned_query[size - i - 1] = c;
	    }
	    size = aligned_target.size();
        for(size_t i = 0; i < (size + 1)/ 2; ++i){
            unsigned char c = TABLE[aligned_target[i]];
            aligned_target[i] = TABLE[aligned_target[size - i - 1]];
            aligned_target[size - i - 1] = c;
        }
    }
    Rearrange(aligned_query, aligned_target);
    // if(aligned_query.size() > 0) TrimEnds();
}

AlignmentResult AlignmentCache::GetAlignment(idx_t qid, idx_t tid, bool sameDirection){
    AlignmentResult result;
	auto iter = cache_.find(ToKey(qid, tid));
    if(iter != cache_.end()){
        result = iter->second;
        cache_.erase(iter);
        result.Swap(sameDirection);
    }
    return result;
}

void AlignmentCache::SetAlignment(idx_t qid, idx_t tid, const AlignmentResult &result){
	cache_[ToKey(qid, tid)] = result;
}

bool AlignmentCache::HasAlignment(idx_t qid, idx_t tid) const {
	return cache_.find(ToKey(qid, tid)) != cache_.end();
}

} // end namespace ns_banded_sw
