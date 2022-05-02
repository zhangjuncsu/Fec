#include "Fec_extension.h"

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
			  bool recuse_long_indel)
{
	const BOOL small_edlib = onc_align(query,
					   query_start,
					   query_size,
					   target,
					   target_start,
					   target_size,
					   align_data,
					   kOcaBlockSize,
					   min_align_size,
                       ONC_TAIL_MATCH_LEN_LONG);
	BOOL r = small_edlib;
	if (r) {
		int qbeg = oca_query_start(*align_data);
		int qend = oca_query_end(*align_data);
		int lhang = qbeg;
		int rhang = query_size - qend;
		if (lhang + rhang > 500) r = FALSE;
	}
	if (r) {
		m5qoff(m5) = oca_query_start(*align_data);
		m5qend(m5) = oca_query_end(*align_data);
		m5soff(m5) = oca_target_start(*align_data);
		m5send(m5) = oca_target_end(*align_data);
		strncpy(m5qaln(m5), align_data->query_align.s, align_data->query_align.l);
		strncpy(m5saln(m5), align_data->target_align.s, align_data->target_align.l);
		m5ident(m5) = oca_ident_perc(*align_data);
		return TRUE;
	}
	
	if(recuse_long_indel) {
		const BOOL dalign = ocda_go(query, 
									query_start,
									query_size,
									target,
									target_start,
									target_size,
									dalign_data,
									min_align_size);
		if (dalign) {
			const BOOL large_edlib = edlib_go(query,
												ocda_query_start(*dalign_data),
												ocda_query_end(*dalign_data),
												target,
												ocda_target_start(*dalign_data),
												ocda_target_end(*dalign_data),
												falign_data,
												ocda_distance(*dalign_data),
												min_align_size,
												TRUE,
												4);
			if (large_edlib) {
				m5qoff(m5) = falign_data->qoff;
				m5qend(m5) = falign_data->qend;
				m5soff(m5) = falign_data->toff;
				m5send(m5) = falign_data->tend;
				strncpy(m5qaln(m5), falign_data->query_align.s, falign_data->query_align.l);
				strncpy(m5saln(m5), falign_data->target_align.s, falign_data->target_align.l);
				m5ident(m5) = falign_data->ident_perc;
				return TRUE;
			}
		}
	}

	if (small_edlib) {
		m5qoff(m5) = oca_query_start(*align_data);
		m5qend(m5) = oca_query_end(*align_data);
		m5soff(m5) = oca_target_start(*align_data);
		m5send(m5) = oca_target_end(*align_data);
		strncpy(m5qaln(m5), align_data->query_align.s, align_data->query_align.l);
		strncpy(m5saln(m5), align_data->target_align.s, align_data->target_align.l);
		m5ident(m5) = oca_ident_perc(*align_data);
		return TRUE;
	}
	
	return FALSE;
}