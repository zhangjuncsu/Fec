/*******************************************************************************************
 *
 *  Local alignment module.  Routines for finding local alignments given a seed position,
 *    representing such an l.a. with its interval and a set of pass-thru points, so that
 *    a detailed alignment can be efficiently computed on demand.
 *
 *  All routines work on a numeric representation of DNA sequences, i.e. 0 for A, 1 for C,
 *    2 for G, and 3 for T.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#ifndef _A_MODULE

#define _A_MODULE

///#include "DB.h"

#ifdef __cplusplus
extern "C" {
#endif 

#include <stdint.h>
#include <stdio.h>

#include "kstring.h"

typedef int8_t		int8;
typedef uint8_t		uint8;
typedef int16_t 	int16;
typedef uint16_t	uint16;
typedef uint32_t 	uint32;
typedef int64_t 	int64;
typedef uint64_t	uint64;

#define TRACE_XOVR 125   //  If the trace spacing is not more than this value, then can
                         //    and do compress traces pts to 8-bit unsigned ints

/*** INTERACTIVE vs BATCH version
     The defined constant INTERACTIVE (set in DB.h) determines whether an interactive or
       batch version of the routines in this library are compiled.  In batch mode, routines
       print an error message and exit.  In interactive mode, the routines place the error
       message in EPLACE (also defined in DB.h) and return an error value, typically NULL
       if the routine returns a pointer, and an unusual integer value if the routine returns
       an integer.
     Below when an error return is described, one should understand that this value is returned
       only if the routine was compiled in INTERACTIVE mode.
***/


/*** PATH ABSTRACTION:
     Coordinates are *between* characters where 0 is the tick just before the first char,
     1 is the tick between the first and second character, and so on.  Our data structure
     is called a Path refering to its conceptualization in an edit graph.
     A local alignment is specified by the point '(abpos,bbpos)' at which its path in
     the underlying edit graph starts, and the point '(aepos,bepos)' at which it ends.
     In otherwords A[abpos+1..aepos] is aligned to B[bbpos+1..bepos] (assuming X[1] is
     the *first* character of X).
     There are 'diffs' differences in an optimal local alignment between the beginning and
     end points of the alignment (if computed by Compute_Trace), or nearly so (if computed
     by Local_Alignment).  
     Optionally, a Path can have additional information about the exact nature of the
     aligned substrings if the field 'trace' is not NULL.  Trace points to either an
     array of integers (if computed by a Compute_Trace routine), or an array of unsigned
     short integers (if computed by Local_Alignment).
     If computed by Local_Alignment 'trace' points at a list of 'tlen' (always even) short
     values:
            d_0, b_0, d_1, b_1, ... d_n-1, b_n-1, d_n, b_n
     to be interpreted as follows.  The alignment from (abpos,bbpos) to (aepos,bepos)
     passes through the n trace points for i in [1,n]:
            (a_i,b_i) where a_i = floor(abpos/TS)*TS + i*TS
                        and b_i = bbpos + (b_0 + b_1 + b_i-1)
     where also let a_0,b_0 = abpos,bbpos and a_(n+1),b_(n+1) = aepos,bepos.  That is, the
     interior (i.e. i != 0 and i != n+1) trace points pass through every TS'th position of
     the aread where TS is the "trace spacing" employed when finding the alignment (see
     New_Align_Spec).  Typically TS is 100.  Then d_i is the number of differences in the
     portion of the alignment between (a_i,b_i) and (a_i+1,b_i+1).  These trace points allow
     the Compute_Trace routines to efficiently compute the exact alignment between the two
     reads by efficiently computing exact alignments between consecutive pairs of trace points.
     Moreover, the diff values give one an idea of the quality of the alignment along every
     segment of TS symbols of the aread.
     If computed by a Compute_Trace routine, 'trace' points at a list of 'tlen' integers
     < i1, i2, ... in > that encodes an exact alignment as follows.  A negative number j
     indicates that a dash should be placed before A[-j] and a positive number k indicates
     that a dash should be placed before B[k], where A and B are the two sequences of the
     overlap.  The indels occur in the trace in the order in which they occur along the
     alignment.  For a good example of how to "decode" a trace into an alignment, see the
     code for the routine Print_Alignment.
***/

typedef struct
  { void     *trace;
    int       tlen;
    int       diffs;
    int       abpos, bbpos;
    int       aepos, bepos;
  } Path;


/*** ALIGNMENT ABSTRACTION:
     An alignment is modeled by an Alignment record, which in addition to a *pointer* to a
     'path', gives pointers to the A and B sequences, their lengths, and indicates whether
     the B-sequence needs to be complemented ('comp' non-zero if so).  The 'trace' pointer
     of the 'path' subrecord can be either NULL, a list of pass-through points, or an exact
     trace depending on what routines have been called on the record.
     One can (1) compute a trace, with Compute_Trace, either from scratch if 'path.trace' = NULL,
     or using the sequence of pass-through points in trace, (2) print an ASCII representation
     of an alignment, or (3) reverse the roles of A and B, and (4) complement a sequence
     (which is a reversible process).
     If the alignment record shows the B sequence as complemented, *** THEN IT IS THE
     RESPONSIBILITY OF THE CALLER *** to make sure that bseq points at a complement of
     the sequence before calling Compute_Trace or Print_Alignment.  Complement_Seq complements
     the sequence a of length n.  The operation does the complementation/reversal in place.
     Calling it a second time on a given fragment restores it to its original state.
     With the introduction of the DAMAPPER, we need to code chains of alignments between a
     pair of sequences.  The alignments of a chain are expected to be found in order either on
     a file or in memory, where the START_FLAG marks the first alignment and the NEXT_FLAG all
     subsequent alignmenst in a chain.  A chain of a single LA is marked with the START_FLAG.
     The BEST_FLAG marks one of the best chains for a pair of sequences.  The convention is
     that either every record has either a START- or NEXT-flag, or none of them do (e.g. as
     produced by daligner), so one can always check the flags of the first alignment to see
     whether or not the chain concept applies to a given collection or not.
***/

#define COMP_FLAG  0x1
#define ACOMP_FLAG 0x2   //  A-sequence is complemented, not B !  Only Local_Alignment notices

#define COMP(x)   ((x) & COMP_FLAG)
#define ACOMP(x)  ((x) & ACOMP_FLAG)

#define START_FLAG 0x4   //  LA is the first of a chain of 1 or more la's
#define NEXT_FLAG  0x8   //  LA is the next segment of a chain.
#define BEST_FLAG  0x10  //  This is the start of the best chain

#define CHAIN_START(x)  ((x) & START_FLAG)
#define CHAIN_NEXT(x)   ((x) & NEXT_FLAG)
#define BEST_CHAIN(x)   ((x) & BEST_FLAG)

typedef struct
  { Path   *path;
    uint32  flags;        /* Pipeline status and complementation flags          */
    char   *aseq;         /* Pointer to A sequence                              */
    char   *bseq;         /* Pointer to B sequence                              */
    int     alen;         /* Length of A sequence                               */
    int     blen;         /* Length of B sequence                               */
  } Alignment;

void Complement_Seq(char *a, int n);

  /* Many routines like Local_Alignment, Compute_Trace, and Print_Alignment need working
     storage that is more efficiently reused with each call, rather than being allocated anew
     with each call.  Each *thread* can create a Work_Data object with New_Work_Data and this
     object holds and retains the working storage for routines of this module between calls
     to the routines.  If enough memory for a Work_Data is not available then NULL is returned.
     Free_Work_Data frees a Work_Data object and all working storage held by it.
  */

  typedef void Work_Data;

  Work_Data *New_Work_Data();

  void       Free_Work_Data(Work_Data *work);

  /* Local_Alignment seeks local alignments of a quality determined by a number of parameters.
     These are coded in an Align_Spec object that can be created with New_Align_Spec and
     freed with Free_Align_Spec when no longer needed.  There are 4 essential parameters:
     ave_corr:    the average correlation (1 - 2*error_rate) for the sought alignments.  For Pacbio
                    data we set this to .70 assuming an average of 15% error in each read.
     trace_space: the spacing interval for keeping trace points and segment differences (see
                    description of 'trace' for Paths above)
     freq[4]:     a 4-element vector where afreq[0] = frequency of A, f(A), freq[1] = f(C),
                    freq[2] = f(G), and freq[3] = f(T).  This vector is part of the header
                    of every HITS database (see db.h).
     If an alignment cannot reach the boundary of the d.p. matrix with this condition (i.e.
     overlap), then the last/first 30 columns of the alignment are guaranteed to be
     suffix/prefix positive at correlation ave_corr * g(freq) where g is an empirically
     measured function that increases from 1 as the entropy of freq decreases.  If memory is
     unavailable or the freq distribution is too skewed then NULL is returned.
     You can get back the original parameters used to create an Align_Spec with the simple
     utility functions below.
  */

  typedef void Align_Spec;

  Align_Spec *New_Align_Spec(double ave_corr, int trace_space, float *freq);

  void        Free_Align_Spec(Align_Spec *spec);

  int    Trace_Spacing      (Align_Spec *spec);
  double Average_Correlation(Align_Spec *spec);
  float *Base_Frequencies   (Align_Spec *spec);

  /* Local_Alignment finds the longest significant local alignment between the sequences in
     'align' subject to:
       (a) the alignment criterion given by the Align_Spec 'spec',
       (b) it passes through one of the points (anti+k)/2,(anti-k)/2 for k in [low,hgh] within
             the underlying dynamic programming matrix (i.e. the points on diagonals low to hgh
             on anti-diagonal anti or anti-1 (depending on whether the diagonal is odd or even)),
       (c) if lbord >= 0, then the alignment is always above diagonal low-lbord, and
       (d) if hbord >= 0, then the alignment is always below diagonal hgh+hbord.
     The path record of 'align' has its 'trace' filled from the point of view of an overlap
     between the aread and the bread.  In addition a Path record from the point of view of the
     bread versus the aread is returned by the function, with this Path's 'trace' filled in
     appropriately.  The space for the returned path and the two 'trace's are in the working
     storage supplied by the Work_Data packet and this space is reused with each call, so if
     one wants to retain the bread-path and the two trace point sequences, then they must be
     copied to user-allocated storage before calling the routine again.  NULL is returned in
     the event of an error.
     Find_Extension is a variant of Local_Alignment that simply finds a local alignment that
     either ends (if prefix is non-zero) or begins (if prefix is zero) at the point
     (anti+diag)/2,(anti-diag)/2).  All other parameters are as before.  It returns a non-zero
     value only when INTERACTIVE is on and it cannot allocate the memory it needs.
     Only the path and trace with respect to the aread is returned.  This routine is experimental
     and may not persist in later versions of the code.
  */

  Path *Local_Alignment(Alignment *align, Work_Data *work, Align_Spec *spec,
                        int low, int hgh, int anti, int lbord, int hbord);

  int   Find_Extension(Alignment *align, Work_Data *work, Align_Spec *spec,    //  experimental !!
                       int diag, int anti, int lbord, int hbord, int prefix);

  /* Given a legitimate Alignment object, Compute_Trace_X computes an exact trace for the alignment.
     If 'path.trace' is non-NULL, then it is assumed to be a sequence of pass-through points
     and diff levels computed by Local_Alignment.  In either case 'path.trace' is set
     to point at an integer array within the storage of the Work_Data packet encoding an
     exact optimal trace from the start to end points.  If the trace is needed beyond the
     next call to a routine that sets it, then it should be copied to an array allocated
     and managed by the caller.
     Compute_Trace_ALL does not require a sequence of pass-through points, as it computes the
     best alignment between (path->abpos,path->bbpos) and (path->aepos,path->bepos) in the
     edit graph between the sequences.  Compute_Trace_PTS computes a trace by computing the
     trace between successive pass through points.  It is much, much faster than Compute_Trace_ALL
     but at the tradeoff of not necessarily being optimal as pass-through points are not all
     perfect.  Compute_Trace_MID computes a trace by computing the trace between the mid-points
     of alignments between two adjacent pairs of pass through points.  It is generally twice as
     slow as Compute_Trace_PTS, but it produces nearer optimal alignments.  All these routines
     return 1 if an error occurred and 0 otherwise.
  */

#define LOWERMOST -1   //   Possible modes for "mode" parameter below)
#define GREEDIEST  0
#define UPPERMOST  1

  int Compute_Trace_ALL(Alignment *align, Work_Data *work);
  int Compute_Trace_PTS(Alignment *align, Work_Data *work, int trace_spacing, int mode);
  int Compute_Trace_MID(Alignment *align, Work_Data *work, int trace_spacing, int mode);

  /* Compute_Trace_IRR (IRR for IRRegular) computes a trace for the given alignment where
     it assumes the spacing between trace points between both the A and B read varies, and
     futher assumes that the A-spacing is given in the short integers normally occupied by
     the differences in the alignment between the trace points.  This routine is experimental
     and may not persist in later versions of the code.
  */

  int Compute_Trace_IRR(Alignment *align, Work_Data *work, int mode);   //  experimental !!

  /* Alignment_Cartoon prints an ASCII representation of the overlap relationhip between the
     two reads of 'align' to the given 'file' indented by 'indent' space.  Coord controls
     the display width of numbers, it must be not less than the width of any number to be
     displayed.
     If the alignment trace is an exact trace, then one can ask Print_Alignment to print an
     ASCII representation of the alignment 'align' to the file 'file'.  Indent the display
     by "indent" spaces and put "width" columns per line in the display.  Show "border"
     characters of sequence on each side of the aligned region.  If upper is non-zero then
     display bases in upper case.  If coord is greater than 0, then the positions of the
     first character in A and B in the given row is displayed with a field width given by
     coord's value.
     Print_Reference is like Print_Alignment but rather than printing exaclty "width" columns
     per segment, it prints "block" characters of the A sequence in each segment.  This results
     in segments of different lengths, but is convenient when looking at two alignments involving
     A as segments are guaranteed to cover the same interval of A in a segment.
     Both Print routines return 1 if an error occurred (not enough memory), and 0 otherwise.
     Flip_Alignment modifies align so the roles of A and B are reversed.  If full is off then
     the trace is ignored, otherwise the trace must be to a full alignment trace and this trace
     is also appropriately inverted.
  */

  void Alignment_Cartoon(FILE *file, Alignment *align, int indent, int coord);

  int  Print_Alignment(FILE *file, Alignment *align, Work_Data *work,
                       int indent, int width, int border, int upper, int coord);

  int  Print_Reference(FILE *file, Alignment *align, Work_Data *work,
                       int indent, int block, int border, int upper, int coord);

  void Flip_Alignment(Alignment *align, int full);
  
  int Get_Alignment(Alignment *align, Work_Data *ework, kstring_t* qaln, kstring_t* taln);
  
  void
  Print_Align_String(FILE* file, kstring_t* qaln, kstring_t* taln);


/*** OVERLAP ABSTRACTION:
     Externally, between modules an Alignment is modeled by an "Overlap" record, which
     (a) replaces the pointers to the two sequences with their ID's in the HITS data bases,
     (b) does not contain the length of the 2 sequences (must fetch from DB), and
     (c) contains its path as a subrecord rather than as a pointer (indeed, typically the
     corresponding Alignment record points at the Overlap's path sub-record).  The trace pointer
     is always to a sequence of trace points and can be either compressed (uint8) or
     uncompressed (uint16).  One can read and write binary records of an "Overlap".
***/

typedef struct {
  Path    path;         /* Path: begin- and end-point of alignment + diffs    */
  uint32  flags;        /* Pipeline status and complementation flags          */
  int     aread;        /* Id # of A sequence                                 */
  int     bread;        /* Id # of B sequence                                 */
} Overlaps;


  /* Read_Overlap reads the next Overlap record from stream 'input', not including the trace
     (if any), and without modifying 'ovl's trace pointer.  Read_Trace reads the ensuing trace
     into the memory pointed at by the trace field of 'ovl'.  It is assumed to be big enough to
     accommodate the trace where each value take 'tbytes' bytes (1 if uint8 or 2 if uint16).
     Write_Overlap write 'ovl' to stream 'output' followed by its trace vector (if any) that
     occupies 'tbytes' bytes per value.  It returns non-zero if there was an error writing.
     Print_Overlap prints an ASCII version of the contents of 'ovl' to stream 'output'
     where the trace occupes 'tbytes' per value and the print out is indented from the left
     margin by 'indent' spaces.
     Compress_TraceTo8 converts a trace fo 16-bit values to 8-bit values in place, and
     Decompress_TraceTo16 does the reverse conversion.
     Check_Trace_Points checks that the number of trace points is correct and that the sum
     of the b-read displacements equals the b-read alignment interval, assuming the trace
     spacing is 'tspace'.  It reports an error message if there is a problem and 'verbose'
     is non-zero.  The 'ovl' came from the file names 'fname'.
  */

  int Read_Overlap(FILE *input, Overlaps *ovl);
  int Read_Trace(FILE *innput, Overlaps *ovl, int tbytes);

  int  Write_Overlap(FILE *output, Overlaps *ovl, int tbytes);
  void Print_Overlap(FILE *output, Overlaps *ovl, int tbytes, int indent);

  void Compress_TraceTo8(Overlaps *ovl);
  void Decompress_TraceTo16(Overlaps *ovl);

  int  Check_Trace_Points(Overlaps *ovl, int tspace, int verbose, char *fname);
  
#ifdef __cplusplus
}
#endif

#endif // _A_MODULE