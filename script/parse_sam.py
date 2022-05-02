import re
import argparse

def parse_cigar(cigar: str, ref_pos: int, cov: list):
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    match, mismatch, insertion, deletion, clip = 0, 0, 0, 0, 0
    index = ref_pos - 1
    for item in pattern.findall(cigar):
        size, op = int(item[0]), item[1]
        if op == 'M' or op == '=':
            match += size
            for i in range(size):
                cov[index + i] += 1
            index += size
        elif op == 'X':
            mismatch += size
            for i in range(size):
                cov[index + i] += 1
            index += size
        elif op == 'I':
            insertion += size
        elif op == 'D':
            deletion += size
            index += size
        elif op == 'S' or op == 'H':
            clip += size
    return match, mismatch, insertion, deletion, clip, cov

def parse_sam(sam_file: str):
    match, mismatch, insertion, deletion, clip = 0, 0, 0, 0, 0
    coverages = {}
    processed = set()
    for line in open(sam_file):
        if '@SQ' in line:
            items = line.split('\t')
            name = items[1].strip()[3:]
            length = int(items[2].strip()[3:])
            coverages[name] = [0] * length
        if line.startswith('@'):
            continue
        items = line.strip().split('\t')
        read_name = items[0].strip()
        if read_name not in processed:
            processed.add(read_name)
        else:
            continue
        flag = int(items[1].strip())
        if flag & 0x100:
            continue
        ref_name = items[2].strip()
        if ref_name != '*':
            assert (ref_name in coverages), 'the length of reference {} is not record'.format(ref_name)
        ref_pos = int(items[3].strip())
        # mapQ = int(items[4].strip())
        cigar = items[5].strip()
        if flag & 0x4 and ref_name == '*':
            continue
        M, X, I, D, C, coverages[ref_name] = parse_cigar(cigar, ref_pos, coverages[ref_name])
        match += M
        mismatch += X
        insertion += I
        deletion += D
        clip += C
    return match, mismatch, insertion, deletion, clip, coverages

def stat_coverage(coverages: dict):
    ref_length = 0
    cov = 0
    for _, coverage in coverages.items():
        ref_length += len(coverage)
        cov += sum(c > 0 for c in coverage)
    fraction = cov / ref_length if ref_length > 0 else 0
    return fraction

def main():
    parser = argparse.ArgumentParser('collect statistics in sam file')
    parser.add_argument('-f', help='sam file')

    args = parser.parse_args()

    match, mismatch, insertion, deletion, clip, coverages = parse_sam(args.f)
    fraction = stat_coverage(coverages)
    print('match\tmismatch\tinsertion\tdeletion\tclip\tfraction\tidentity')
    print('{}\t{}\t{}\t{}\t{}\t{:.6f}\t{:.6f}'.format(match, mismatch, insertion, deletion, clip, fraction, match / (match + mismatch + insertion)))

if __name__ == '__main__':
    main()