#!/usr/bin/env python3
"""
compute_score.py — High-performance Python replacement for compute_score.pl

Computes ICU, CC, and CAI codon scores compatible with the COOL codon
optimization tool (Chung & Lee, BMC Syst Biol 2012, 6:134).

Designed for large-scale batch processing (600K+ sequences) using:
  • NumPy vectorised operations (replaces all Perl hash-based loops)
  • multiprocessing.Pool for parallel sequence scoring

Usage:
    python compute_score.py

Output format (identical to original Perl scores.txt):
    Label    IC    CC    CAI    (tab-separated, 14 decimal places)
"""

import sys
import os
import multiprocessing as mp
import time
import numpy as np

# ============================================================
#  Codon index encoding
#  Nucleotide:  A=0  C=1  G=2  T=3=U
#  codon_idx = b0*16 + b1*4 + b2   →  range [0, 63]
#  ALL_CODONS[i] gives the 3-char codon string for index i
# ============================================================

_BASES     = 'ACGT'
ALL_CODONS = [b0 + b1 + b2 for b0 in _BASES for b1 in _BASES for b2 in _BASES]
CODON2IDX  = {c: i for i, c in enumerate(ALL_CODONS)}   # 'ACG' → 17

# Fast lookup table: ord(nucleotide) → 0/1/2/3  (−1 for unknown)
_NT_LUT              = np.full(256, -1, dtype=np.int8)
for _c, _v in [('A', 0), ('C', 1), ('G', 2), ('T', 3), ('U', 3)]:
    _NT_LUT[ord(_c)] = _v
    _NT_LUT[ord(_c.lower())] = _v    # lowercase too

# Amino-acid index  (21 entries: 20 aa + stop '*')
ALL_AA  = sorted(['A','C','D','E','F','G','H','I','K','L',
                  'M','N','P','Q','R','S','T','V','W','Y','*'])
AA2IDX  = {a: i for i, a in enumerate(ALL_AA)}
N_AA    = len(ALL_AA)   # 21
STOP_I  = AA2IDX['*']
MET_I   = AA2IDX['M']


# ============================================================
#  Reference table loading  (mirrors comp_bias.pl)
# ============================================================

def _normalize_codon(s: str) -> str:
    """Upper-case, U→T, strip whitespace."""
    return s.strip().upper().replace('U', 'T')


def load_trans_table(path: str):
    """
    Parse trans_table.txt:  <codon>\t<aa>
    Returns
        c2aa : np.ndarray  shape (64,)  int8
            c2aa[i] = AA index for codon i, or -1 if not in table
    """
    c2aa = np.full(64, -1, dtype=np.int8)
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            codon = _normalize_codon(parts[0])
            aa    = parts[1].strip().upper().replace('.', '*')
            if codon not in CODON2IDX or aa not in AA2IDX:
                continue
            c2aa[CODON2IDX[codon]] = AA2IDX[aa]
    return c2aa


def build_icu_table(count_codon_path: str, count_aa_path: str,
                    c2aa: np.ndarray) -> np.ndarray:
    """
    Replicate comp_icu_table() from comp_bias.pl.

    icu_freq[c] = codon_count[c] / aa_count[ aa_of_c ]
                = 0  if aa_count is 0

    Returns icu_freq : np.ndarray shape (64,) float64
    """
    codon_count = np.zeros(64,   dtype=np.float64)
    aa_count    = np.zeros(N_AA, dtype=np.float64)

    with open(count_codon_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            codon = _normalize_codon(parts[0])
            if codon not in CODON2IDX:
                continue
            cidx = CODON2IDX[codon]
            cnt  = float(parts[1])
            codon_count[cidx] = cnt
            aidx = int(c2aa[cidx])
            if aidx >= 0:
                aa_count[aidx] += cnt

    # Vectorised division
    denom    = aa_count[c2aa.clip(0)]          # shape (64,)
    has_aa   = (c2aa >= 0) & (denom > 0)
    icu_freq = np.zeros(64, dtype=np.float64)
    icu_freq[has_aa] = codon_count[has_aa] / denom[has_aa]
    return icu_freq


def build_cc_table(count_codonp_path: str, c2aa: np.ndarray) -> np.ndarray:
    """
    Replicate comp_cc_table() from comp_bias.pl.

    cc_freq[c1*64+c2] = codonpair_count[c1,c2] / aapair_count[aa1,aa2]
    Only defined when c2aa[c1] is NOT stop (*).

    Returns cc_freq : np.ndarray shape (4096,) float64
    """
    # Pre-build: for each pair index (c1,c2) → amino-acid-pair index
    # -1 means "excluded" (c1 codes for stop, or c1/c2 unknown)
    aap_lut = np.full(4096, -1, dtype=np.int32)
    for c1 in range(64):
        aa1 = int(c2aa[c1])
        if aa1 < 0 or aa1 == STOP_I:
            continue
        for c2 in range(64):
            aa2 = int(c2aa[c2])
            if aa2 < 0:
                continue
            aap_lut[c1 * 64 + c2] = aa1 * N_AA + aa2

    codonp_count = np.zeros(4096,        dtype=np.float64)
    aap_count    = np.zeros(N_AA * N_AA, dtype=np.float64)

    with open(count_codonp_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            pair_str = _normalize_codon(parts[0])
            if len(pair_str) != 6:
                continue
            c1s, c2s = pair_str[:3], pair_str[3:]
            if c1s not in CODON2IDX or c2s not in CODON2IDX:
                continue
            pidx = CODON2IDX[c1s] * 64 + CODON2IDX[c2s]
            cnt  = float(parts[1])
            codonp_count[pidx] = cnt
            aip = int(aap_lut[pidx])
            if aip >= 0:
                aap_count[aip] += cnt

    # Vectorised computation of cc_freq
    valid_pidx   = np.where(aap_lut >= 0)[0]
    denom_vals   = aap_count[aap_lut[valid_pidx]]
    nonzero_mask = denom_vals > 0

    cc_freq = np.zeros(4096, dtype=np.float64)
    sel     = valid_pidx[nonzero_mask]
    cc_freq[sel] = codonp_count[sel] / denom_vals[nonzero_mask]

    return cc_freq, aap_lut


def build_tables(ref_dir: str) -> dict:
    """
    Load and precompute all reference tables from ref_dir.
    Returns a dict that is passed to workers (read-only).
    """
    def p(name):
        return os.path.join(ref_dir, name)

    c2aa      = load_trans_table(p('trans_table.txt'))
    icu_freq  = build_icu_table(p('count_codon.txt'), p('count_aa.txt'), c2aa)
    cc_freq, aap_lut = build_cc_table(p('count_codonp.txt'), c2aa)

    # Most-frequent ICU value per amino acid  (denominator for CAI)
    most_freq_ic_val = np.zeros(N_AA, dtype=np.float64)
    for c in range(64):
        aidx = int(c2aa[c])
        if aidx >= 0 and icu_freq[c] > most_freq_ic_val[aidx]:
            most_freq_ic_val[aidx] = icu_freq[c]

    # num_ic = scalar keys %icu_freq in Perl  (= 64, all codons present)
    num_ic = int(np.sum(c2aa >= 0))   # codons with a valid aa mapping

    # num_cc = scalar keys %cc_freq in Perl
    #  = number of (c1,c2) pairs where c1 is NOT stop  (regardless of count)
    num_cc = int(np.sum(aap_lut >= 0))

    return {
        'c2aa'            : c2aa,              # (64,)  int8
        'icu_freq'        : icu_freq,           # (64,)  float64
        'cc_freq'         : cc_freq,            # (4096,) float64
        'aap_lut'         : aap_lut,            # (4096,) int32
        'most_freq_ic_val': most_freq_ic_val,   # (21,)  float64
        'num_ic'          : num_ic,
        'num_cc'          : num_cc,
    }


# ============================================================
#  FASTA parser
# ============================================================

def parse_fasta(path: str):
    """
    Generator → (header: str, seq: str)
    seq is uppercase, U→T, whitespace stripped.
    Supports single-line and multi-line FASTA.
    """
    header    = None
    seq_parts = []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\r\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header    = line[1:]
                seq_parts = []
            else:
                seq_parts.append(
                    line.upper().replace('U', 'T')
                        .replace(' ', '').replace('\t', '')
                )
    if header is not None:
        yield header, ''.join(seq_parts)


# ============================================================
#  Sequence encoder
# ============================================================

def seq_to_codon_idx(seq: str, n3: int) -> np.ndarray:
    """
    Convert first n3 characters of seq (n3 must be divisible by 3)
    to an int16 array of codon indices 0-63.
    Unknown nucleotides produce index -1.
    """
    raw  = np.frombuffer(seq[:n3].encode('ascii', errors='replace'), dtype=np.uint8)
    nts  = _NT_LUT[raw]                          # (n3,)  int8, -1 for unknown
    arr3 = nts.reshape(-1, 3)                    # (n3//3, 3)
    cidx = (arr3[:, 0].astype(np.int16) * 16
            + arr3[:, 1].astype(np.int16) * 4
            + arr3[:, 2].astype(np.int16))       # codon indices
    # Mark codons with any unknown NT as -1
    bad  = np.any(arr3 < 0, axis=1)
    cidx[bad] = -1
    return cidx


# ============================================================
#  Core scoring function
#  (exact replication of process_seq in compute_score.pl)
# ============================================================

def compute_scores(header: str, seq: str, T: dict):
    """
    Compute (header, ICU, CC, CAI) for one sequence.
    T = tables dict from build_tables().
    Returns None if the sequence is too short or empty.

    Mirrors Perl functions:
        comp_ic_counts, comp_cc_counts, comp_icu, comp_cc, comp_cai
    """
    c2aa             = T['c2aa']
    icu_freq         = T['icu_freq']
    cc_freq          = T['cc_freq']
    aap_lut          = T['aap_lut']
    most_freq_ic_val = T['most_freq_ic_val']
    num_ic           = T['num_ic']
    num_cc           = T['num_cc']

    # ── 1. trim to multiple of 3 ──────────────────────────────────────
    n = len(seq)
    n -= n % 3
    if n < 6:
        return None        # Perl: "too short to be considered"

    # ── 2. encode codons ───────────────────────────────────────────────
    cidxs = seq_to_codon_idx(seq, n)        # int16, -1 for unknown
    nc    = len(cidxs)
    if nc < 2:
        return None

    aa_arr = np.where(cidxs >= 0, c2aa[cidxs.clip(0)].astype(np.int16), -1)

    # ── 3. strip multiple trailing stops ──────────────────────────────
    # Perl: if seq ends with ≥2 stop codons, remove all but last
    i_last = nc - 1
    while i_last >= 0 and int(aa_arr[i_last]) == STOP_I:
        i_last -= 1
    n_trailing_stops = nc - 1 - i_last       # codons after last non-stop codon
    if n_trailing_stops > 1:
        # keep non-stop codons + exactly one stop, matching Perl:
        # substr($seq, 0, 3 - length($1))  where $1 = all trailing stops
        keep = nc - (n_trailing_stops - 1)
        cidxs  = cidxs[:keep]
        aa_arr = aa_arr[:keep]
        nc = keep
    if nc < 2:
        return None

    # ── 4. ICU counts ─────────────────────────────────────────────────
    # Perl: first codon counted only if it encodes Met (M)
    first_aa = int(aa_arr[0])
    ic_start = 0 if first_aa == MET_I else 1

    ic_cidxs = cidxs[ic_start:]
    ic_valid  = ic_cidxs[ic_cidxs >= 0]           # remove unknown codons
    ic_aa     = aa_arr[ic_start:]
    ic_aa_val = ic_aa[(ic_cidxs >= 0) & (ic_aa >= 0)]

    c_counts  = np.bincount(ic_valid.astype(np.int64), minlength=64).astype(np.float64)
    aa_counts = np.bincount(ic_aa_val.astype(np.int64), minlength=N_AA).astype(np.float64)

    # ── 5. ICU score  (comp_icu) ──────────────────────────────────────
    # For each codon c (64 total):
    #   if aa_counts[aa_of_c] > 0:
    #       score += |c_counts[c] / aa_counts[aa_of_c]  −  icu_freq[c]|
    #   else:
    #       score += icu_freq[c]
    aa_denom  = aa_counts[c2aa.clip(0)]   # aa_counts for aa of each codon
    valid_ic  = (c2aa >= 0) & (aa_denom > 0)

    icu_terms         = np.empty(64, dtype=np.float64)
    icu_terms[ valid_ic] = np.abs(c_counts[valid_ic] / aa_denom[valid_ic]
                                  - icu_freq[valid_ic])
    icu_terms[~valid_ic]  = np.where(c2aa[~valid_ic] >= 0, icu_freq[~valid_ic], 0.0)
    score_icu = -np.sum(icu_terms) / num_ic

    # ── 6. CC counts  (comp_cc_counts) ───────────────────────────────
    # Perl logic (index i over @codons_array[1:end-1]):
    #   - include (codons_array[i], codons_array[i+1]) if:
    #       both indices exist (valid codon)
    #       aa[i] != stop
    #       aa[i+1] != stop  OR  i+1 == last codon index
    #   - First pair treated same as ICU: only included if first codon is Met

    c1_all = cidxs[:-1]           # int16
    c2_all = cidxs[1:]            # int16
    aa1    = aa_arr[:-1]
    aa2    = aa_arr[1:]
    n_pairs = len(c1_all)

    c1_ok  = c1_all >= 0
    c2_ok  = c2_all >= 0
    aa1_ok = (aa1 >= 0) & (aa1 != STOP_I)
    aa2_ok = aa2 >= 0
    # c2 may be stop only on the very last pair (Perl: i+1 < $#codons_array excludes stop)
    aa2_stop_ok = (aa2 == STOP_I) & (np.arange(n_pairs) == n_pairs - 1)
    aa2_valid   = aa2_ok & ((aa2 != STOP_I) | aa2_stop_ok)

    pair_mask = c1_ok & c2_ok & aa1_ok & aa2_valid

    # Met rule for first pair
    if n_pairs > 0 and first_aa != MET_I:
        pair_mask[0] = False

    p_c1 = c1_all[pair_mask].astype(np.int64)
    p_c2 = c2_all[pair_mask].astype(np.int64)
    p_aa1 = aa1[pair_mask].astype(np.int64)
    p_aa2 = aa2[pair_mask].astype(np.int64)

    pair_idx_arr = p_c1 * 64 + p_c2
    cc_counts    = np.bincount(pair_idx_arr, minlength=4096).astype(np.float64)

    aap_idx_arr  = p_aa1 * N_AA + p_aa2
    aap_counts   = np.bincount(aap_idx_arr, minlength=N_AA * N_AA).astype(np.float64)

    # ── 7. CC score  (comp_cc) ─────────────────────────────────────────
    # For each pair p in cc table (where aap_lut[p] >= 0):
    #   aap = aap_lut[p]
    #   if aap_counts[aap] > 0:
    #       score += |cc_counts[p] / aap_counts[aap]  −  cc_freq[p]|
    #   else:
    #       score += cc_freq[p]
    valid_p   = aap_lut >= 0                          # bool (4096,)
    aap_here  = aap_counts[aap_lut.clip(0)]           # (4096,)

    has_aap   = valid_p & (aap_here > 0)
    cc_terms  = np.empty(4096, dtype=np.float64)
    cc_terms[has_aap]   = np.abs(cc_counts[has_aap] / aap_here[has_aap]
                                 - cc_freq[has_aap])
    cc_terms[~has_aap]  = np.where(valid_p[~has_aap], cc_freq[~has_aap], 0.0)
    score_cc  = -np.sum(cc_terms) / num_cc

    # ── 8. CAI score  (comp_cai) ──────────────────────────────────────
    # Perl: iterate all codons in full sequence (not just from ic_start)
    all_valid_c = cidxs[cidxs >= 0]
    all_aa_v    = c2aa[all_valid_c.clip(0)]
    cai_mask    = all_aa_v >= 0

    cai_c   = all_valid_c[cai_mask]
    cai_aa  = all_aa_v[cai_mask].astype(np.int64)

    mfv     = most_freq_ic_val[cai_aa]      # shape (m,)
    icu_v   = icu_freq[cai_c]               # shape (m,)

    # Guard: avoid log(0)
    safe_icu = np.where(icu_v  > 0, icu_v,  1e-300)
    safe_mfv = np.where(mfv    > 0, mfv,    1.0)

    log_sum = np.sum(np.log(safe_icu / safe_mfv))
    m       = len(cai_c)
    score_cai = float(np.exp(log_sum / m)) if m > 0 else 0.0

    return header, float(score_icu), float(score_cc), float(score_cai)


# ============================================================
#  Multiprocessing helpers
# ============================================================

# Global tables reference (set in each worker via initializer)
_WORKER_TABLES = None


def _worker_init(tables_dict: dict):
    global _WORKER_TABLES
    _WORKER_TABLES = tables_dict


def _worker_score(record: tuple):
    """Called in a worker process for one (header, seq) record."""
    header, seq = record
    try:
        result = compute_scores(header, seq, _WORKER_TABLES)
        if result is None:
            return None
        h, icu, cc, cai = result
        return f"{h}\t{icu:.14f}\t{cc:.14f}\t{cai:.14f}"
    except Exception as exc:
        return f"# ERROR [{header}]: {exc}"


# ============================================================
#  Entry point
# ============================================================

def main():
    # --- Parameters ---
    file_path = 'sequence(DNA).txt'
    output_file = 'scores.txt'
    
    # Using default values instead of arguments
    jobs = 0        # 0 = all CPU cores
    chunk = 300     # tuning parameter
    ref_dir_input = None
    quiet = False

    # ── locate reference files ───────────────────────────────────────
    if ref_dir_input:
        ref_dir = os.path.abspath(ref_dir_input)
    else:
        ref_dir = os.path.dirname(os.path.abspath(__file__))

    if not quiet:
        print(f"[INFO] Reference dir : {ref_dir}")
        print(f"[INFO] Input         : {file_path}")
        print(f"[INFO] Output        : {output_file}")

    # ── load reference tables ────────────────────────────────────────
    if not quiet:
        print("[INFO] Loading reference tables...")
    T0 = time.perf_counter()
    tables = build_tables(ref_dir)
    if not quiet:
        print(f"[INFO] num_ic={tables['num_ic']}  num_cc={tables['num_cc']}  "
              f"(table load: {time.perf_counter()-T0:.2f}s)")

    # ── parse FASTA ──────────────────────────────────────────────────
    if not quiet:
        print("[INFO] Parsing FASTA sequences...")
    records = list(parse_fasta(file_path))
    total   = len(records)
    if not quiet:
        print(f"[INFO] {total:,} sequences loaded.")

    if total == 0:
        print("[WARN] No sequences found in input file. Exiting.")
        return

    # ── determine parallelism ────────────────────────────────────────
    n_jobs = jobs if jobs > 0 else mp.cpu_count()
    n_jobs = max(1, min(n_jobs, total))

    if not quiet:
        print(f"[INFO] Workers: {n_jobs}   Chunk size: {chunk}")

    # ── process & write ──────────────────────────────────────────────
    t_start = time.perf_counter()

    with open(output_file, 'w', buffering=1 << 20) as out_fh:
        out_fh.write("Label\tIC\tCC\tCAI\n")
        written = 0
        done    = 0

        def _report(done, rate=None):
            if quiet:
                return
            if rate and rate > 0:
                eta = (total - done) / rate
                print(f"\r[PROG] {done:>9,}/{total:,}  "
                      f"({done/total*100:5.1f}%)  "
                      f"{rate:>8,.0f} seq/s  "
                      f"ETA {int(eta//60):02d}:{int(eta%60):02d}",
                      end='', flush=True)
            else:
                print(f"\r[PROG] {done:>9,}/{total:,}  "
                      f"({done/total*100:5.1f}%)",
                      end='', flush=True)

        if n_jobs == 1:
            # ── single-process (easier porting/debugging) ────────────
            for header, seq in records:
                result = compute_scores(header, seq, tables)
                if result is not None:
                    h, icu, cc, cai = result
                    out_fh.write(f"{h}\t{icu:.14f}\t{cc:.14f}\t{cai:.14f}\n")
                    written += 1
                done += 1
                if done % 5000 == 0:
                    elapsed = time.perf_counter() - t_start
                    _report(done, done / elapsed if elapsed > 0 else None)
        else:
            # ── multi-process ────────────────────────────────────────
            with mp.Pool(
                processes=n_jobs,
                initializer=_worker_init,
                initargs=(tables,),
            ) as pool:
                for line in pool.imap(_worker_score, records, chunksize=chunk):
                    if line is not None:
                        out_fh.write(line + '\n')
                        written += 1
                    done += 1
                    if done % 5000 == 0:
                        elapsed = time.perf_counter() - t_start
                        _report(done, done / elapsed if elapsed > 0 else None)

    elapsed = time.perf_counter() - t_start
    rate    = total / elapsed if elapsed > 0 else float('inf')
    if not quiet:
        print(f"\n[DONE] {total:,} sequences processed in "
              f"{elapsed:.1f}s  ({rate:,.0f} seq/s)")
        print(f"[DONE] {written:,} results written → {output_file}")


if __name__ == '__main__':
    mp.freeze_support()     # Required for Windows .exe via PyInstaller
    main()
