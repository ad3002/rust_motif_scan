//! rust_motif_scan — fast, zero-dependency IUPAC degenerate-motif scanner.
//!
//! Scans a FASTA for an IUPAC motif on BOTH strands and writes BED. With a
//! `--canonical` reference (or the built-in `--cenpb` preset) it additionally
//! scores each hit against a reference motif by Hamming distance, decomposed
//! into left-arm / spacer / right-arm regions.
//!
//! Zero dependencies: builds with `cargo build --release` or even plain
//! `rustc -O src/main.rs`. Multithreaded across FASTA records (std::thread).
//!
//! See README.md for the full column spec and the CENP-B box rationale.

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

// ---- IUPAC / base masks (A=1,C=2,G=4,T=8) ----------------------------------

/// IUPAC pattern char -> allowed-base mask.
fn iupac_mask(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1, b'C' => 2, b'G' => 4, b'T' => 8, b'U' => 8,
        b'R' => 1 | 4, b'Y' => 2 | 8, b'S' => 2 | 4, b'W' => 1 | 8,
        b'K' => 4 | 8, b'M' => 1 | 2,
        b'B' => 2 | 4 | 8, b'D' => 1 | 4 | 8, b'H' => 1 | 2 | 8, b'V' => 1 | 2 | 4,
        b'N' => 1 | 2 | 4 | 8,
        _ => 0,
    }
}

/// Genome base -> single-bit mask; anything not exactly A/C/G/T => 0 (no match).
#[inline]
fn genome_mask(b: u8) -> u8 {
    match b { b'A' | b'a' => 1, b'C' | b'c' => 2, b'G' | b'g' => 4, b'T' | b't' => 8, _ => 0 }
}

#[inline]
fn comp_mask(m: u8) -> u8 {
    ((m & 1) << 3) | ((m & 8) >> 3) | ((m & 2) << 1) | ((m & 4) >> 1)
}

#[inline]
fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|b| match b.to_ascii_uppercase() {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', x => x,
    }).collect()
}

// ---- FASTA -----------------------------------------------------------------

struct Record { name: String, seq: Vec<u8> }

fn read_fasta(path: &str) -> std::io::Result<Vec<Record>> {
    let f = File::open(path)?;
    let mut rdr = BufReader::with_capacity(1 << 20, f);
    let mut recs = Vec::new();
    let mut line = Vec::new();
    let mut cur: Option<Record> = None;
    loop {
        line.clear();
        if rdr.read_until(b'\n', &mut line)? == 0 { break; }
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) { line.pop(); }
        if line.first() == Some(&b'>') {
            if let Some(r) = cur.take() { recs.push(r); }
            let hdr = &line[1..];
            let end = hdr.iter().position(|c| *c == b' ' || *c == b'\t').unwrap_or(hdr.len());
            cur = Some(Record { name: String::from_utf8_lossy(&hdr[..end]).into_owned(), seq: Vec::new() });
        } else if let Some(r) = cur.as_mut() {
            r.seq.extend_from_slice(&line);
        }
    }
    if let Some(r) = cur.take() { recs.push(r); }
    Ok(recs)
}

// ---- config ----------------------------------------------------------------

struct Config {
    name: String,
    pat_fwd: Vec<u8>,
    pat_rev: Vec<u8>,
    /// canonical reference masks (same length as motif) for scoring, if any
    canon: Option<Vec<u8>>,
    /// (left_end, right_start) cut points for arm/spacer decomposition
    arms: (usize, usize),
}

/// Score one hit window (forward genomic substring `win`) on the given strand.
/// Returns the BED line (without trailing newline).
fn format_hit(cfg: &Config, chrom: &str, start: usize, end: usize, plus: bool, win: &[u8]) -> String {
    let l = cfg.pat_fwd.len();
    // canonical-orientation sequence (5'->3' of the motif on its own strand)
    let canon_seq: Vec<u8> = if plus {
        win.iter().map(|b| b.to_ascii_uppercase()).collect()
    } else {
        revcomp(win)
    };
    let strand = if plus { '+' } else { '-' };
    let seq = String::from_utf8_lossy(&canon_seq);

    match &cfg.canon {
        None => format!("{}\t{}\t{}\t{}\t0\t{}\t{}", chrom, start, end, cfg.name, strand, seq),
        Some(cmask) => {
            let (le, rs) = cfg.arms;
            let (mut mis, mut hamm, mut hl, mut hr) = (0u32, 0u32, 0u32, 0u32);
            for k in 0..l {
                if genome_mask(canon_seq[k]) & cmask[k] != 0 {
                    mis += 1;
                } else {
                    hamm += 1;
                    if k < le { hl += 1; } else if k >= rs { hr += 1; }
                }
            }
            let score = (mis as usize * 1000 / l) as u32;
            let larm = &seq[0..le];
            let spacer = &seq[le..rs];
            let rarm = &seq[rs..l];
            format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                chrom, start, end, cfg.name, score, strand,
                mis, hamm, hl, hr, seq, larm, spacer, rarm)
        }
    }
}

fn scan_record(rec: &Record, cfg: &Config, out: &mut String) -> (u64, u64) {
    let l = cfg.pat_fwd.len();
    let s = &rec.seq;
    if s.len() < l { return (0, 0); }
    let gm: Vec<u8> = s.iter().map(|b| genome_mask(*b)).collect();
    let (mut fwd, mut rev) = (0u64, 0u64);
    let last = s.len() - l;
    for i in 0..=last {
        let mut okf = true;
        for k in 0..l { if gm[i + k] & cfg.pat_fwd[k] == 0 { okf = false; break; } }
        if okf { fwd += 1; out.push_str(&format_hit(cfg, &rec.name, i, i + l, true, &s[i..i + l])); out.push('\n'); }
        let mut okr = true;
        for k in 0..l { if gm[i + k] & cfg.pat_rev[k] == 0 { okr = false; break; } }
        if okr { rev += 1; out.push_str(&format_hit(cfg, &rec.name, i, i + l, false, &s[i..i + l])); out.push('\n'); }
    }
    (fwd, rev)
}

// ---- CENP-B flex: two-anchor variable-spacer (edit-distance) scan ----------
//
// Finds left-arm core TTCG ... <variable spacer> ... right-arm core CGGG, each
// arm allowing up to `e` substitutions. Unlike the fixed 17-bp scan this sees
// spacer-LENGTH variants (indels) — the boxes fuzznuc structurally misses.

const FLEX_LCORE: &[u8] = b"TTCG";   // left arm core (canonical), + strand
const FLEX_RCORE: &[u8] = b"CGGG";   // right arm core
const FLEX_LCORE_M: &[u8] = b"CCCG"; // revcomp(CGGG): left anchor of a - strand box on the fwd seq
const FLEX_RCORE_M: &[u8] = b"CGAA"; // revcomp(TTCG): right anchor of a - strand box

struct FlexCfg { name: String, e: u32, smin: usize, smax: usize, max_arm_hd: u32 }

#[inline]
fn hamm4(w: &[u8], exp: &[u8]) -> u32 {
    let mut h = 0;
    for k in 0..4 { if w[k].to_ascii_uppercase() != exp[k] { h += 1; } }
    h
}

/// Format one flex hit (forward box bytes `win`); None if it spans a non-ACGT base
/// or its total arm Hamming exceeds `max_arm_hd`.
fn flex_emit(name: &str, chrom: &str, start: usize, end: usize, plus: bool, win: &[u8], max_arm_hd: u32) -> Option<String> {
    let canon: Vec<u8> = if plus { win.iter().map(|b| b.to_ascii_uppercase()).collect() } else { revcomp(win) };
    let l = canon.len();
    if l < 8 { return None; }
    if canon.iter().any(|b| genome_mask(*b) == 0) { return None; } // no box across a gap/N
    let lh = hamm4(&canon[0..4], FLEX_LCORE);
    let rh = hamm4(&canon[l - 4..l], FLEX_RCORE);
    if lh + rh > max_arm_hd { return None; } // total-arm-error cap
    let spacer_len = l - 8;
    let larm = String::from_utf8_lossy(&canon[0..4]);
    let rarm = String::from_utf8_lossy(&canon[l - 4..l]);
    let spacer = String::from_utf8_lossy(&canon[4..l - 4]);
    let seq = String::from_utf8_lossy(&canon);
    let score = ((8 - lh - rh) as usize * 1000 / 8) as u32;
    let strand = if plus { '+' } else { '-' };
    Some(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        chrom, start, end, name, score, strand, l, larm, lh, spacer_len, spacer, rarm, rh, seq))
}

fn flex_scan_record(rec: &Record, cfg: &FlexCfg, out: &mut String) -> (u64, u64) {
    let s = &rec.seq;
    let n = s.len();
    let min_box = 8 + cfg.smin;
    if n < min_box { return (0, 0); }
    let (mut fwd, mut rev) = (0u64, 0u64);
    let last_i = n - min_box;
    for i in 0..=last_i {
        // + strand: TTCG at i, CGGG at i+4+gap
        if hamm4(&s[i..i + 4], FLEX_LCORE) <= cfg.e {
            for gap in cfg.smin..=cfg.smax {
                let j = i + 4 + gap;
                if j + 4 > n { break; }
                if hamm4(&s[j..j + 4], FLEX_RCORE) <= cfg.e {
                    if let Some(line) = flex_emit(&cfg.name, &rec.name, i, j + 4, true, &s[i..j + 4], cfg.max_arm_hd) {
                        out.push_str(&line); out.push('\n'); fwd += 1;
                    }
                }
            }
        }
        // - strand: CCCG at i, CGAA at i+4+gap (forward representation of a minus box)
        if hamm4(&s[i..i + 4], FLEX_LCORE_M) <= cfg.e {
            for gap in cfg.smin..=cfg.smax {
                let j = i + 4 + gap;
                if j + 4 > n { break; }
                if hamm4(&s[j..j + 4], FLEX_RCORE_M) <= cfg.e {
                    if let Some(line) = flex_emit(&cfg.name, &rec.name, i, j + 4, false, &s[i..j + 4], cfg.max_arm_hd) {
                        out.push_str(&line); out.push('\n'); rev += 1;
                    }
                }
            }
        }
    }
    (fwd, rev)
}

fn run_flex(fasta: &str, cfg: FlexCfg, threads: usize, header: bool) {
    let t0 = Instant::now();
    let recs = Arc::new(read_fasta(fasta).unwrap_or_else(|e| { eprintln!("read error: {e}"); std::process::exit(1); }));
    let n = recs.len();
    let nthreads = threads.max(1).min(n.max(1));
    eprintln!("loaded {} records ({} bp); cenpb-flex arm-tol={} max-arm-hd={} spacer={}-{} on {} threads",
        n, recs.iter().map(|r| r.seq.len()).sum::<usize>(), cfg.e, cfg.max_arm_hd, cfg.smin, cfg.smax, nthreads);
    let stdout = std::io::stdout();
    let mut w = BufWriter::with_capacity(1 << 20, stdout.lock());
    if header {
        writeln!(w, "#chrom\tstart\tend\tname\tscore\tstrand\tbox_len\tlarm\tlarm_hamm\tspacer_len\tspacer\trarm\trarm_hamm\tseq").unwrap();
    }
    let cfg = Arc::new(cfg);
    let mut handles = Vec::new();
    for t in 0..nthreads {
        let recs = Arc::clone(&recs); let cfg = Arc::clone(&cfg);
        handles.push(thread::spawn(move || {
            let mut buf = String::new(); let (mut f, mut r) = (0u64, 0u64); let mut idx = t;
            while idx < recs.len() { let (a, b) = flex_scan_record(&recs[idx], &cfg, &mut buf); f += a; r += b; idx += nthreads; }
            (buf, f, r)
        }));
    }
    let (mut tf, mut tr) = (0u64, 0u64);
    for h in handles { let (buf, f, r) = h.join().unwrap(); w.write_all(buf.as_bytes()).unwrap(); tf += f; tr += r; }
    w.flush().unwrap();
    eprintln!("DONE cenpb-flex total={} fwd={} rev={} elapsed={:.2}s", tf + tr, tf, tr, t0.elapsed().as_secs_f64());
}

// ---- chains: box landmarks -> per-contig chains, optical-map-style match -----
//
// Treats CENP-B boxes as ordered landmarks along each contig (like restriction
// sites in optical mapping). A chain element = (inter-box distance, orientation,
// integrity class). Cross-genome conservation is then a CHAIN match, NOT a
// coordinate liftover — robust to rearrangements. `chain-match` here is a
// PROTOTYPE: bucketed-distance k-gram seeding (generalizes GCP-Centeny model3);
// the production affine optical-map aligner is alphasplitter's align_and_cigar.

use std::collections::HashMap;

fn integ_class(mis: u32) -> char { if mis >= 17 { 'C' } else if mis >= 15 { 'N' } else { 'd' } }

/// Read a --cenpb box BED (path or /dev/stdin) -> contig -> sorted [(start,strand,mis)].
fn load_boxes(path: &str) -> HashMap<String, Vec<(usize, char, u32)>> {
    let f = File::open(path).unwrap_or_else(|e| { eprintln!("open {path}: {e}"); std::process::exit(1); });
    let mut rdr = BufReader::with_capacity(1 << 20, f);
    let mut by: HashMap<String, Vec<(usize, char, u32)>> = HashMap::new();
    let mut line = String::new();
    loop {
        line.clear();
        if rdr.read_line(&mut line).unwrap_or(0) == 0 { break; }
        let t = line.trim_end();
        if t.is_empty() || t.starts_with('#') { continue; }
        let c: Vec<&str> = t.split('\t').collect();
        if c.len() < 7 { continue; }
        let (start, strand, mis) = match (c[1].parse::<usize>(), c[5].chars().next(), c[6].parse::<u32>()) {
            (Ok(s), Some(st), Ok(m)) => (s, st, m), _ => continue,
        };
        by.entry(c[0].to_string()).or_default().push((start, strand, mis));
    }
    for v in by.values_mut() { v.sort_by_key(|x| x.0); }
    by
}

fn run_chains(path: &str) {
    let by = load_boxes(path);
    let stdout = std::io::stdout();
    let mut w = BufWriter::with_capacity(1 << 20, stdout.lock());
    writeln!(w, "#contig\tidx\tstart\tstrand\tinteg\tdist_next").unwrap();
    let mut contigs: Vec<&String> = by.keys().collect();
    contigs.sort();
    for contig in contigs {
        let v = &by[contig];
        for i in 0..v.len() {
            let dist_next: i64 = if i + 1 < v.len() { (v[i + 1].0 - v[i].0) as i64 } else { -1 };
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", contig, i, v[i].0, v[i].1, integ_class(v[i].2), dist_next).unwrap();
        }
    }
    w.flush().unwrap();
    eprintln!("chains: {} contigs", by.len());
}

/// Per-box symbol packing (bucketed dist-to-next, strand, integrity) — the user's
/// 3-part chain element. Distance ALONE is non-specific in tandem arrays (periodic
/// spacing); strand + integrity break the symmetry.
#[inline]
fn elem_sym(v: &[(usize, char, u32)], i: usize, bucket: usize) -> u64 {
    let d = ((v[i + 1].0 - v[i].0) / bucket) as u64;
    let s = if v[i].1 == '-' { 1u64 } else { 0 };
    let ig = match integ_class(v[i].2) { 'C' => 0u64, 'N' => 1, _ => 2 };
    d * 8 + s * 4 + ig
}

/// k-grams over (k-1) consecutive 3-part symbols.
fn elem_kgrams(v: &[(usize, char, u32)], k: usize, bucket: usize) -> Vec<Vec<u64>> {
    if v.len() < k { return Vec::new(); }
    let syms: Vec<u64> = (0..v.len() - 1).map(|i| elem_sym(v, i, bucket)).collect();
    let w = k - 1;
    if syms.len() < w { return Vec::new(); }
    (0..=syms.len() - w).map(|i| syms[i..i + w].to_vec()).collect()
}

fn run_chain_match(a: &str, b: &str, k: usize, bucket: usize, min_shared: usize, max_gram_freq: usize) {
    let ba = load_boxes(a);
    let bb = load_boxes(b);
    // index B's k-grams (occurrence list per gram)
    let mut idx: HashMap<Vec<u64>, Vec<String>> = HashMap::new();
    for (contig, v) in &bb {
        for g in elem_kgrams(v, k, bucket) { idx.entry(g).or_default().push(contig.clone()); }
    }
    // drop non-specific (repetitive) grams that occur in too many places
    let dropped = idx.iter().filter(|(_, c)| c.len() > max_gram_freq).count();
    idx.retain(|_, c| c.len() <= max_gram_freq);
    let stdout = std::io::stdout();
    let mut w = BufWriter::new(stdout.lock());
    writeln!(w, "#contigA\tcontigB\tshared_kgrams").unwrap();
    let mut total = 0u64;
    let mut contigs: Vec<&String> = ba.keys().collect(); contigs.sort();
    for ca in contigs {
        let mut counts: HashMap<&str, usize> = HashMap::new();
        for g in elem_kgrams(&ba[ca], k, bucket) {
            if let Some(cbs) = idx.get(&g) { for cb in cbs { *counts.entry(cb.as_str()).or_default() += 1; } }
        }
        let mut hits: Vec<(&str, usize)> = counts.into_iter().filter(|x| x.1 >= min_shared).collect();
        hits.sort_by(|x, y| y.1.cmp(&x.1));
        for (cb, n) in hits.into_iter().take(5) { writeln!(w, "{}\t{}\t{}", ca, cb, n).unwrap(); total += 1; }
    }
    w.flush().unwrap();
    eprintln!("chain-match: A={} B={} contigs k={} bucket={} max_gram_freq={} (dropped {} repetitive grams) -> {} pairs (>= {} shared)",
        ba.len(), bb.len(), k, bucket, max_gram_freq, dropped, total, min_shared);
}

/// Liftover-free conservation: for each reference chain k-gram, count how many
/// query genomes contain it. Output BED: contig start end cons_count. The
/// reference position lets downstream stratify by repeat/chromatin/position axis.
fn run_chain_conservation(ref_path: &str, query_list: &str, k: usize, bucket: usize) {
    let refb = load_boxes(ref_path);
    // reference gram instances with genomic position (start of the window's first box)
    let mut ref_pos: Vec<(String, usize, usize)> = Vec::new(); // contig, start, end
    let mut ref_gram: Vec<Vec<u64>> = Vec::new();
    let mut contigs: Vec<&String> = refb.keys().collect(); contigs.sort();
    for contig in contigs {
        let v = &refb[contig];
        if v.len() < k { continue; }
        let syms: Vec<u64> = (0..v.len() - 1).map(|i| elem_sym(v, i, bucket)).collect();
        let w = k - 1;
        for i in 0..=syms.len() - w {
            ref_pos.push((contig.clone(), v[i].0, v[i + w].0));
            ref_gram.push(syms[i..i + w].to_vec());
        }
    }
    let mut cons = vec![0u32; ref_gram.len()];
    // load query list
    let ql = std::fs::read_to_string(query_list).unwrap_or_else(|e| { eprintln!("query list {query_list}: {e}"); std::process::exit(1); });
    let queries: Vec<&str> = ql.lines().map(|s| s.trim()).filter(|s| !s.is_empty()).collect();
    let nq = queries.len();
    eprintln!("ref grams={} queries={} k={} bucket={}", ref_gram.len(), nq, k, bucket);
    for (qi, q) in queries.iter().enumerate() {
        let qb = load_boxes(q);
        let mut qset: std::collections::HashSet<Vec<u64>> = std::collections::HashSet::new();
        for (_, v) in &qb { for g in elem_kgrams(v, k, bucket) { qset.insert(g); } }
        for (idx, g) in ref_gram.iter().enumerate() { if qset.contains(g) { cons[idx] += 1; } }
        if (qi + 1) % 50 == 0 { eprintln!("  {}/{} queries", qi + 1, nq); }
    }
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::with_capacity(1 << 20, stdout.lock());
    writeln!(wtr, "#contig\tstart\tend\tcons_count\tn_query").unwrap();
    for (i, (c, s, e)) in ref_pos.iter().enumerate() {
        writeln!(wtr, "{}\t{}\t{}\t{}\t{}", c, s, e, cons[i], nq).unwrap();
    }
    wtr.flush().unwrap();
    eprintln!("DONE chain-conservation: {} ref grams scored over {} queries", ref_gram.len(), nq);
}

// ---- CLI -------------------------------------------------------------------

const CENPB_MOTIF: &str = "NTTCGNNNNANNCGGGN";
const CENPB_CANON: &str = "YTTCGTTGGAARCGGGA";

fn usage() -> ! {
    eprintln!(
"rust_motif_scan — IUPAC motif scanner (both strands), BED output.

USAGE:
  rust_motif_scan <fasta> --motif <IUPAC> [options]
  rust_motif_scan <fasta> --cenpb [options]
  rust_motif_scan <fasta> --cenpb-flex [--spacer MIN-MAX] [--arm-tol E] [options]
  rust_motif_scan chains <boxes.cenpb.bed>            # box BED -> per-contig chains
  rust_motif_scan chain-match <A.bed> <B.bed> [--k 17] [--bucket 5] [--min-shared 1]
                                                     # optical-map-style k-gram chain match (prototype)

OPTIONS:
  --motif <IUPAC>     search pattern (IUPAC; scans both strands)
  --cenpb             preset: --motif {m} --canonical {c} --arms 5,12 --name ECS
  --cenpb-flex        two-anchor variable-spacer scan: TTCG ... <spacer> ... CGGG,
                      each arm allowing E substitutions; finds spacer-LENGTH variants
                      (indels) the fixed scan misses. Columns: chrom start end name
                      score strand box_len larm larm_hamm spacer_len spacer rarm rarm_hamm seq
  --spacer <MIN-MAX>  flex spacer-length range in bp (default 3-15)
  --arm-tol <E>       flex per-arm substitution tolerance (default 1)
  --max-arm-hd <H>    flex cap on TOTAL arm Hamming larm_hamm+rarm_hamm (default 2*E);
                      e.g. --arm-tol 1 --max-arm-hd 1 keeps intact + single-arm-1-mut boxes
  --canonical <IUPAC> reference motif for Hamming scoring (same length as --motif)
  --arms <L,S>        arm cut points for decomposition: left=[0,L) spacer=[L,S) right=[S,len)
  --name <STR>        BED name column (default: the motif, or ECS for --cenpb)
  --threads <N>       worker threads (default: all cores)
  --header            print a #-prefixed column header line
  -h, --help          this help

OUTPUT (BED, tab-separated):
  plain:           chrom start end name 0 strand seq
  with --canonical: chrom start end name score strand mis hamm hammL hammR seq larm spacer rarm
    score = round(mis/len*1000); mis/hamm = IUPAC-aware matches/mismatches to canonical;
    hammL/hammR = mismatches within the left/right arm; seq = canonical orientation (revcomp if -).",
        m = CENPB_MOTIF, c = CENPB_CANON);
    std::process::exit(2);
}

fn main() {
    let argv: Vec<String> = env::args().collect();
    if argv.len() < 2 { usage(); }

    // chain subcommands (input is a --cenpb box BED, not a FASTA)
    match argv[1].as_str() {
        "chains" => { run_chains(&argv.get(2).cloned().unwrap_or_else(|| usage())); return; }
        "chain-match" => {
            let a = argv.get(2).cloned().unwrap_or_else(|| usage());
            let b = argv.get(3).cloned().unwrap_or_else(|| usage());
            let (mut k, mut bucket, mut min_shared, mut max_gram_freq) = (17usize, 5usize, 1usize, 50usize);
            let mut i = 4;
            while i < argv.len() {
                match argv[i].as_str() {
                    "--k" => { i += 1; k = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or(17); }
                    "--bucket" => { i += 1; bucket = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or(5); }
                    "--min-shared" => { i += 1; min_shared = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or(1); }
                    "--max-gram-freq" => { i += 1; max_gram_freq = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or(50); }
                    _ => {}
                }
                i += 1;
            }
            run_chain_match(&a, &b, k.max(2), bucket.max(1), min_shared.max(1), max_gram_freq.max(1));
            return;
        }
        "chain-conservation" => {
            let r = argv.get(2).cloned().unwrap_or_else(|| usage());
            let ql = argv.get(3).cloned().unwrap_or_else(|| usage());
            let (mut k, mut bucket) = (11usize, 5usize);
            let mut i = 4;
            while i < argv.len() {
                match argv[i].as_str() {
                    "--k" => { i += 1; k = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or(11); }
                    "--bucket" => { i += 1; bucket = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or(5); }
                    _ => {}
                }
                i += 1;
            }
            run_chain_conservation(&r, &ql, k.max(2), bucket.max(1));
            return;
        }
        _ => {}
    }

    let mut fasta: Option<String> = None;
    let mut motif: Option<String> = None;
    let mut canonical: Option<String> = None;
    let mut name: Option<String> = None;
    let mut arms: Option<(usize, usize)> = None;
    let mut threads = thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
    let mut header = false;
    let mut cenpb = false;
    let mut flex = false;
    let mut spacer: Option<(usize, usize)> = None;
    let mut arm_tol: Option<u32> = None;
    let mut max_arm_hd: Option<u32> = None;

    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "-h" | "--help" => usage(),
            "--cenpb" => cenpb = true,
            "--cenpb-flex" => flex = true,
            "--header" => header = true,
            "--spacer" => {
                i += 1;
                let v: Vec<usize> = argv.get(i).unwrap_or_else(|| usage()).split('-').filter_map(|s| s.parse().ok()).collect();
                if v.len() != 2 { usage(); }
                spacer = Some((v[0], v[1]));
            }
            "--arm-tol" => { i += 1; arm_tol = Some(argv.get(i).and_then(|s| s.parse().ok()).unwrap_or_else(|| usage())); }
            "--max-arm-hd" => { i += 1; max_arm_hd = Some(argv.get(i).and_then(|s| s.parse().ok()).unwrap_or_else(|| usage())); }
            "--motif" => { i += 1; motif = Some(argv.get(i).cloned().unwrap_or_else(|| usage())); }
            "--canonical" => { i += 1; canonical = Some(argv.get(i).cloned().unwrap_or_else(|| usage())); }
            "--name" => { i += 1; name = Some(argv.get(i).cloned().unwrap_or_else(|| usage())); }
            "--threads" => { i += 1; threads = argv.get(i).and_then(|s| s.parse().ok()).unwrap_or_else(|| usage()); }
            "--arms" => {
                i += 1;
                let v: Vec<usize> = argv.get(i).unwrap_or_else(|| usage()).split(',').filter_map(|s| s.parse().ok()).collect();
                if v.len() != 2 { usage(); }
                arms = Some((v[0], v[1]));
            }
            s if !s.starts_with('-') && fasta.is_none() => fasta = Some(s.to_string()),
            _ => { eprintln!("unknown arg: {}", argv[i]); usage(); }
        }
        i += 1;
    }

    if flex {
        let (smin, smax) = spacer.unwrap_or((3, 15));
        let e = arm_tol.unwrap_or(1);
        let mhd = max_arm_hd.unwrap_or(2 * e); // default: no extra cap beyond per-arm tol
        let nm = name.clone().unwrap_or_else(|| "ECS".to_string());
        let fa = fasta.clone().unwrap_or_else(|| usage());
        run_flex(&fa, FlexCfg { name: nm, e, smin, smax, max_arm_hd: mhd }, threads, header);
        return;
    }

    if cenpb {
        motif.get_or_insert_with(|| CENPB_MOTIF.to_string());
        canonical.get_or_insert_with(|| CENPB_CANON.to_string());
        arms.get_or_insert((5, 12));
        name.get_or_insert_with(|| "ECS".to_string());
    }

    let fasta = fasta.unwrap_or_else(|| usage());
    let motif = motif.unwrap_or_else(|| { eprintln!("error: --motif or --cenpb required"); usage(); });
    let mname = name.unwrap_or_else(|| motif.clone());

    let pat_fwd: Vec<u8> = motif.bytes().map(iupac_mask).collect();
    if pat_fwd.iter().any(|m| *m == 0) { eprintln!("error: --motif has a non-IUPAC char"); std::process::exit(2); }
    let pat_rev: Vec<u8> = pat_fwd.iter().rev().map(|m| comp_mask(*m)).collect();

    let canon = match &canonical {
        None => None,
        Some(c) => {
            if c.len() != motif.len() { eprintln!("error: --canonical len {} != --motif len {}", c.len(), motif.len()); std::process::exit(2); }
            let cm: Vec<u8> = c.bytes().map(iupac_mask).collect();
            if cm.iter().any(|m| *m == 0) { eprintln!("error: --canonical has a non-IUPAC char"); std::process::exit(2); }
            Some(cm)
        }
    };
    let arms = arms.unwrap_or((0, motif.len())); // no decomposition unless set

    let cfg = Arc::new(Config { name: mname, pat_fwd, pat_rev, canon, arms });

    let t0 = Instant::now();
    let recs = Arc::new(read_fasta(&fasta).unwrap_or_else(|e| { eprintln!("read error: {e}"); std::process::exit(1); }));
    let n = recs.len();
    let nthreads = threads.max(1).min(n.max(1));
    eprintln!("loaded {} records ({} bp); motif {} ({}) on {} threads",
        n, recs.iter().map(|r| r.seq.len()).sum::<usize>(), motif,
        if cfg.canon.is_some() { "scored" } else { "plain" }, nthreads);

    let stdout = std::io::stdout();
    let mut w = BufWriter::with_capacity(1 << 20, stdout.lock());
    if header {
        if cfg.canon.is_some() {
            writeln!(w, "#chrom\tstart\tend\tname\tscore\tstrand\tmis\thamm\thammL\thammR\tseq\tlarm\tspacer\trarm").unwrap();
        } else {
            writeln!(w, "#chrom\tstart\tend\tname\tscore\tstrand\tseq").unwrap();
        }
    }

    let mut handles = Vec::new();
    for t in 0..nthreads {
        let recs = Arc::clone(&recs); let cfg = Arc::clone(&cfg);
        handles.push(thread::spawn(move || {
            let mut buf = String::new(); let (mut f, mut r) = (0u64, 0u64);
            let mut idx = t;
            while idx < recs.len() { let (a, b) = scan_record(&recs[idx], &cfg, &mut buf); f += a; r += b; idx += nthreads; }
            (buf, f, r)
        }));
    }
    let (mut tf, mut tr) = (0u64, 0u64);
    for h in handles { let (buf, f, r) = h.join().unwrap(); w.write_all(buf.as_bytes()).unwrap(); tf += f; tr += r; }
    w.flush().unwrap();
    eprintln!("DONE motif={} total={} fwd={} rev={} elapsed={:.2}s", motif, tf + tr, tf, tr, t0.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cenpb_cfg() -> Config {
        let pat_fwd: Vec<u8> = CENPB_MOTIF.bytes().map(iupac_mask).collect();
        let pat_rev: Vec<u8> = pat_fwd.iter().rev().map(|m| comp_mask(*m)).collect();
        let canon: Vec<u8> = CENPB_CANON.bytes().map(iupac_mask).collect();
        Config { name: "ECS".into(), pat_fwd, pat_rev, canon: Some(canon), arms: (5, 12) }
    }

    #[test]
    fn revcomp_canonical() {
        assert_eq!(revcomp(b"CTTCGTTGGAAACGGGA"), b"TCCCGTTTCCAACGAAG".to_vec());
    }

    #[test]
    fn comp_mask_roundtrip() {
        for c in [b'A', b'C', b'G', b'T', b'R', b'Y', b'N'] {
            let m = iupac_mask(c);
            assert_eq!(comp_mask(comp_mask(m)), m);
        }
    }

    #[test]
    fn perfect_plus_hit() {
        let cfg = cenpb_cfg();
        let line = format_hit(&cfg, "chr1", 10, 27, true, b"CTTCGTTGGAAACGGGA");
        let f: Vec<&str> = line.split('\t').collect();
        // chrom start end name score strand mis hamm hammL hammR seq larm spacer rarm
        assert_eq!(f[0], "chr1");
        assert_eq!(f[4], "1000");          // score = 17/17*1000
        assert_eq!(f[5], "+");
        assert_eq!(f[6], "17");            // mis
        assert_eq!(f[7], "0");             // hamm
        assert_eq!(f[8], "0");             // hammL
        assert_eq!(f[9], "0");             // hammR
        assert_eq!(f[10], "CTTCGTTGGAAACGGGA");
        assert_eq!(f[11], "CTTCG");        // larm [0,5)
        assert_eq!(f[12], "TTGGAAA");      // spacer [5,12)
        assert_eq!(f[13], "CGGGA");        // rarm [12,17)
    }

    #[test]
    fn minus_hit_is_revcomp() {
        let cfg = cenpb_cfg();
        // forward window = revcomp of the perfect box -> reported on - strand, seq back to canonical
        let line = format_hit(&cfg, "chr1", 0, 17, false, b"TCCCGTTTCCAACGAAG");
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[5], "-");
        assert_eq!(f[6], "17");            // mis = perfect
        assert_eq!(f[10], "CTTCGTTGGAAACGGGA");
    }

    #[test]
    fn flex_perfect_box() {
        // canonical core-to-core: TTCG + TTGGAAA(7) + CGGG = 15 bp
        let line = flex_emit("ECS", "c", 0, 15, true, b"TTCGTTGGAAACGGG", 8).unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        // chrom start end name score strand box_len larm larm_hamm spacer_len spacer rarm rarm_hamm seq
        assert_eq!(f[4], "1000");   // score (arms perfect)
        assert_eq!(f[6], "15");     // box_len
        assert_eq!(f[7], "TTCG");   // larm
        assert_eq!(f[8], "0");      // larm_hamm
        assert_eq!(f[9], "7");      // spacer_len
        assert_eq!(f[11], "CGGG");  // rarm
        assert_eq!(f[12], "0");     // rarm_hamm
    }

    #[test]
    fn flex_short_spacer_and_arm_mut() {
        // 6-bp spacer (indel) + a left-arm mutation TTCG->TTGG
        let line = flex_emit("ECS", "c", 0, 14, true, b"TTGGTTGAAACGGG", 8).unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[9], "6");      // spacer_len = 14 - 8
        assert_eq!(f[7], "TTGG");   // larm
        assert_eq!(f[8], "1");      // larm_hamm (C->G)
        assert_eq!(f[12], "0");     // rarm_hamm
    }

    #[test]
    fn flex_minus_revcomp() {
        // revcomp(TTCGTTGGAAACGGG) on the fwd seq -> reported on - strand, canonical back
        let fwd = revcomp(b"TTCGTTGGAAACGGG");
        let line = flex_emit("ECS", "c", 0, 15, false, &fwd, 8).unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[5], "-");
        assert_eq!(f[13], "TTCGTTGGAAACGGG"); // canonical orientation
        assert_eq!(f[8], "0");
        assert_eq!(f[12], "0");
    }

    #[test]
    fn flex_max_arm_hd_filters_both_mutated() {
        // TTGG(larm 1mut) + TTGGAAA + CTGG(rarm 1mut) => total arm HD = 2
        assert!(flex_emit("ECS", "c", 0, 15, true, b"TTGGTTGGAAACTGG", 1).is_none()); // capped at 1
        assert!(flex_emit("ECS", "c", 0, 15, true, b"TTGGTTGGAAACTGG", 2).is_some()); // allowed at 2
    }

    #[test]
    fn chain_kgrams_and_class() {
        assert_eq!(integ_class(17), 'C');
        assert_eq!(integ_class(15), 'N');
        assert_eq!(integ_class(10), 'd');
        // boxes at 0,100,205,300; bucket 5; symbol = dist*8 + strand*4 + integ
        // sym0 d20,+,C=160; sym1 d21,+,C=168; sym2 d19,-,N=157
        let v = vec![(0usize, '+', 17u32), (100, '+', 17), (205, '-', 16), (300, '+', 14)];
        let g = elem_kgrams(&v, 3, 5); // k=3 -> 2-symbol grams
        assert_eq!(g.len(), 2);
        assert_eq!(g[0], vec![160, 168]);
        assert_eq!(g[1], vec![168, 157]);
    }

    #[test]
    fn spacer_mismatch_scored() {
        let cfg = cenpb_cfg();
        // pos9 G->C mismatch (inside spacer)
        let line = format_hit(&cfg, "c", 0, 17, true, b"CTTCGTTGCAAACGGGA");
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[6], "16");            // mis
        assert_eq!(f[7], "1");             // hamm
        assert_eq!(f[8], "0");             // hammL (arm intact)
        assert_eq!(f[9], "0");             // hammR (arm intact)
    }
}
