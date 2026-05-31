# rust_motif_scan

Fast, **zero-dependency** IUPAC degenerate-motif scanner. Scans a FASTA on **both strands**, writes **BED**, and — given a canonical reference — scores every hit by **Hamming distance** to that reference, decomposed into **left-arm / spacer / right-arm** regions.

Built for the genomics case where you scan a degenerate motif across many large genomes and want both the occurrences *and* a per-hit integrity score. Ships a dedicated **`--cenpb`** mode for the centromeric **CENP-B box**, whose two conserved "arms" flanking a variable spacer motivated the arm-aware output.

- **No dependencies.** `cargo build --release`, or even `rustc -O src/main.rs`.
- **Multithreaded** across FASTA records (`std::thread`).
- **Fast:** scans the 3.1 Gbp T2T-CHM13v2.0 genome for a 17-mer degenerate motif, both strands, in **~8 s** on 192 threads (vs EMBOSS `fuzznuc` ~42 s, `seqkit locate` ~6.5 min — identical hit counts).

## Install

```bash
git clone https://github.com/ad3002/rust_motif_scan.git
cd rust_motif_scan
cargo build --release           # binary at target/release/rust_motif_scan
# or, no cargo:
rustc -O src/main.rs -o rust_motif_scan
```

## Usage

```bash
# General: any IUPAC motif, both strands -> BED6 + sequence
rust_motif_scan genome.fa --motif GGWGGW > hits.bed

# Score hits against a canonical reference (same length) with arm decomposition
rust_motif_scan genome.fa --motif NTTCGNNNNANNCGGGN \
    --canonical YTTCGTTGGAARCGGGA --arms 5,12 --name ECS > ecs.bed

# CENP-B box preset (the two lines above, built in)
rust_motif_scan genome.fa --cenpb --header > ecs.bed

# Edit-distance / two-anchor mode: TTCG ... <variable spacer> ... CGGG, each arm
# allowing E substitutions — finds spacer-LENGTH variants the fixed scan misses
rust_motif_scan genome.fa --cenpb-flex --spacer 3-15 --arm-tol 1 --header > boxes_flex.bed
```

### `--cenpb-flex` (two-anchor, variable spacer)

The fixed `--cenpb` scan is a 17-bp window — substitutions only. `--cenpb-flex` instead anchors on the two conserved arm cores `TTCG` … `CGGG` (each allowing `--arm-tol` substitutions) separated by a spacer of length `--spacer MIN-MAX`, so it also reports **spacer-length (indel) variants**. Columns:

```
chrom start end name score strand box_len larm larm_hamm spacer_len spacer rarm rarm_hamm seq
```

`box_len = 8 + spacer_len`; `larm_hamm`/`rarm_hamm` = substitutions in the left/right arm core (≤ arm-tol); `score = round((8−larm_hamm−rarm_hamm)/8 · 1000)`; `seq`/arms/spacer in canonical orientation. Answers e.g. "how many intact-arm boxes have a non-canonical spacer length?" (`larm_hamm=0 & rarm_hamm=0 & spacer_len≠7`) and the arm-mutation spectrum.

`--arm-tol 1` alone is permissive (both arms may carry a substitution → mostly composition noise). Use **`--max-arm-hd`** to cap the *total* arm Hamming: `--arm-tol 1 --max-arm-hd 1` keeps only **intact + single-arm-single-mutation** boxes (the clean, anchored set). Note: in low-complexity / repeat regions one left anchor can pair with several right anchors — aggregate (or take the nearest pair) downstream.

### Options

| flag | meaning |
|---|---|
| `--motif <IUPAC>` | search pattern (IUPAC alphabet; scanned on both strands) |
| `--cenpb` | preset = `--motif NTTCGNNNNANNCGGGN --canonical YTTCGTTGGAARCGGGA --arms 5,12 --name ECS` |
| `--canonical <IUPAC>` | reference motif (same length as `--motif`) enabling Hamming scoring |
| `--arms <L,S>` | arm cut points: left = `[0,L)`, spacer = `[L,S)`, right = `[S,len)` |
| `--name <STR>` | BED name column (default: the motif string; `ECS` under `--cenpb`) |
| `--threads <N>` | worker threads (default: all cores) |
| `--header` | emit a `#`-prefixed column header line |
| `-h, --help` | help |

## Output (BED, tab-separated)

**Plain** (no `--canonical`):

```
chrom  start  end  name  0  strand  seq
```

**Scored** (`--canonical` / `--cenpb`):

```
chrom  start  end  name  score  strand  mis  hamm  hammL  hammR  seq  larm  spacer  rarm
```

- `start` 0-based, `end` exclusive (half-open BED).
- `strand` `+` / `-`; **`seq` is the canonical-orientation motif** — the reverse complement is taken for `-` hits so every row reads 5'→3' in the motif's own frame.
- `mis` / `hamm` = IUPAC-aware matches / mismatches of `seq` vs `--canonical` over all positions (`hamm = len − mis`).
- `hammL` / `hammR` = mismatches within the left arm `[0,L)` / right arm `[S,len)`.
- `score` = `round(mis / len * 1000)` (BED6-valid 0–1000; higher = closer to canonical).
- `larm` / `spacer` / `rarm` = the canonical-orientation substrings.

## The CENP-B box mode

The CENP-B box is a 17-bp centromeric motif bound by CENP-B; **9 nucleotides are essential for binding** and sit in two conserved blocks — `TTCG` (positions 2–5) and `CGGG` (positions 13–16), each carrying a CpG — flanking a more variable **spacer**. `--cenpb` searches the degenerate docking pattern `NTTCGNNNNANNCGGGN` (both strands) and scores each hit against the canonical CENP-B box **`YTTCGTTGGAARCGGGA`** (Y=C/T, R=A/G).

Because the search fixes the 9 docking nucleotides, the **arms come out conserved by construction** (low `hammL`/`hammR`) and the deviation concentrates in the spacer + flanks — which is exactly the biology the arm decomposition makes visible. Default arms `5,12` give left = `[0,5)`, spacer = `[5,12)`, right = `[12,17)`.

> Note: the raw genome-wide count of `NTTCGNNNNANNCGGGN` (~166 k in CHM13) is dominated by centromeric α-satellite. "Ectocentromeric sites" (chromosome-arm occurrences) are a downstream subset obtained by excluding centromeric regions — apply your own region filter to `ecs.bed`.

## Semantics & correctness

- A genome base matches only if it is exactly `A/C/G/T`; any other byte (`N`, gap, IUPAC ambiguity in the genome) matches nothing, so **no hit is ever called across an assembly gap**.
- Pattern uses the full IUPAC alphabet (`N`=ACGT, `R`=AG, `Y`=CT, …). Hit counts are identical to EMBOSS `fuzznuc -complement` and `seqkit locate -d` (validated on CHM13: 166 059 hits for `NTTCGNNNNANNCGGGN`).

## License

MIT © Aleksey Komissarov
