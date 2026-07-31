# VeroMass Aligner

Standalone Python/Tkinter desktop tool for LC-MS chromatographic peak
detection and cross-sample retention-time alignment/correspondence. Reads a
folder of LC-MS runs and produces a single aligned feature table across all
samples: chromatographic peak picking -> correspondence grouping -> retention
time correction -> re-grouping -> Excel export.

Pure Python — no R runtime dependency, no rpy2.

## Accepted inputs (v1.3.0)

`.mzML`, `.mzXML`, `.raw` (Thermo), `.mgf`, and `.zip` archives of any of these
(extracted recursively — `expand_inputs`, handles nested zips such as a dataset
zip whose entries are per-sample zips each containing one mzML). The folder
scanner picks up all of these. The **Convert tab** runs `expand_inputs` +
`convert_raw`; the **Align tab** (with `skip_conversion=True`) skips both and
goes straight to peak detection on already-converted files.

- **MGF** spectra are ALREADY peak-picked (one precursor per BEGIN/END IONS
  block). There is no MS1 chromatographic axis, so the MGF path skips ROI
  building / peak detection entirely — each spectrum becomes one peak directly
  (`peaks_from_mgf`). This is by far the fastest path.
- MGF fragment lists are parsed LAZILY (`read_mgf_spectra` stores the raw
  `frag_text`; `_get_pairs` / `_parse_frag_text` parse floats only for the
  ~one-spectrum-per-feature actually chosen for MS2 output). A typical MGF has
  ~240 fragments/spectrum, so eager parsing would convert millions of floats up
  front — lazy parsing cut real 2-file MGF runs from ~9.6s to ~3s.

## Two-tab architecture (Convert | Align)

The app has a **two-tab notebook** — Convert and Align are now separate jobs,
each with its own thread, queue, and log. See `[[two-tab-architecture-verified]]`.

### Convert tab (Tab 0)
- **Pipeline**: `scan_input_files` → `expand_inputs` (zip extraction) → `convert_raw` (RAW→mzML, concurrent threads) + copy passthrough MS files into output folder
- **Output**: self-contained folder of mzML-ready files
- **Controls**: Start / Stop / Open Output Folder (no Pause/Reset)
- **Auto-link**: on completion, sets the Align tab's input folder to the conversion output

### Align tab (Tab 1)
- **Pipeline**: `scan_input_files` → `run_alignment(skip_conversion=True)` → detection → grouping → RT correction → re-grouping → MS2 → export
- **Input**: mzML/mzXML/mgf files (already converted — no expansion/RAW conversion)
- **Controls**: Start / Pause / Stop / Reset / Open Output Folder
- **Parameters**: all 8 fields seeded from `PeakDetectionParams()` / `GroupingParams()` defaults

### Testing
Run `python VeroMass_Aligner.py`. Test Convert first (with `.raw` or `.zip` files),
then switch to Align tab (input auto-fills). Both logs are independent per tab.

## Entry point

`VeroMass_Aligner.py` — single file, same convention as sibling tools
(`mgf-extractor`, `MoleculeID_Processor`, `phyto-crossmatcher`): Tkinter dark
UI, `_ensure_deps()` runtime pip bootstrap, background-thread processing with
`queue.Queue` progress polling, launched via `Start_VeroMass_Aligner_Windows.bat`.

## Pipeline (in file order)

1. **File I/O** — `read_ms_scans_mzml` / `read_ms_scans_mzxml`: return
   `(ms1_scans, ms2_scans)`. MS1 scans keep full centroid m/z+intensity
   arrays (needed for chromatographic peak detection, not just TIC/base-peak).
   MS2 scans keep precursor m/z + fragment pairs for later fragment lookup.
   `.raw` files are converted to mzML first via `convert_raw()`.
2. **Raw conversion** (`convert_raw` / `_run_with_heartbeat`) —
   **ThermoRawFileParser only** (msConvert is NOT probed, since target
   machines are not assumed to have ProteoWizard installed and the probe
   itself could block). TRFP is auto-downloaded (pinned v1.4.5) and cached
   under `~/VeroMass_Aligner/ThermoRawFileParser/`. Conversion runs via
   `Popen` with a heartbeat log line every 10s and a hard 30-min timeout, so a
   slow/large file never looks like a frozen UI and never hangs forever.
3. **Peak detection** (`detect_chrom_peaks` / `_build_rois` /
   `_detect_peaks_in_roi`) — builds m/z-tolerance ROIs across consecutive
   scans, then does Gaussian-smoothed multi-scale local-maxima picking along
   each ROI's RT-intensity trace, producing per-file chromatographic peaks
   (m/z/RT apex + bounds, height, area, S/N).
4. **Correspondence** (`group_peaks_density` since v1.8.0; legacy `group_peaks`
   retained) — coarse m/z gap-partitioning, then intensity-weighted KDE density
   clustering over RT within each slice. Produces cross-sample features with
   size/m.z/RT/width/shape-distance stats.
5. **RT correction** (`compute_rt_correction` / `apply_rt_correction`) —
   anchor features present in most samples get a per-sample loess RT correction
   curve (with linear extrapolation past the anchor range), applied to all peaks.
6. **Re-grouping** — `group_peaks_density` re-run on RT-corrected peaks for the
   final aligned feature table.
7. **MS2 attachment** (`attach_ms2`) — for each feature, MS2 is looked up
   **only in its representative sample** (the sample with the highest peak
   height for that feature), matched by precursor m/z within grouping
   tolerance and RT within the feature window. This is deliberate: it removes
   redundant identical MS2 for the same m/z across samples and keeps the
   output file small. Fragments are stored as a base-peak-normalised top-20
   string.
8. **Export** (`write_feature_table`) — `.xlsx` with three sheets:
   - `Features` — one row per aligned feature: `feature, Size, Charge, Mass,
     m.z, RT, Base.Peak, m.z.Width, RT.Height, m.z.Min, m.z.Max, RT.Min,
     RT.Max, Shape.Distance, MS2.Sample, MS2.RT, MS2.Fragments`.
   - `Peaks` — one row per detected peak per sample per feature: the full list
     of detected m/z, RT, and intensities (Height + Area) across every
     chromatogram used in the alignment.
   - `Intensities` — wide per-sample intensity matrix (one column per sample).

## Naming rule

Do NOT use "xcms" / "XCMS" anywhere in code, filenames, UI text, comments, or
docs. This tool is VeroMass-branded. (An earlier draft was named "XCMS
Aligner"; fully renamed 2026-07-20.)

## Detection recall — the v1.4.0 correctness fix (READ THIS FIRST)

A scientist reported the tool "does not capture all" features. He was right.
Measured on the real `L60_ddms_neg` run (2,279 MS1 scans) against a deliberately
conservative ground truth (apex >= 5x its own local baseline with a clean rise
and fall over 3 points each side — peaks any scientist accepts on sight):

| version | peaks | recall |
|---|---|---|
| v1.3.0 (as shipped to the scientist) | 9,382 | **24.3%** |
| + ROI temporal contiguity + baseline noise | 13,710 | 66.3% |
| **v1.4.0 (also drop prominence, min width 1s)** | **43,189** | **98.0%** |

Three real defects, all in peak detection:

1. **ROI temporal contiguity was lost in the v1.2.0 vectorized rewrite.** The
   fast rewrite grouped every point sharing an m/z across the WHOLE run into
   one ROI. Measured: only 13 of 4,000 ROIs were a single contiguous elution;
   the median held 5 separate segments, the worst 248. A compound's real 10 s
   peak was merged with all scattered background at that m/z.
   Fix: split each m/z bucket wherever the scan index jumps > `max_gap_scans`.
2. **The noise floor was self-referential.** `max(noise, median(inten)*0.5)`
   derived "noise" from the trace's own median, so a well-sampled, genuinely
   eluting compound defined its own noise floor and was rejected at S/N ~2.3 —
   real abundant peaks discarded for being consistently present.
   Fix: baseline = p25 of the elution window, sigma = MAD of the lower half
   (x1.4826), S/N = (height - baseline) / sigma.
3. **The `prominence` gate was the single largest killer (29% of all real
   peaks).** Topographic prominence is measured against neighbouring maxima,
   not the baseline, so genuine peaks were rejected purely for eluting near
   another peak — and it was redundant with the S/N test. Removed. Width floor
   also lowered 5s -> 1s (real median peak width here is ~3.8 s; the surviving
   width distribution had been clipped exactly at the old cutoff, proving it
   was cutting into real peaks).

**Do not reintroduce a prominence filter, a median-based noise floor, or a
whole-run ROI without re-running the recall benchmark.** Over-detection is far
less harmful than under-detection here — downstream grouping (`min_frac_samples`)
and cloud matching filter further, but a feature never detected is unrecoverable.

## v1.5.0 — intensity correctness + the real slowness cause

Reported on two ~200 MB Thermo `.raw` files (`312_R1`, `573_R1`, 19-min runs,
~2,150 MS1 scans, ~1,450 centroids/scan): "took way longer, fewer peaks and
lower intensities." Two distinct real defects, both fixed:

1. **Intensity was under-reported — the tool returned the SMOOTHED value as the
   peak height** (`height = smoothed[pk]`). Smoothing is a detection aid for
   *locating* a peak; reporting it as the intensity understates every peak.
   Measured on real data: median **17% low**, worst decile **29% low**, total
   summed ion intensity only **85% of true**. Now reports the true measured
   apex (max raw intensity inside the peak bounds) — this is exactly xcms's
   `maxo` (maximum observed intensity) convention, versus `into` for the
   integrated area. Side effect: peak counts rose ~40% (78k -> 112k per file)
   because the S/N test now uses the true apex instead of a suppressed value.
2. **`attach_ms2` was O(features x MS2 scans)** — a linear scan of every MS2
   scan for every feature (~50k x ~6k = ~300M Python comparisons). It did not
   finish in 5+ minutes and scaled directly with the much larger feature counts
   the v1.4.0 recall fixes produce. Now indexes each sample's MS2 by precursor
   m/z once and binary-searches: **minutes -> 0.5s**.

Verified end-to-end on those two raw files (warm conversion cache):
**140 s total**, 111,889 + 119,125 peaks, **49,168 features**, 231,014 peak
rows, detection recall **96.5% / 97.7%**. Grouping drops zero peaks.

Remaining known cost: the Excel export is ~66 s of that 140 s (231k-row Peaks
sheet, ~2.6M cells via openpyxl). Optimizing it means either the `xlsxwriter`
engine (not currently installed) or emitting the large Peaks table as CSV.

## v1.5.1 — the GUI silently overrode every tuning fix (READ BEFORE TUNING)

The v1.4.0/v1.5.0 detection fixes were verified by calling
`PeakDetectionParams()` directly — but **every real run goes through the GUI**,
whose parameter fields were hardcoded literals that OVERRIDE the dataclass
defaults on `_start()`. A leftover `peak min (s) = 5` there meant users kept
getting the old over-aggressive width gate no matter what the code default said.

Measured on `312_R1.raw`, same version, same ROIs (60,724), only that field
differing:

| GUI `peak min (s)` | peaks | summed intensity |
|---|---|---|
| 5 (stale literal) | 5,636 | 3.70e9 |
| 1 (dataclass default) | **111,889** | **4.80e10** |

20x the peaks, 13x the intensity. The user's delivered workbook had 12,293 peak
rows where a correct run produces 231,014.

Fix: the UI fields are now **seeded from `PeakDetectionParams()` /
`GroupingParams()`** instead of repeating literals, so they cannot drift again.
**Never reintroduce hardcoded parameter literals in `_build_ui`** — and when
benchmarking a detection change, run it through the GUI path (or at minimum
construct params the way `_start()` does), not just the dataclass defaults,
or you will validate a configuration no user ever runs.

## v1.6.0 — large per-peak tables go to CSV

Confirmed from a real user run: the Excel write was 44 s of a ~2m17s run, all
of it serializing the 231k-row Peaks sheet. Benchmarked on that exact shape:

| writer | 231k rows |
|---|---|
| openpyxl | 72.9 s |
| xlsxwriter | 50.0 s (not worth a new dependency) |
| **CSV** | **2.1 s** |

So when the per-peak table exceeds `PEAKS_XLSX_MAX_ROWS` (100,000) it is written
as a companion `aligned_peaks.csv` and the Peaks worksheet is omitted; the
workbook keeps Features + Intensities. Measured export **44 s -> 24 s**.
`write_feature_table` now returns `(xlsx_path, peaks_csv_path_or_None)` and the
run log names the CSV. Small studies are unaffected — they keep the single-file
workbook with all three sheets.

## v1.7.0 — chain-clustering bug corrupted EVERY exported RT (partially fixed)

Reported by the user: aligned this tool's output against real xcms results on
`312_R1.raw`/`573_R1.raw` and found systematic RT mismatches vs xcms, not just
missing peaks. Root-caused precisely (not guessed): **peak detection was
completely correct** — raw per-file `detect_chrom_peaks` output matched xcms's
apex RT and intensity almost exactly for every example checked (e.g.
`713.229`: xcms RT 6.633, raw detected apex 6.625, height `200,674,448` —
matching the exported `Base.Peak` to the exact digit). The corruption was
downstream, in `group_peaks`'s correspondence clustering:

- **The bug**: both the m/z-bucketing and RT-sub-clustering steps used
  **single-linkage/chain clustering** — checking each point only against its
  immediate neighbor (`cur[-1]`), with no bound on total cluster width. A
  chain of closely-spaced points can link a cluster far wider than the
  intended tolerance. Measured: a feature reported `Size=7` from a
  **2-sample** alignment (should be at most 2) — that contamination pulled
  the cluster's height-weighted consensus RT away from the true apex, which
  then became a bad anchor for `compute_rt_correction`, which then applied a
  **spurious +0.24 to +0.35 min shift to every peak in both samples** — not
  just the contaminated ones. This is why comparing against xcms looked like
  "missing peaks" system-wide: even correctly-detected, correctly-isolated
  peaks came out with the wrong RT because the correction curve itself was
  poisoned.
- **The fix**: bound both clustering steps to a **fixed window from the
  seed/first member**, never chain-expanding (`hi_mz` no longer grows past
  the seed's own tolerance; the RT sub-cluster bound is measured from
  `cur[0]`, not `cur[-1]`). Verified: systemic RT error dropped **~5x** (0.25-
  0.35 min -> 0.06-0.07 min) on the same real files, and previously
  wrongly-merged decoy peaks are now correctly kept separate.
- **NOT fully fixed — a second, distinct issue remains**: bucket assignment
  is still **order-dependent**. Demonstrated directly: a real 344.133 Da
  compound's two same-sample peaks (312_R1 `mz=344.13429`, 573_R1
  `mz=344.13205`, corrected RTs only ~0.126 min apart — well inside the 0.25
  min bandwidth) merge correctly in an isolated 3-peak repro, but do NOT
  merge in the real ~230k-peak run. Mechanism: an unrelated peak elsewhere in
  the dense dataset can form its own overlapping m/z window first and claim
  one of the trio before it gets its proper turn in the single left-to-right
  pass — a real cross-sample match can still be split apart in dense regions.
  **This needs a real rewrite** (pick the best match per point, not
  first-available-in-scan-order — e.g. mutual-nearest-neighbor or a proper
  density-based method), not a small patch. Deliberately deferred as a
  separate, scoped follow-up rather than rushed into this pass.
- **Do not re-introduce chain/single-linkage clustering** (checking only
  `cur[-1]`/growing a window indefinitely) anywhere in `group_peaks` — that
  is the exact class of bug this fix addressed, and it silently corrupts
  results that look fine in isolation but break under real peak density.

## v1.8.0 — density-based correspondence + real loess RT correction

The order-dependent "bucket stealing" left open at v1.7.0 was reworked properly.
`group_peaks` (seed-window scheme) is retained but no longer the default;
correspondence now runs through `group_peaks_density` / `_density_cluster_bin`,
the same family of algorithm xcms's `PeakDensityParam` uses:

- Coarse m/z partitioning via gap-splitting (safe here — fine separation is
  delegated to the RT-density step, not the gap threshold). **Do NOT gap-split
  the RT step** — a mid-development attempt to gap-split BOTH axes is
  mathematically single-linkage/chain clustering and measured markedly worse
  (54.3% match, `size>2` reappeared). Reverted.
- Within each m/z slice: build an intensity-weighted Gaussian KDE over RT
  (bandwidth `rt_bw_sec`), take the global density maximum, capture points
  within `rt_capture_mult * bw`, enforce one-peak-per-sample. No seed, no
  greedy first-mover — assignment falls out of the aggregate density shape, so
  processing order provably cannot change the result.
- `compute_rt_correction` replaced its global degree-2 polyfit (which *clamped*
  RT at the anchor-range edges — worst exactly where correction is needed most)
  with a real loess (`_loess_curve`, tricube kernel, span=0.2, matching xcms's
  `smooth="loess"`) plus *linear extrapolation* past the anchor range.

Verified against the real reference report (60,512 features, two real
positive-mode files): match 73.1%->73.5%, RT bias +0.065->+0.0001 min, RT std
0.145->0.134, illegal `size>2` in a 2-sample run 32,063->0. Cost: correspondence
~10s->~115s per pass (accepted — accuracy was the priority).

## v1.8.1 — silent-drop correctness fix + the accuracy ceiling was mapped

Independently re-verified all v1.8.0 numbers to the digit (isolated path AND the
real `run_alignment` end-to-end path, incl. the exported workbook: 73.5%,
size2=8,923, size>2=0). Then two things were found and one fixed:

1. **Silent-drop bug in `_density_cluster_bin` (FIXED).** The peeling loop
   removed ALL captured peaks from the pool but only emitted the per-sample
   *winners*; same-sample non-winners were silently discarded — real detected
   chromatographic peaks that never reached the output, contradicting the
   function's own docstring (step 4: "return the rest to the pool"). Measured
   discard: **2,912 peaks (1.3%) at `rt_capture_mult=0.1`, 86,770 (37.6%) at
   0.5.** This radius-dependent bleed is exactly why narrower `rt_capture_mult`
   always won in tuning — it was not a true optimum, just less dropped signal.
   Fix: remove only the kept winners; leave non-winners in the pool for a later
   density maximum. Re-verified: match 73.51%->73.54%, features
   219,128->221,989 (recovered signal), `size>2` still 0. The residual reason
   0.1 still edges 0.5 even after the fix (73.54 vs 72.08) is *legitimate* —
   a too-wide radius folds genuinely distinct nearby compounds into one density
   maximum — so **`rt_capture_mult=0.1` stays the default.**

2. **Correspondence is at its detection ceiling — the match rate is NOT a
   correspondence problem.** Of the 16,027 unmatched reference features, **99.2%
   are detection misses** (no raw peak within 30ppm/0.5min); only 0.8% (123) are
   lost in grouping. Raw-peak ceiling @30ppm/0.5 = 73.7%, grouping achieves
   73.5%. Further correspondence work will not move the match rate.

3. **The "73.5%" metric is heavily coincidence-inflated.** A decoy control
   (permute reference RT, keep m/z, re-match) shows a **31.9% coincidence floor**
   at the tightest tolerance — because the tool emits 219k features vs the
   reference's 60k (3.6x). True agreement above chance is ~42%. Widening the
   match tolerance recovers *no* real signal (true signal falls). The lever is
   **detection over-production** (mostly narrow spikes the reference's
   `peakwidth=c(9.4,32)` rejects), which is a recall/precision tradeoff the
   scientist must weigh — do NOT unilaterally tighten detection (the tool exists
   because a scientist wanted MORE recall).

**Reference-file column semantics (previously unresolved):** `rtmed` is in
MINUTES; `mean` = exact arithmetic mean of the two per-sample columns; `npeaks`
= count of *actually detected* peaks (the per-sample intensity matrix is
`fillChromPeaks`-backfilled, so a feature can be npeaks==1 yet have both columns
populated); no systematic m/z bias (signed ppm err +0.00). The `maxint` vs
per-sample-column relationship is *not* a constant transform (consistent with
apex-height `maxo` vs integrated-area `into`, but a bright row implies an
impossible 11,000-scan peak width) — still needs the scientist's export-script
definition to fully close. The paper's `RT < 1020s` cutoff does NOT explain the
misses (rt<17min matches 73.4%, rt>=17 matches 77.5%) — that lead is closed.

## Known limitations / not yet implemented

- RT correction uses a real loess (`_loess_curve`, tricube-weighted local
  linear regression, span=0.2) with linear extrapolation past the anchor range,
  as of v1.8.0 (was previously a global degree-2 polyfit that clamped at the
  edges).
- No isotope/adduct-aware charge detection — `Charge` is hardcoded to 1 and
  `Mass` assumes `[M+H]+`.
- **Peak-detection performance (v1.2.0 rewrite — the important one)**: the
  ROI builder is now FULLY VECTORIZED. Instead of the incremental centWave-style
  loop (carry a set of "active" ROIs, match every m/z point in every scan
  against that set — O(scans × active_ROIs), which ballooned to ~70k active
  ROIs and ~90-104s per 2,300-scan file), `_build_rois` flattens all MS1
  points, drops sub-prefilter noise in one mask, sorts once by m/z, and splits
  into ROIs on m/z gaps (`np.diff`/`np.where`). `_detect_peaks_in_roi` runs
  `find_peaks` on the native scan trace (no fixed-grid upsampling).
  Measured: real 446-scan file 5.6s -> 0.6s; simulated 2,300-scan file ~90s ->
  ~1.4s (~65x). Memory is O(total_points) flat arrays — NOT the dense
  O(n_bins × n_scans) matrix that OOM-crashed an earlier implementation
  (`reference_lcms_test_datasets` memory).
  - Peak detection is left SERIAL: it's CPU-bound Python and the GIL makes a
    ThreadPool *slower* than serial (measured 5.2s threaded vs 3.9s serial for
    4 mzML files). True parallelism there would need processes; not worth it
    now that each file is ~1s.
  - RAW->mzML conversions ARE pre-run concurrently (ThreadPool) — they're
    external ThermoRawFileParser subprocesses (~30s each) that don't hold the
    GIL, so N conversions collapse to ~one conversion's wall-time.
  - **Always benchmark peak detection on a real full-length run, not just
    synthetic or small files** — the original bottleneck was invisible on
    small/synthetic input.
- No automated test suite — verify by running the `.bat` launcher against a
  real folder of mzML files and inspecting `aligned_features.xlsx`.

## Real-data parser gotcha (fixed, keep in mind)

`read_ms_scans_mzml` originally used `elem_a or elem_b` to pick the `<binary>`
element. ElementTree Elements with no child elements are falsy under `bool()`,
so a real `<binary>` (text-only) node evaluated False and the `or` silently
discarded it — every real mzML returned 0 scans, with no error. Fixed with an
explicit `is None` check. The same anti-pattern still exists in
`MoleculeID_Processor.py` (a spawned task was filed to fix it there).
**Always test with real mzML, not just synthetic data** — this bug class is
invisible on synthetic input.

## Test data

Real mzML test files: `IBD_small.zip` (10-sample set) — see project memory
`reference_lcms_test_datasets`. That zip is itself a zip-of-zips (each sample
is a `.zip` containing one `.mzML`), which the recursive `expand_inputs` now
handles directly. Real MGF test files: the `L49/L50_ddms_*` set (also in that
memory) — ~6k spectra each, good for exercising the fast MGF path.

Note: MGF has no continuous RT axis, so it can't be chromatographically
peak-detected — but the tool handles it by treating each pre-picked spectrum
as a peak directly (it does NOT run centWave on MGF).
