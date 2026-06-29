use clap::Parser;
use realfft::RealFftPlanner;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

const BOTTOM_F: f64 = 20.0;
const TOP_F: f64 = 14000.0;
/// Base EQ point density before CV‑adaptive culling (~1000 max points).
const BASE_RESOLUTION: f64 = 8000.0;

/// Generate EQ correction curves from brown-noise WAV measurements.
///
/// Record brown noise through a reference system to create a target WAV,
/// then record through the headphones with varying mic position to create
/// measurement WAVs.  The tool computes the FFT of each pass, derives
/// per‑frequency mean and variance, and interpolates the correction
/// across frequencies where variance is high.
#[derive(Parser, Debug)]
#[command(name = "autoeq", version, about)]
struct Args {
    /// One or more measurement WAV files (headphones, brown noise, identical length)
    #[arg(short = 'm', long = "measurement", required = true, num_args = 1..)]
    measurements: Vec<PathBuf>,

    /// Target WAV file (reference system recording, same length/sample-rate)
    #[arg(short = 't', long = "target", required = true)]
    target: PathBuf,

    /// Background noise WAV file for spectral subtraction (same length/sample-rate).
    /// Recorded with no signal playing — ambient noise only.
    #[arg(short = 'n', long = "noise")]
    noise: Option<PathBuf>,

    /// Output prefix (generates <prefix>_desktop.csv and <prefix>_mobile.csv)
    #[arg(short = 'o', long = "output", default_value = "autoeq")]
    output: String,

    /// High-frequency rolloff: FREQ,DB_PER_OCTAVE (e.g. "4000,-6"). Repeatable.
    #[arg(long = "high-rolloff", value_name = "FREQ,DB_OCTAVE", allow_hyphen_values = true)]
    high_rolloff: Vec<String>,

    /// Low-frequency rolloff: FREQ,DB_PER_OCTAVE (e.g. "100,-12"). Repeatable.
    #[arg(long = "low-rolloff", value_name = "FREQ,DB_OCTAVE", allow_hyphen_values = true)]
    low_rolloff: Vec<String>,

    /// Boost in dB to apply to the final EQ [default: 0.0]
    #[arg(short = 'b', long = "boost", default_value = "0.0")]
    boost: f64,

    /// Coefficient of variation target.  Determines EQ point density:
    /// spacing scales as CV / target‑cv.  Lower target → denser points.
    /// Consistent measurements (low CV) naturally get more points;
    /// noisy ones get fewer.  [default: 0.12]
    #[arg(long = "target-cv", default_value = "0.12")]
    target_cv: f64,

    /// Cutoff frequency for the harmonic bass enhancer's target compensation.
    /// Below cutoff the target is rebuilt as a mirror of the above‑cutoff
    /// response, scaled so the Chebyshev harmonics (h₂·x², h₃·x³) reproduce
    /// the desired level at 2f / 3f.  Sub‑fc/3 is set to −40 dB (silent).
    /// Above cutoff the target is left untouched.
    #[arg(long = "bass-enhancer-cutoff")]
    bass_enhancer_cutoff: Option<f64>,
}

fn main() {
    let args = Args::parse();

    // ── 1. Read & validate all WAV files ──────────────────────────────────
    let wavs: Vec<WavData> = args
        .measurements
        .iter()
        .map(|p| read_wav(p))
        .collect();
    let target_wav = read_wav(&args.target);

    let sample_rate = wavs[0].sample_rate;
    for w in &wavs[1..] {
        assert_eq!(
            w.sample_rate, sample_rate,
            "WAV sample rate mismatch — all files must use the same rate"
        );
    }
    assert_eq!(
        target_wav.sample_rate, sample_rate,
        "Target WAV sample rate must match measurement WAVs"
    );



    // ── 2. Welch FFT each measurement WAV → per-bin mean + variance     ──
    let all_stats: Vec<Vec<BinStats>> = wavs
        .iter()
        .map(|w| welch_stats(&w.samples, sample_rate))
        .collect();

    // ── 3. Pool per-bin stats across measurement WAVs                   ──
    let bin_count = all_stats[0].len();
    for s in &all_stats[1..] {
        assert_eq!(s.len(), bin_count, "WAV freq ranges differ — sample rates must match");
    }

    let pooled: Vec<BinStats> = (0..bin_count)
        .map(|i| {
            let total_n: u64 = all_stats.iter().map(|s| s[i].count).sum();
            let total_sum: f64 = all_stats
                .iter()
                .map(|s| s[i].mean * s[i].count as f64)
                .sum();
            // Pool CVs by weighted average
            let pooled_cv: f64 = all_stats
                .iter()
                .map(|s| s[i].cv * s[i].count as f64)
                .sum::<f64>()
                / total_n as f64;
            BinStats {
                freq: all_stats[0][i].freq,
                count: total_n,
                mean: total_sum / total_n as f64,
                cv: pooled_cv,
            }
        })
        .collect();

    let total_windows = pooled[0].count;
    eprintln!(
        "{} windows across {} measurement WAV(s)",
        total_windows, all_stats.len(),
    );

    // ── 4. Build full‑resolution mean Graph + CV Graph                  ──
    let measurement_full = Graph {
        points: pooled.iter().map(|s| Point { x: s.freq, y: s.mean }).collect(),
    };
    let mut cv_full = Graph {
        points: pooled.iter().map(|s| Point { x: s.freq, y: s.cv }).collect(),
    };

    // ── 5. Noise: spectral subtraction + combine CV into point gen      ──
    let noise_full_res: Option<Graph> = args.noise.as_ref().map(|noise_path| {
        let noise_wav = read_wav(noise_path);
        assert_eq!(noise_wav.sample_rate, sample_rate);
        let noise_stats = welch_stats(&noise_wav.samples, sample_rate);

        // Combine noise CV into measurement CV — unstable background
        // degrades correction confidence at those frequencies.
        for i in 0..cv_full.points.len().min(noise_stats.len()) {
            cv_full.points[i].y = cv_full.points[i].y.max(noise_stats[i].cv);
        }

        Graph {
            points: noise_stats.iter().map(|s| Point { x: s.freq, y: s.mean }).collect(),
        }
    });

    // ── 6. Generate adaptive EQ points (CV now includes noise CV)       ──
    let eq_pts = eq_points_adaptive(&cv_full, args.target_cv);
    eprintln!("{} EQ points (adaptive)", eq_pts.len());

    // ── 7. Resample measurement + spectral subtraction                   ──
    let mut measurement_mean = measurement_full.resample(eq_pts.clone());

    if let Some(ref noise_full) = noise_full_res {
        let noise_mean = noise_full.resample(eq_pts.clone());
        for i in 0..measurement_mean.points.len() {
            let m = measurement_mean.points[i].y;
            let n = noise_mean.points[i].y;
            measurement_mean.points[i].y = (m - n).max(m * 0.01);
        }
    }

    // ── 8. Target: Welch → full‑res Graph → resample → shape          ──
    let target_stats = welch_stats(&target_wav.samples, sample_rate);
    let target_full = Graph {
        points: target_stats.iter().map(|s| Point { x: s.freq, y: s.mean }).collect(),
    };
    let mut target = target_full.resample(eq_pts.clone());

    // Manual rolloffs
    for spec in &args.high_rolloff {
        let (freq, db) = parse_rolloff(spec);
        target = target.rolloff(freq, db, RolloffMode::High);
    }
    for spec in &args.low_rolloff {
        let (freq, db) = parse_rolloff(spec);
        target = target.rolloff(freq, db, RolloffMode::Low);
    }
    // ── 9. Correction = target / measurement_mean                        ──
    let mut comp = Graph {
        points: eq_pts
            .iter()
            .map(|&freq| Point {
                x: freq,
                y: target.point(freq) / measurement_mean.point(freq),
            })
            .collect(),
    };

    // Bass enhancer psychoacoustic compensation
    if let Some(cutoff) = args.bass_enhancer_cutoff {
        comp = comp.bass_enhancer_preprocess(cutoff);
    }

    // ── 10. Write outputs                                                 ──
    let desktop_path = format!("{}_desktop.csv", args.output);
    let mut outfile = File::create(&desktop_path).unwrap();
    outfile
        .write_all(comp.to_jamesdsp_eq(args.boost).as_bytes())
        .unwrap();
    eprintln!("Wrote {}", desktop_path);

    let mobile_path = format!("{}_mobile.csv", args.output);
    let mut outfile_mobile = File::create(&mobile_path).unwrap();
    outfile_mobile
        .write_all(comp.to_jamesdsp_eq_android(args.boost).as_bytes())
        .unwrap();
    eprintln!("Wrote {}", mobile_path);
}

// ── WAV I/O ──────────────────────────────────────────────────────────────────

struct WavData {
    samples: Vec<f64>,
    sample_rate: u32,
}

/// Read a mono or stereo WAV, returning mono f64 samples (stereo = channel avg).
/// Handles 8/16/24/32-bit integer and 32-bit float formats.
fn read_wav(path: &Path) -> WavData {
    let mut reader = hound::WavReader::open(path).expect("Could not open WAV file");
    let spec = reader.spec();
    let sample_rate = spec.sample_rate;
    let channels = spec.channels;
    let bits = spec.bits_per_sample;
    let fmt = spec.sample_format;

    let samples: Vec<f64> = match channels {
        1 => read_samples(&mut reader, bits, fmt),
        2 => {
            let all = read_samples(&mut reader, bits, fmt);
            let left: Vec<f64> = all.iter().copied().step_by(2).collect();
            let right: Vec<f64> = all.iter().copied().skip(1).step_by(2).collect();
            left.iter()
                .zip(right.iter())
                .map(|(l, r)| (l + r) / 2.0)
                .collect()
        }
        n => panic!("Unsupported channel count: {}", n),
    };

    WavData { samples, sample_rate }
}

fn read_samples(
    reader: &mut hound::WavReader<std::io::BufReader<std::fs::File>>,
    bits: u16,
    fmt: hound::SampleFormat,
) -> Vec<f64> {
    match (bits, fmt) {
        (16, hound::SampleFormat::Int) => reader
            .samples::<i16>()
            .map(|s| s.expect("WAV read error") as f64 / 32768.0)
            .collect(),
        (24, hound::SampleFormat::Int) => reader
            .samples::<i32>()
            .map(|s| s.expect("WAV read error") as f64 / 8388608.0)
            .collect(),
        (32, hound::SampleFormat::Int) => reader
            .samples::<i32>()
            .map(|s| s.expect("WAV read error") as f64 / 2147483648.0)
            .collect(),
        (8, hound::SampleFormat::Int) => reader
            .samples::<i8>()
            .map(|s| s.expect("WAV read error") as f64 / 128.0)
            .collect(),
        (32, hound::SampleFormat::Float) => reader
            .samples::<f32>()
            .map(|s| s.expect("WAV read error") as f64)
            .collect(),
        _ => panic!("Unsupported sample format: {:?} {}bit", fmt, bits),
    }
}

// ── FFT ──────────────────────────────────────────────────────────────────────

/// FFT size for Welch overlapped-window averaging.  Must be a power of two.
const WELCH_FFT_SIZE: usize = 16384;
/// Overlap fraction (0.5 = 50% overlap).
const WELCH_OVERLAP: f64 = 0.5;

/// Per‑bin statistics from Welch windows.
struct BinStats {
    freq: f64,
    count: u64,    // number of windows that contributed
    mean: f64,     // raw average magnitude (for correction curve)
    cv: f64,       // level‑normalized coefficient of variation (for reliability)
}

/// Run Welch's method and return per‑bin mean (raw, for correction) and
/// per‑bin stats normalized by each window's own mean level (for CV).
///
/// Normalization removes common‑mode gain variation between windows
/// (mic distance drift) so CV only reflects frequency‑specific variance.
fn welch_stats(samples: &[f64], sample_rate: u32) -> Vec<BinStats> {
    let fft_size = WELCH_FFT_SIZE;
    let step = (fft_size as f64 * (1.0 - WELCH_OVERLAP)) as usize;
    let bin_count = fft_size / 2 + 1;

    // Hann window
    let window: Vec<f64> = (0..fft_size)
        .map(|n| {
            0.5 * (1.0 - (2.0 * std::f64::consts::PI * n as f64 / (fft_size - 1) as f64).cos())
        })
        .collect();

    let mut planner = RealFftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(fft_size);

    // Raw magnitudes → mean spectrum for correction
    let mut raw_sum = vec![0.0f64; bin_count];
    // Level‑normalized magnitudes → CV of spectral shape
    let mut norm_sum = vec![0.0f64; bin_count];
    let mut norm_sum_sq = vec![0.0f64; bin_count];
    let mut window_count = 0u64;

    let mut offset = 0;
    while offset + fft_size <= samples.len() {
        let mut input = vec![0.0f64; fft_size];
        for i in 0..fft_size {
            input[i] = samples[offset + i] * window[i];
        }
        let mut output = fft.make_output_vec();
        fft.process(&mut input, &mut output).expect("FFT failed");

        // Gather magnitudes for this window
        let mags: Vec<f64> = output.iter()
            .map(|c| c.norm() / fft_size as f64)
            .collect();

        // Window mean (drops out overall gain variation)
        let mean_mag = mags.iter().sum::<f64>() / bin_count as f64;

        for i in 0..bin_count {
            raw_sum[i] += mags[i];
            if mean_mag > 0.0 {
                let nm = mags[i] / mean_mag;
                norm_sum[i] += nm;
                norm_sum_sq[i] += nm * nm;
            } else {
                norm_sum[i] += 1.0;
                norm_sum_sq[i] += 1.0;
            }
        }
        window_count += 1;
        offset += step;
    }

    let n = window_count as f64;
    let freqs: Vec<f64> = (0..bin_count)
        .map(|i| i as f64 * sample_rate as f64 / fft_size as f64)
        .collect();

    (0..bin_count)
        .filter(|&i| freqs[i] >= BOTTOM_F && freqs[i] <= TOP_F)
        .map(|i| {
            let mean = raw_sum[i] / n;
            let norm_mean = norm_sum[i] / n;
            let norm_var = (norm_sum_sq[i] / n) - (norm_mean * norm_mean);
            let cv = if norm_mean > 0.0 {
                norm_var.sqrt() / norm_mean
            } else {
                f64::INFINITY
            };
            BinStats {
                freq: freqs[i],
                count: window_count,
                mean: mean.max(0.0),
                cv,
            }
        })
        .collect()
}

// ── EQ point generation ──────────────────────────────────────────────────────

fn parse_rolloff(spec: &str) -> (f64, f64) {
    let mut parts = spec.split(',');
    let freq: f64 = parts
        .next()
        .expect("Rolloff spec missing frequency")
        .parse()
        .expect("Invalid rolloff frequency");
    let db: f64 = parts
        .next()
        .expect("Rolloff spec missing dB/octave")
        .parse()
        .expect("Invalid rolloff dB/octave");
    (freq, db)
}

// ── Enums ────────────────────────────────────────────────────────────────────

enum RolloffMode {
    High,
    Low,
}

// ── dB conversions ───────────────────────────────────────────────────────────

fn db_to_ratio(db: f64) -> f64 {
    10.0_f64.powf(db / 20.0)
}

fn ratio_to_db(ratio: f64) -> f64 {
    ratio.log10() * 20.0
}

/// Generate EQ points with density determined by per‑frequency CV.
///
/// Base density follows equal‑loudness contours.  Spacing is scaled by
/// `cv / target_cv` — clean measurements (low CV) get dense points;
/// noisy ones get sparse points with no artificial floor.
fn eq_points_adaptive(cv_graph: &Graph, target_cv: f64) -> Vec<f64> {
    let base_density = |freq: f64| -> f64 {
        let equal_loudness = Graph {
            points: vec![
                Point { x: 20.0, y: 109.0 },
                Point { x: 80.0, y: 82.0 },
                Point { x: 400.0, y: 62.0 },
                Point { x: 1000.0, y: 60.0 },
                Point { x: 1500.0, y: 64.0 },
                Point { x: 2500.0, y: 57.0 },
                Point { x: 4000.0, y: 57.0 },
                Point { x: 8500.0, y: 73.0 },
                Point { x: 15000.0, y: 72.0 },
                Point { x: 19000.0, y: 68.0 },
                Point { x: 30000.0, y: 130.0 },
            ],
        };
        BASE_RESOLUTION / equal_loudness.point(freq)
    };

    let mut freq = BOTTOM_F;
    let mut out = vec![];
    while freq < TOP_F {
        out.push(freq);
        let base_step = freq / base_density(freq) * 2.0 / std::f64::consts::E;
        let cv = cv_graph.point(freq).max(0.001);
        freq += base_step * (cv / target_cv);
    }
    out.push(TOP_F);
    out
}

// ── Butterworth filter response ──────────────────────────────────────────────

/// 2nd-order Butterworth low‑pass magnitude response at frequency f with cutoff fc.
fn butterworth_lp(f: f64, fc: f64) -> f64 {
    let w = f / fc;
    1.0 / (1.0 + w * w * w * w).sqrt()
}

/// 2nd-order Butterworth high‑pass magnitude response at frequency f with cutoff fc.
fn butterworth_hp(f: f64, fc: f64) -> f64 {
    let w = f / fc;
    (w * w) / (1.0 + w * w * w * w).sqrt()
}

// ── Graph ────────────────────────────────────────────────────────────────────

#[derive(Clone)]
struct Graph {
    points: Vec<Point>,
}

impl Graph {
    fn to_jamesdsp_eq(&self, boost: f64) -> String {
        let mut out = String::new();
        let dbs = self.map(ratio_to_db);
        let max = dbs.max();
        for point in dbs.points {
            out.push_str(
                format!("{:.3}\t{:.3}\n", point.x, (point.y - max + boost).min(0.0)).as_str(),
            )
        }
        out
    }

    fn to_jamesdsp_eq_android(&self, boost: f64) -> String {
        let mut out = "GraphicEQ: ".to_string();
        let dbs = self.map(ratio_to_db);
        let max = dbs.max();
        for point in dbs.points {
            out.push_str(
                format!("{:.3} {:.3};", point.x, (point.y - max + boost).min(0.0)).as_str(),
            )
        }
        out
    }

    fn rolloff(&self, freq: f64, rolloff_per_octave: f64, mode: RolloffMode) -> Self {
        let atten = db_to_ratio(rolloff_per_octave);
        Graph {
            points: self
                .points
                .iter()
                .map(|point| {
                    let ratio = match mode {
                        RolloffMode::High => point.x / freq,
                        RolloffMode::Low => freq / point.x,
                    };
                    let octave = ratio.max(1.0).log2();
                    Point {
                        x: point.x,
                        y: point.y * atten.powf(octave),
                    }
                })
                .collect(),
        }
    }

    /// Preprocess a compensation curve for harmonic bass enhancer input.
    ///
    /// The enhancer high‑passes input at fc, low‑passes the bass into two
    /// bands (T₃ LP at fc/2, T₂ LP at fc), and waveshapes each with
    /// Chebyshev polynomials to generate harmonics *above* cutoff.  The
    /// EQ's job below cutoff is to provide the correct *input level* for
    /// the waveshaper — nothing passes through directly.
    ///
    /// Both the input LP and output HP (2nd‑order Butterworth) are
    /// compensated so each harmonic lands exactly on the target level.
    /// Drive includes 1/H_LP(f) and 1/ⁿ√H_HP(harmonic) to cancel the
    /// filter response for each path.
    ///
    ///   T₃ [fc/3, fc/2] → 3f:  ∛(comp(3f)/H_HP(3f)) / H_LP_T3(f)
    ///   T₂ [fc/2, fc]   → 2f:  √(comp(2f)/H_HP(2f)) / H_LP_T2(f)
    ///
    /// A crossfade across [fc/2, fc] blends from the compensated T₂‑drive
    /// to comp(fc) — the direct compensation *at the cutoff* — so the
    /// curve meets the untouched region with no step.
    fn bass_enhancer_preprocess(&self, fc: f64) -> Self {
        let f_lo = fc / 3.0;
        let f_mid = fc / 2.0;
        let comp_fc = self.point(fc);

        Graph {
            points: self
                .points
                .iter()
                .map(|point| {
                    let f = point.x;
                    let y = if f < f_lo {
                        // Below fc/3: no Chebyshev harmonic passes the output HP.
                        0.0
                    } else if f < f_mid {
                        // T₃ band [fc/3, fc/2]: drives 3rd harmonic at 3f.
                        let hp = butterworth_hp(3.0 * f, fc);
                        let lp = butterworth_lp(f, fc / 2.0);
                        (self.point(3.0 * f) / hp).cbrt() / lp
                    } else if f < fc {
                        // T₂ band [fc/2, fc]: compensated drive, crossfaded
                        // to comp(fc) for a seamless boundary at cutoff.
                        let hp = butterworth_hp(2.0 * f, fc);
                        let lp = butterworth_lp(f, fc);
                        let drive = (self.point(2.0 * f) / hp).sqrt() / lp;
                        let t = (f - f_mid) / (fc - f_mid);
                        drive * (1.0 - t) + comp_fc * t
                    } else {
                        // Above cutoff — untouched
                        point.y
                    };
                    Point { x: f, y }
                })
                .collect(),
        }
    }

    fn point(&self, x: f64) -> f64 {
        if let Some(lower) = self
            .points
            .windows(2)
            .enumerate()
            .find(|(_, fs)| fs[0].x <= x && fs[1].x >= x)
        {
            let lower = lower.0;
            let p1 = self.points[lower];
            let p2 = self.points[lower + 1];
            let slope = (p2.y - p1.y) / (p2.x - p1.x);
            p1.y + (x - p1.x) * slope
        } else if x < self.points[0].x {
            self.points[0].y
        } else {
            self.points.last().unwrap().y
        }
    }

    fn map(&self, tform: impl Fn(f64) -> f64) -> Graph {
        Graph {
            points: self
                .points
                .iter()
                .map(|point| Point {
                    x: point.x,
                    y: tform(point.y),
                })
                .collect(),
        }
    }

    /// Resample to target x positions using 20-sample RMS averaging per point.
    fn resample(&self, target: Vec<f64>) -> Graph {
        Graph {
            points: target
                .iter()
                .enumerate()
                .map(|(i, point)| {
                    let left = if i > 0 { target[i - 1] } else { *point };
                    let right = target.get(i + 1).unwrap_or(point);
                    let x1 = (point + left) / 2.0;
                    let x2 = (point + right) / 2.0;

                    let mut sum = 0_f64;
                    for i in 0..20 {
                        let x = x1 + (x2 - x1) * i as f64 / 20.0;
                        sum += self.point(x).powf(2.0);
                    }
                    Point {
                        x: *point,
                        y: (sum / 20.0).sqrt(),
                    }
                })
                .collect(),
        }
    }

    fn max(&self) -> f64 {
        self.points
            .iter()
            .map(|point| point.y)
            .reduce(|a, b| a.max(b))
            .unwrap()
    }
}

#[derive(Clone, Copy, Debug)]
struct Point {
    x: f64,
    y: f64,
}
