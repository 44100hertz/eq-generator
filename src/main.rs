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

    /// Output prefix (generates <prefix>_desktop.csv and <prefix>_mobile.csv)
    #[arg(short = 'o', long = "output", default_value = "autoeq")]
    output: String,

    /// High-frequency rolloff: FREQ,DB_PER_OCTAVE (e.g. "4000,-6"). Repeatable.
    #[arg(long = "high-rolloff", value_name = "FREQ,DB_OCTAVE", allow_hyphen_values = true)]
    high_rolloff: Vec<String>,

    /// Low-frequency rolloff: FREQ,DB_PER_OCTAVE (e.g. "100,-12"). Repeatable.
    #[arg(long = "low-rolloff", value_name = "FREQ,DB_OCTAVE", allow_hyphen_values = true)]
    low_rolloff: Vec<String>,

    /// Crispy peak amplitude (e.g. 0.8)
    #[arg(long = "crispy-peak")]
    crispy_peak: Option<f64>,

    /// Boost in dB to apply to the final EQ [default: 0.0]
    #[arg(short = 'b', long = "boost", default_value = "0.0")]
    boost: f64,

    /// Coefficient of variation target.  Determines EQ point density:
    /// spacing scales as CV / target‑cv.  Lower target → denser points.
    /// Consistent measurements (low CV) naturally get more points;
    /// noisy ones get fewer.  [default: 0.12]
    #[arg(long = "target-cv", default_value = "0.12")]
    target_cv: f64,

    /// Harman target curve intensity (0 = flat, 1 = full Harman 2018 over-ear).
    /// Applied above 122 Hz only; bass shelf is controlled separately.
    /// [default: 0.0]
    #[arg(long = "harman-intensity", default_value = "0.0")]
    harman_intensity: f64,

    /// Bass shelf level in dB below 122 Hz.  Positive = boost, negative = cut.
    /// Independent from Harman intensity.  [default: 0.0]
    #[arg(long = "bass-level", default_value = "0.0", allow_hyphen_values = true)]
    bass_level: f64,
}

fn main() {
    let args = Args::parse();

    // ── 1. Read & validate all WAV files ──────────────────────────────────
    let mut wavs: Vec<WavData> = args
        .measurements
        .iter()
        .map(|p| read_wav(p))
        .collect();
    let target_wav = read_wav(&args.target);

    let sample_rate = wavs[0].sample_rate;
    let len = wavs[0].samples.len();
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

    // Trim all to shortest length if they differ slightly
    let min_len = wavs.iter()
        .map(|w| w.samples.len())
        .chain(std::iter::once(target_wav.samples.len()))
        .min()
        .unwrap();
    for w in &mut wavs {
        w.samples.truncate(min_len);
    }
    let target_wav = WavData {
        samples: target_wav.samples[..min_len].to_vec(),
        sample_rate: target_wav.sample_rate,
    };
    if min_len != len {
        eprintln!("Trimmed all WAVs to {} samples (shortest file)", min_len);
    }

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
    let cv_full = Graph {
        points: pooled.iter().map(|s| Point { x: s.freq, y: s.cv }).collect(),
    };

    // ── 5. Generate adaptive EQ points                                    ──
    let eq_pts = eq_points_adaptive(&cv_full, args.target_cv);
    eprintln!("{} EQ points (adaptive)", eq_pts.len());

    // ── 6. Resample to adaptive EQ points                                 ──
    let measurement_mean = measurement_full.resample(eq_pts.clone());

    // ── 7. Target: Welch → full‑res Graph → resample → shape          ──
    let target_stats = welch_stats(&target_wav.samples, sample_rate);
    let target_full = Graph {
        points: target_stats.iter().map(|s| Point { x: s.freq, y: s.mean }).collect(),
    };
    let mut target = target_full.resample(eq_pts.clone());

    // Apply Harman 2018 over‑ear curve (above 122 Hz only)
    if args.harman_intensity > 0.0 {
        target = target.apply_harman(args.harman_intensity);
    }
    // Bass shelf below 122 Hz (independent from Harman)
    if args.bass_level.abs() > 0.01 {
        target = target.bass_shelf(122.0, args.bass_level);
    }
    // Manual rolloffs
    for spec in &args.high_rolloff {
        let (freq, db) = parse_rolloff(spec);
        target = target.rolloff(freq, db, RolloffMode::High);
    }
    for spec in &args.low_rolloff {
        let (freq, db) = parse_rolloff(spec);
        target = target.rolloff(freq, db, RolloffMode::Low);
    }
    // Crispy peak (dB‑domain, adds a sinc‑shaped boost at ~11 kHz)
    if let Some(db) = args.crispy_peak {
        target = target.crispy_peak(db);
    }

    // ── 8. Correction = target / measurement_mean                        ──
    let comp = Graph {
        points: eq_pts
            .iter()
            .map(|&freq| Point {
                x: freq,
                y: target.point(freq) / measurement_mean.point(freq),
            })
            .collect(),
    };

    // ── 9. Write outputs                                                  ──
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

// ── Interpolation for unreliable EQ points ───────────────────────────────────

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

// ── Harman target curve ──────────────────────────────────────────────────────

/// Harman 2018 over‑ear target, dB relative to 0 at 122 Hz.
/// Data: jaakkopasanen/AutoEq Harman over-ear 2018.csv.
/// Below 122 Hz is flat — bass shelf is applied separately via `--bass-level`.
fn harman_curve() -> Graph {
    Graph {
        points: vec![
            Point { x: 20.0, y: 0.00 },
            Point { x: 122.0, y: 0.00 },
            Point { x: 122.3, y: 0.02 },
            Point { x: 142.0, y: -0.69 },
            Point { x: 164.9, y: -1.27 },
            Point { x: 191.4, y: -1.82 },
            Point { x: 222.2, y: -2.09 },
            Point { x: 258.0, y: -1.97 },
            Point { x: 299.5, y: -1.67 },
            Point { x: 347.8, y: -1.34 },
            Point { x: 403.7, y: -1.12 },
            Point { x: 468.7, y: -0.94 },
            Point { x: 544.2, y: -0.69 },
            Point { x: 631.8, y: -0.45 },
            Point { x: 733.4, y: -0.26 },
            Point { x: 851.5, y: -0.15 },
            Point { x: 988.6, y: -0.02 },
            Point { x: 1147.7, y: 0.40 },
            Point { x: 1332.5, y: 1.24 },
            Point { x: 1546.9, y: 2.37 },
            Point { x: 1795.9, y: 3.92 },
            Point { x: 2085.0, y: 5.71 },
            Point { x: 2420.7, y: 7.16 },
            Point { x: 2810.3, y: 8.09 },
            Point { x: 3262.7, y: 8.54 },
            Point { x: 3787.9, y: 8.40 },
            Point { x: 4397.6, y: 7.58 },
            Point { x: 5105.5, y: 6.28 },
            Point { x: 5927.3, y: 5.21 },
            Point { x: 6881.4, y: 3.96 },
            Point { x: 7989.1, y: 2.44 },
            Point { x: 9275.1, y: 0.55 },
            Point { x: 10768.1, y: -2.00 },
            Point { x: 12501.4, y: -5.02 },
            Point { x: 14513.7, y: -7.50 },
            Point { x: 16850.0, y: -10.65 },
            Point { x: 19562.3, y: -20.79 },
            Point { x: 20000.0, y: -23.00 },
        ],
    }
}

// ── Graph ────────────────────────────────────────────────────────────────────

#[derive(Clone)]
struct Graph {
    points: Vec<Point>,
}

const PI: f64 = std::f64::consts::PI;

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

    /// Crispy peak: adds a sinc‑shaped boost at ~11 kHz.
    /// Amplitude is now in dB (e.g. 0.8 = +0.8 dB at the peak).
    fn crispy_peak(&self, db: f64) -> Self {
        // Convert to dB, add the peak, convert back
        let dbs = self.map(ratio_to_db);
        Graph {
            points: dbs
                .points
                .iter()
                .map(|point| {
                    let fkhz = point.x / 1000.0;
                    let sinc = |p: f64| {
                        if p < -PI || p > PI {
                            0.0
                        } else {
                            p.sin() / p
                        }
                    };
                    let boost = sinc(fkhz * 2.0 - 22.0) * db;
                    Point {
                        x: point.x,
                        y: db_to_ratio(point.y + boost),
                    }
                })
                .collect(),
        }
    }

    /// Apply the Harman 2018 over‑ear target curve (above 122 Hz only).
    /// `intensity` blends between flat (0.0) and full Harman (1.0).
    fn apply_harman(&self, intensity: f64) -> Self {
        let harman = harman_curve();
        Graph {
            points: self
                .points
                .iter()
                .map(|point| {
                    let h_db = if point.x >= 122.0 {
                        harman.point(point.x) * intensity
                    } else {
                        0.0
                    };
                    Point {
                        x: point.x,
                        y: point.y * db_to_ratio(h_db),
                    }
                })
                .collect(),
        }
    }

    /// Constant gain shelf below `freq` Hz.  Smooth transition over one
    /// octave above `freq` (raised cosine blend).
    fn bass_shelf(&self, freq: f64, db: f64) -> Self {
        let f_hi = freq * 2.0; // one octave transition
        let gain = db_to_ratio(db);
        Graph {
            points: self
                .points
                .iter()
                .map(|point| {
                    let factor = if point.x <= freq {
                        gain
                    } else if point.x >= f_hi {
                        1.0
                    } else {
                        // raised cosine blend
                        let t = (point.x - freq) / (f_hi - freq);
                        let blend = 0.5 + 0.5 * (PI * t).cos();
                        gain * blend + 1.0 * (1.0 - blend)
                    };
                    Point {
                        x: point.x,
                        y: point.y * factor,
                    }
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
