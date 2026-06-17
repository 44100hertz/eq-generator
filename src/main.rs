use std::fs::{self, File};
use std::io::Write;

const BOTTOM_F: f64 = 20.0;
const TOP_F: f64 = 14000.0;
const RESOLUTION: f64 = 800.0;

fn main() {
    let measurements = [
        // "ALDI/1.txt",
        // "ALDI/2.txt",
        "kenrad/L2.txt",
        "kenrad/R2.txt",
        // "corolla/2.txt",
        // "corolla/3.txt",
        // "corolla/4.txt",
        // "corolla/5.txt",
    ]
    .iter()
    .map(|file| Graph::from_audacity_spectogram(file).resample(eq_points(), AverageMode::Rms))
    .collect::<Vec<_>>();

    let floor = Graph::from_audacity_spectogram("jamo/noise.txt")
       .resample(eq_points(), AverageMode::Mean);
    let target = Graph::from_audacity_spectogram("kenrad/target-brown.txt")
        .resample(eq_points(), AverageMode::Rms)
        // .rolloff(200.0, -9.0, RolloffMode::Low)
        // .rolloff(5000.0, -6.0, RolloffMode::High)
        // .rolloff(1000.0, 5.0, RolloffMode::High);

        // Half-harmon + anti-subbass + crispy peak
        .rolloff(4000.0, -6.0, RolloffMode::High)
        .rolloff(1000.0, 2.5, RolloffMode::High)
        .rolloff(300.0, 1.0, RolloffMode::High)
        .rolloff(100.0, -12.0, RolloffMode::Low)
        .rolloff(180.0, 3.0, RolloffMode::Low)
        .crispyPeak(0.8)
        // Harmon curve + Anti-subbass
//         .rolloff(4000.0, -9.0, RolloffMode::High)
//         .rolloff(1000.0, 5.0, RolloffMode::High)
         //.rolloff(140.0, 5.0, RolloffMode::Low)
         //.rolloff(100.0, -12.0, RolloffMode::Low)
         //.rolloff(180.0, -2.0, RolloffMode::Low)
         .sum(&floor);
        // 
//        .rolloff(4000.0, -9.0, RolloffMode::High)
//        .rolloff(1000.0, 5.0, RolloffMode::High)
//        .rolloff(140.0, 5.0, RolloffMode::Low)
//        .rolloff(100.0, -12.0, RolloffMode::Low);
        // .rolloff(200.0, 0.5, RolloffMode::High);
        
        // .rolloff(200.0, -3.0, RolloffMode::Low);
        // .sum(&floor);

    let rms = rms(measurements);

    let comp = Graph {
        points: eq_points()
            .into_iter()
            .map(|freq| Point {
                x: freq,
                y: target.point(freq) / rms.point(freq),
            })
            .collect(),
    };

    let mut outfile = File::create("autoeq_desktop.csv").unwrap();
    outfile
        .write_all(comp.to_jamesdsp_eq(0.0).as_bytes())
        .unwrap();

    let mut outfile_mobile = File::create("autoeq_mobile.csv").unwrap();
    outfile_mobile
        .write_all(comp.to_jamesdsp_eq_android(0.0).as_bytes())
        .unwrap();
}

enum AverageMode {
    Mean,
    Rms,
}

enum RolloffMode {
    High,
    Low,
}

// Assumes that graphs all have the same indices...
// And that all y values are negative, no less than 120
fn rms(graphs: Vec<Graph>) -> Graph {
    Graph {
        points: (0..graphs[0].points.len())
            .map(|i| {
                let sum_squares = graphs
                    .iter()
                    .fold(0.0, |acc, graph| acc + graph.points[i].y.powf(2.0));
                let rms_sum = (sum_squares / graphs.len() as f64).sqrt();
                Point {
                    x: graphs[0].points[i].x,
                    y: rms_sum,
                }
            })
            .collect(),
    }
}

// fn deviation(graphs: Vec<Graph>) -> Graph {
//     Graph {
//         points: (0..graphs[0].points.len())
//             .map(|i| {
//                 let mean = graphs
//                     .iter()
//                     .fold(0.0, |acc, graph| acc + graph.points[i].y);

//                 let deviation = graphs
//                     .iter()
//                     .fold(0.0, |acc, graph| acc + (mean - graph.points[i].y).powf(2.0) / graphs.len() as f64)
//                     .sqrt();

//                 let rms_sum = (sum_squares / graphs.len() as f64).sqrt();
//                 Point {
//                     x: graphs[0].points[i].x,
//                     y: rms_sum,
//                 }
//             })
//             .collect(),
//     }
// }

fn db_to_ratio(db: f64) -> f64 {
    10.0_f64.powf(db / 20.0)
}
fn db_to_linear(db: f64) -> f64 {
    10.0_f64.powf((db + 60.0) / 20.0)
}
fn ratio_to_db(ratio: f64) -> f64 {
    ratio.log10() * 20.0
}

fn eq_points() -> Vec<f64> {
    // how many points per octave?
    fn density(freq: f64) -> f64 {
        // based on a rough tracing of equal loudness countour
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
                Point { x: 15000.0, y: 72.0, },
                Point { x: 19000.0, y: 68.0, },
                Point { x: 30000.0, y: 130.0, },
            ],
        };
        let dbs = equal_loudness.point(freq);
        RESOLUTION / dbs
    }
    let mut freq = BOTTOM_F;
    let mut out = vec![];
    while freq < TOP_F {
        out.push(freq);
        freq += freq / density(freq) * 2.0 / std::f64::consts::E;
    }
    out.push(TOP_F);
    out
}

struct Graph {
    points: Vec<Point>,
}

const PI: f64 = std::f64::consts::PI;

impl Graph {
    fn from_audacity_spectogram(path: &str) -> Self {
        let file = fs::read_to_string(path).expect("Could not read file");
        let mut lines = file.split("\n");
        // skip over headers
        lines.next().expect("Empty file!");
        Graph {
            points: lines
                .filter(|line| !line.is_empty())
                .map(|line| {
                    let mut fields = line.split("\t");
                    Point {
                        x: fields.next().unwrap().parse().unwrap(),
                        y: db_to_linear(fields.next().unwrap().parse().unwrap()),
                    }
                })
                .collect(),
        }
    }

    fn to_jamesdsp_eq(&self, boost: f64) -> String {
        let mut out = String::new();
        let dbs = self.map(ratio_to_db);
        let max = dbs.max();
        for point in dbs.points {
            out.push_str(format!("{:.3}\t{:.3}\n", point.x, (point.y - max + boost).min(0.0)).as_str())
        }
        out
    }

    fn to_jamesdsp_eq_android(&self, boost: f64) -> String {
        let mut out = "GraphicEQ: ".to_string();
        let dbs = self.map(ratio_to_db);
        let max = dbs.max();
        for point in dbs.points {
            out.push_str(format!("{:.3} {:.3};", point.x, (point.y - max + boost).min(0.0)).as_str())
        }
        out
    }

    fn sum(&self, floor: &Graph) -> Self {
        Graph {
            points: self
            .points
            .iter()
            .zip(floor.points.iter())
            .map(|(point, fpoint)| {
                Point {
                    x: point.x,
                    y: point.y + fpoint.y,
                }
            })
            .collect()
        }
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
            .collect()
        }
    }

    fn crispyPeak(&self, amp: f64) -> Self {
        Graph {
            points: self.points.iter().map(|point| {
                let fkhz = point.x / 1000.0;
                let sinc = |p: f64| { if p < -PI || p > PI  { 0.0 } else { p.sin() / p } };
                Point {
                    x: point.x,
                    y: point.y + sinc(fkhz * 2.0 - 22.0) * amp,
                }
            }).collect()
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
            // lerp
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

    // take a graph with many x points, and make it have the same x points as another
    // by using 20 simple rms sums for each one.
    // each mean is based on the points to the left and right of the target freq.
    fn resample(&self, target: Vec<f64>, averaging: AverageMode) -> Graph {
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
                        sum += match averaging {
                            AverageMode::Rms => self.point(x).powf(2.0),
                            AverageMode::Mean => self.point(x).abs(),
                        }
                    }
                    Point {
                        x: *point,
                        y: match averaging {
                            AverageMode::Rms => (sum / 20.0).sqrt(),
                            AverageMode::Mean => sum / 20.0,
                        }
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
