use std::fs::{self, File};
use std::io::{Write};

const BOTTOM_F: f64 = 50.0;
const TOP_F: f64 = 14000.0;

fn main() {
    let measurements = [
        "spectrum1.txt",
        "spectrum2.txt",
    ].iter().map(|file| Graph::from_audacity_spectogram(file).resample(eq_points())).collect();
    let target = Graph::from_audacity_spectogram("target.txt").resample(eq_points());
    let rms = rms(measurements);
    let comp_raw = Graph {
        points: eq_points().into_iter().map(|freq| {
            Point {
                x: freq,
                y: target.point(freq) - rms.point(freq),
            }
        }).collect()
    };
    let comp_max = comp_raw.max();
    let comp = comp_raw.map(|x: f64| x - comp_max);
    let mut outfile = File::create("autoeq.csv").unwrap();
    outfile.write_all(comp.to_jamesdsp_eq().as_bytes()).unwrap();
}

// Assumes that graphs all have the same indices...
// And that all y values are negative, no less than 120
fn rms(graphs: Vec<Graph>) -> Graph {
    Graph {
        points: (0..graphs[0].points.len()).map(|i| {
            let sum_squares = graphs.iter().fold(0.0, |acc, graph| {
                acc + db_to_linear(graph.points[i].y)
            });
            let rms_sum = (sum_squares / graphs.len() as f64).sqrt();
            Point {
                x: graphs[0].points[i].x,
                y: linear_to_db(rms_sum),
            }
        }).collect()
    }
}

fn db_to_linear(db: f64) -> f64 { 10.0_f64.powf(db / 10.0) }
fn linear_to_db(linear: f64) -> f64 { linear.log10() * 10.0 }

fn eq_points() -> Vec<f64> {
    // how many points per octave?
    fn density(freq: f64) -> f64 {
        // based on a rough tracing of equal loudness countour
        let equal_loudness = Graph {
            points: vec![
                Point{x: 20.0, y: 109.0},
                Point{x: 80.0, y: 82.0},
                Point{x: 400.0, y: 62.0},
                Point{x: 1000.0, y: 60.0},
                Point{x: 1500.0, y: 64.0},
                Point{x: 2500.0, y: 57.0},
                Point{x: 4000.0, y: 57.0},
                Point{x: 8500.0, y: 73.0},
                Point{x: 15000.0, y: 72.0},
                Point{x: 19000.0, y: 68.0},
                Point{x: 30000.0, y: 130.0},
            ],
        };
        let dbs = equal_loudness.point(freq);
        400.0 / dbs
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

// an array of points, frequency to db
struct Graph {
    points: Vec<Point>,
}

impl Graph {
    fn from_audacity_spectogram(path: &str) -> Self {
        let file = fs::read_to_string(path).expect("Could not read file");
        let mut lines = file.split("\n");
        // skip over headers
        lines.next().expect("Empty file!");
        Graph {
            points: lines.filter(|line| !line.is_empty())
                        .map(|line| {
                let mut fields = line.split("\t");
                Point {
                    x: fields.next().unwrap().parse().unwrap(),
                    y: fields.next().unwrap().parse().unwrap(),
                }
            }).collect()
        }
    }

    fn to_jamesdsp_eq(&self) -> String {
        let mut out = String::new();
        for point in self.points.iter() {
            out.push_str(format!("{:.3}\t{:.3}\n", point.x, point.y).as_str())
        }
        out
    }

    fn point(&self, x: f64) -> f64 {
        if let Some(lower) = self.points.windows(2).enumerate().find(|(_, fs)| fs[0].x <= x && fs[1].x >= x) {
            let lower = lower.0;
            let p1 = self.points[lower];
            let p2 = self.points[lower+1];
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
            points: self.points.iter().map(|point| {
                Point {
                    x: point.x,
                    y: tform(point.y),
                }
            }).collect()
        }
    }

    // take a graph with many x points, and make it have the same x points as another
    // by using 20 simple rms sums for each one.
    // each mean is based on the points to the left and right of the target freq.
    fn resample(&self, target: Vec<f64>) -> Graph {
        Graph {
            points: target.iter().enumerate().map(|(i, point)| {
                let left = if i > 0 { target[i-1] } else { *point };
                let right = target.get(i+1).unwrap_or(point);

                let x1 = (point + left) / 2.0;
                let x2 = (point + right) / 2.0;

                let mut sum = 0_f64;
                for i in 0..20 {
                    let x = x1 + (x2 - x1) * i as f64 / 20.0;
                    sum += db_to_linear(self.point(x)).powf(2.0);
                }
                Point {
                    x: *point,
                    y: linear_to_db((sum / 20.0).sqrt()),
                }
            }).collect()
        }
    }

    fn max(&self) -> f64 {
        self.points.iter().map(|point| point.y).reduce(|a,b| a.max(b)).unwrap()
    }
}

#[derive(Clone, Copy, Debug)]
struct Point {
    x: f64,
    y: f64,
}
