Our library is well-tested and we shouldn't be cavalier. However, it needs its logic narrowed as to become more effective and realistic.

- [x] Strip out redundant python DSP pipeline and rely entirely on C. Python is still used to generate coefficients for C, for data visualization, etc.
- [x] Remove all response modeling from our libraries except speaker EQ FFT; the pipeline is the response and we should run sweeps thru our pipeline to generate graphs for review. Also, use said sweeps to measure second and third-order harmonics levels.
- [x] When the above are complete, write a report on the nature and quality of our current speaker configuration system. Shall we create a preset format, etc?

See REPORT.md for the full report.
