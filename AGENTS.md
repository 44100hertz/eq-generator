# Defense is the last defense.

The results have huge peaks that shouldn't be there.
DO: Consider why the peaks exist.
DON'T: Jump straight to capping peaks at a fixed value.

You just pruned low-confidence values, but you expected 300 values and got 1000.
DO: Check if your data is _just that confident_.
DON'T: Simply prune every third value.

# Physics, not math, is the king of the world

You're seeing a speaker needs +0dB at 1khz, and +50dB correction at 16khz.
DO: Assume the microphone or audio interface is rolling off.
DON'T: Try to boost the speaker by 50dB.

The user told you to correct a speaker response, and 10hz is at -20dB.
DO: Ignore 10hz because nobody can even hear it.
DON'T: Analyze how to boost 10hz for an hour, then write a complex extra-steep filter.

Our measurement data is pretty noisy.
DO: Figure out if the "noise" is actually a physical response.
DON'T: Immediately smooth over the noise.

# You're not a student

Don't try to prove your intellect. Don't reinvent the wheel. You're not being tested.
You're being judged only by code quality and functionality.
Code that works perfectly and could be both read and written by a stupid person is ideal.
When it needs to be more complex, of course it should be.

# Definition of "maintainable"

...is how much mental effort it takes to read. That's it.
Minimize cognitive burden, but that doesn't mean minimize functionality.
Functionality of course is tied to purpose. Code with clear purpose is very readable.
Less code is more readable. Dense code is not less code.
Never optimize unless it's hot enough / slow enough to optimize.

# !!!!!! TEST AND ASSERT !!!!!!

The program crashes.
DO: Use step-through debugging, or print debugging in the worst cases.
DON'T: Immediately trace the entire control flow where you _suspect_ the crash is.

The user says something broke in the past (that you didn't touch).
DO: Check the git history.
DON'T: Guess at the exact problem _now_.
Remember: we tend to write code that _looks_ right. Spotting bugs often requires
data to locate it.

There are tests, but something is broken.
DO: Run tests to know what _isn't_ broken.
DON'T: Manually read / check code that already has solid tests on it.

There are no tests, but the equalizer is broken.
DO: Start by writing an FFT analyzer for the equalizer, to know if it works.
DON'T: Identify a random problem, say "Fixed", and call it done.

A function is called in state A, but in theory it could be called in state B which would break it.
DO: Assert that the function is not in state B early on.
DON'T: Trace through every single logical path.

If I say there's a bug, ALWAYS RUN THE FULL PIPELINE, THE REAL PIPELINE END-TO-END.

Don't isolate a part of the pipeline until you've identified the problem by running
what actually runs on the DSP box. THEN you can bisect things down.
