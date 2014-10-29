The program contained in ucsf-roc.c has been compiled for linux
(ucsf-roc-linux.exe) and windows (ucsf-roc-sua.exe) both using the GNU
C compiler as follows: gcc -O3 -o ucsf-roc.exe ucsf-roc.c

The usage is simple. Following is the help listing:

csh-% ./ucsf-roc.exe

ROC computes the roc curve for the data given. Output suitable for gnuplot.

roc [options] posdata negdata outprefix
  -ci   ci_percent n       Compute ci_percent confidence interval (e.g. 95) based on n iterations of resampling
  -fpthresh n val1 val2... Compute TP rates at specified FP percentages. Default: 3 1.0 5.0 10.0


The posdata (example in "pos") are the values for the positive
exemplars in your data set (for molecular screening experiments, these
are your actives). By convention, the desire is for these values to be
larger than those for the values in negdata (example in "neg"). If
your data follow the opposite convention, invert their signs. The
program produces both statistics and data files for making graphs. An
example makes this most clear:


csh-% ./ucsf-roc.exe -fpthresh 3 1 5 10 -ci 95 1000 pos neg foo
PosFile has 47 data lines
NegFile has 1797 data lines
ROC Area: 0.72487
Computing confidence intervals: (100)(200)(300)(400)(500)(600)(700)(800)(900)(1000)
Confidence Interval (95.00%):  0.6653 - 0.7791
Average resampled area: 0.725170

The otuput to standard error lets you know how things progress. In
this example, we are asking for true-positive rates at 3 specific
false positive rates (1, 5, and 10%) along with a 95% confidence
interval for the ROC area.

The file foo-rocstats contains the following, which should be self-explanatory:

TP%_at_1.00_FP%: 4.255
TP%_at_5.00_FP%: 19.149
TP%_at_10.00_FP%: 27.660
ROC_Area: 0.72487
ROC_Area_CI_95.00%: 0.6653 - 0.7791

The file foo-roc contains the ROC plot itself. The files foo-pos and
foo-neg are cumulative histograms of the underlying data in the
positive and negative files. The two PDF files are plots resulting
from running Gnuplot on the file GP-cmds (to generate EPS, then
converted to PDF).

NOTE: The file neg5 is simply the neg file replicated 5 times. You
will see that you get the same results using pos with either neg or
neg5. This is because ROC analysis does not skew with meaningless
manipulations of data values or ratios of set sizes. If you are
proposing a new metric, you should make sure that your metric yields
sensible behavior like this!

