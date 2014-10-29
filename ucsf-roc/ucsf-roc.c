/*
  This implementation of ROC curve, AUC, and confidence-interval
  computation was made public as the result of an ACS National Meeting
  special symposium co-chaired by Anthony Nicholls and Ajay Jain in
  Fall 2007.

  This code may be freely distributed, but changes in the code should
  be clearly marked and described below. It would be helpful if
  significant improvements or error-corrections are sent to
  ajain@jainlab.org so that the central distribution will be kept up
  to date.

  This distribution is maintained at the following URL: http://www.jainlab.org/Public/ucsf-roc.zip

  The code was adapted from Tom Fawcett's pseudo-code in the following
  nice paper:  An Introduction to ROC Analysis. Pattern Recognition Letters 27 (2006) 861-874.

  The area computation was regression-tested against the ACM KDD Cup
  implementation in perf.exe (http://www.sigkdd.org/kddcup/index.php?section=2004&method=soft).
  While this may seem excessive, it turns out that the issue of tied
  values between the scores of positive and negative exemplars can
  have different interpretations. The implementation here follows the
  "trapezoid rule" which gives an unbiased estimate of ROC area in the
  case of ties.
  
  Please consider that a standard legal disclaimer has been inserted
  here. You can thank me in person for not actually putting one in! --Ajay
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define SMALL -10000000
#define BIG 10000000

extern int count_lines(FILE *fd1);
extern void jain_error(char *string);
extern void quicksort(double *array, int p, int r, int *labels);
extern int partition(double *array, int p, int r, int *labels);
extern void swap(double *array, int left, int right, int *labels);
extern double percentile(double *vals, int nvals, double pct);
extern double compute_roc(double *allvals, int *labels, int nallvals, char *path);

int nfpthresh = 3;
double fpthresh[20] = {1.0, 5.0, 10.0};

double trap_area(int x1, int x2, int y1, int y2)
{
  double base = x1-x2;
  double height = (y1+y2)/2.0;

  if (base < 0.0) base = -base;
  return(base*height);
}

void quicksort(double *array, int p, int r, int *labels) { /* p = 0 = left; r = right = nvals -1 on input */
  int pivot;
  if (r > p) {
    pivot = partition(array,p,r,labels);
    quicksort(array,p,pivot-1,labels);
    quicksort(array,pivot+1,r,labels);
  }
}

int partition(double *array, int p, int r, int *labels) {
  int left, right, pivot, pivot_item_carry;
  double pivot_item;
  
  pivot_item = array[p];
  if (labels != NULL) pivot_item_carry = labels[p];

  pivot = left = p;
  right = r;
  while(left < right) {
    while(array[left] <= pivot_item) {
      left++;
      if (left == right) break;
    }
    while(array[right] > pivot_item) {
      right--;
    }
    if (left < right) {
      swap(array, left, right, labels);
    }
  }
  array[p] = array[right];
  array[right] = pivot_item;
  if (labels != NULL) {
    labels[p] = labels[right];
    labels[right] = pivot_item_carry;
  }

  return(right);
}

void swap(double *array, int left, int right, int *labels) {
  double temp = array[left];

  array[left] = array[right];
  array[right] = temp;

  if (labels != NULL) {
    int tempint = labels[left];

    labels[left] = labels[right];
    labels[right] = tempint;
  }
}

double percentile(double *sortvals, int nvals, double pct)
{
  int bin_num;
  double ret_val;

  quicksort(sortvals,0,nvals-1,NULL);

  /* We now have a sorted list */
  bin_num = (int) (pct*(nvals-1.0e-06));
  if (bin_num < 0) bin_num = 0;
  if (bin_num >= nvals) bin_num = nvals-1;

  ret_val = sortvals[bin_num];

  return(ret_val);
}

double random_value(double a, double b)
{
  int r;

  r = rand();

  return((r/((double) RAND_MAX))*(b-a)+a);
}

int main(int argc, char **argv)
{
  int arg,ninput_posvals,ninput_negvals,i;
  int idx, *labels, nallvals;
  double total_area;
  FILE *posfd,*negfd;
  char path[1024];
  int ci_num = 0;
  double ci_percent;
  double min_val, max_val;
  
  double *posvals,*negvals;
  double *origposvals,*orignegvals;
  double *allvals,*sortvals;
  double *bins;
  
  if (argc == 1) {
    fprintf(stderr,"\n");
    fprintf(stderr,"ROC computes the roc curve for the data given. Output suitable for gnuplot.\n\n");
    fprintf(stderr,"roc [options] posdata negdata outprefix\n");
    fprintf(stderr,"  -ci   ci_percent n       Compute ci_percent confidence interval (e.g. 95) based on n iterations of resampling\n");
    fprintf(stderr,"  -fpthresh n val1 val2... Compute TP rates at specified FP percentages. Default: 3 1.0 5.0 10.0\n");
    fprintf(stderr,"\n");
    exit(0);
  }

  for (arg = 1; arg <= argc; ++arg) {
    if (strcmp(argv[arg],"-ci") == 0) {
      sscanf(argv[++arg],"%lf",&ci_percent);
      sscanf(argv[++arg],"%d",&ci_num);
    }
    else if (strcmp(argv[arg],"-fpthresh") == 0) {
      int i;
      sscanf(argv[++arg],"%d",&nfpthresh);
      if (nfpthresh > 20) jain_error("The number of FP thresholds cannot be more than 20.");
      for (i = 0; i < nfpthresh; ++i) {
	sscanf(argv[++arg],"%lf",&(fpthresh[i]));
	if ((fpthresh[i] <= 0.0) || (fpthresh[i] >= 100.0)) {
	  jain_error("FP thresholds must be between 0.0 and 100.0.");
	}
      }
    }
    else {
      break;
    }
  }

  posfd = fopen(argv[arg],"r");
  if (posfd == NULL) {
    fprintf(stderr,"Can't open %s\n",argv[arg]);
    exit(0);
  }

  negfd = fopen(argv[++arg],"r");
  if (negfd == NULL) {
    fprintf(stderr,"Can't open %s\n",argv[arg]);
    exit(0);
  }

  sprintf(path,"%s",argv[++arg]);

  ninput_posvals = count_lines(posfd);
  fprintf(stderr,"PosFile has %d data lines\n",ninput_posvals);
  ninput_negvals = count_lines(negfd);
  fprintf(stderr,"NegFile has %d data lines\n",ninput_negvals);

  posvals = (double *) calloc(ninput_posvals,sizeof(double));
  if (posvals == NULL) jain_error("Can't calloc posvals");
  negvals = (double *) calloc(ninput_negvals,sizeof(double));
  if (negvals == NULL) jain_error("Can't calloc negvals");
  allvals = (double *) calloc(ninput_posvals+ninput_negvals+2,sizeof(double));
  if (posvals == NULL) jain_error("Can't calloc allvals");
  labels = (int *) calloc(ninput_posvals+ninput_negvals+2,sizeof(int));
  if (labels == NULL) jain_error("Can't calloc class labels");

  origposvals = (double *) calloc(ninput_posvals,sizeof(double));
  if (origposvals == NULL) jain_error("Can't calloc origposvals");
  orignegvals = (double *) calloc(ninput_negvals,sizeof(double));
  if (orignegvals == NULL) jain_error("Can't calloc orignegvals");

  for (i = 0; i < ninput_posvals; ++i) {
    fscanf(posfd,"%lf ",&(posvals[i]));
  }
  for (i = 0; i < ninput_negvals; ++i) {
    fscanf(negfd,"%lf ",&(negvals[i]));
  }

  /* Find max and min */
  min_val = BIG;
  max_val = SMALL;

  for (i = 0; i < ninput_posvals; ++i) {
    if (posvals[i] > max_val) max_val = posvals[i];
    if (posvals[i] < min_val) min_val = posvals[i];
    origposvals[i] = posvals[i];
  }
  for (i = 0; i < ninput_negvals; ++i) {
    if (negvals[i] > max_val) max_val = negvals[i];
    if (negvals[i] < min_val) min_val = negvals[i];
    orignegvals[i] = negvals[i];
  }

  /* Copy the data to the allvals array and mark with labels */
  idx = 0;
  for (i = 0; i < ninput_posvals; ++i) {
    allvals[idx] = posvals[i];
    labels[idx] = 1;		/* we are a positive */
    ++idx;
  }
  for (i = 0; i < ninput_negvals; ++i) {
    allvals[idx] = negvals[i];
    labels[idx] = 0;		/* we are a negative */
    ++idx;
  }
  nallvals = idx;

  total_area = compute_roc(allvals,labels,nallvals,path);
  fprintf(stderr,"ROC Area: %.5lf\n",total_area);

  if (ci_num > 0) {		/* We need to compute a confidence interval */
    double *resampled_vals, *areas, ci_low,ci_high, mean_area;
    int try;
    FILE *rocstatsfd;

    resampled_vals = (double *) calloc(nallvals,sizeof(double));
    if (resampled_vals == NULL) jain_error("Can't calloc resampled_vals");

    areas = (double *) calloc(ci_num,sizeof(double));
    if (areas == NULL) jain_error("Can't calloc areas");

    fprintf(stderr,"Computing confidence intervals: ");
    mean_area = 0.0;
    for (try = 0; try < ci_num; ++try) {
      if ((try+1)%(ci_num/10) == 0) fprintf(stderr,"(%d)",try+1);
      /* Resample the positive and negative values */
      idx = 0;
      for (i = 0; i < ninput_posvals; ++i) {
	resampled_vals[idx] = posvals[(int) random_value(0.0,(ninput_posvals - 1.0e-10))];
	labels[idx] = 1;		/* we are a positive */
	++idx;
      }
      for (i = 0; i < ninput_negvals; ++i) {
	resampled_vals[idx] = negvals[(int) random_value(0.0,(ninput_negvals - 1.0e-10))];
	labels[idx] = 0;		/* we are a negative */
	++idx;
      }
      areas[try] = compute_roc(resampled_vals,labels,nallvals,NULL);
      mean_area += areas[try];
    }
    mean_area = mean_area/ci_num;
    ci_low = percentile(areas,ci_num,(100.0-ci_percent)/200.0);
    ci_high = percentile(areas,ci_num,(100.0+ci_percent)/200.0);
    fprintf(stderr,"\nConfidence Interval (%.2lf%%):  %.4lf - %.4lf\n",ci_percent,ci_low,ci_high);
    fprintf(stderr,"Average resampled area: %.6lf\n",mean_area);
    sprintf(path,"%s-rocstats",path);
    rocstatsfd = fopen(path,"a");
    if (rocstatsfd == NULL) {
      jain_error("Can't open output roc stats file.");
    }
    fprintf(rocstatsfd,"ROC_Area_CI_%.2lf%%: %.4lf - %.4lf\n",ci_percent,ci_low,ci_high);
    fclose(rocstatsfd);
  }
}

int count_lines(FILE *fd1)
{
  int n = 0;

  while (fscanf(fd1,"%*[^\n\r]%*[\n\r]") != EOF) ++n;
  rewind(fd1);
  return(n);
}

void jain_error(char *string)
{
  fprintf(stderr,"%s\n",string);
  exit(0);
}

double compute_roc(double *allvals, int *labels, int nallvals, char *prefix)
{
  int i, tp, fp, old_tp, old_fp;
  double n, p, old_val;
  double total_area;
  FILE *rocfd, *posfd, *negfd, *rocstatsfd;

  if (prefix != NULL) {		/* We are dumping the graphs */
    char path[1024];

    sprintf(path,"%s-roc",prefix);
    rocfd = fopen(path,"w");
    if (rocfd == NULL) {
      jain_error("Can't open output roc plot file.");
    }

    sprintf(path,"%s-pos",prefix);
    posfd = fopen(path,"w");
    if (posfd == NULL) {
      jain_error("Can't open output pos hist plot file.");
    }

    sprintf(path,"%s-neg",prefix);
    negfd = fopen(path,"w");
    if (negfd == NULL) {
      jain_error("Can't open output neg hist plot file.");
    }

    sprintf(path,"%s-rocstats",prefix);
    rocstatsfd = fopen(path,"w");
    if (rocstatsfd == NULL) {
      jain_error("Can't open output roc stats file.");
    }
  }

  n = p = 0.0;
  for (i = 0; i < nallvals; ++i) {
    if (labels[i] == 1) {
      p += 1.0;
    }
    if (labels[i] == 0) {
      n += 1.0;
    }
  }

  quicksort(allvals,0,nallvals-1,labels);  /* sort the values fast */

  total_area = 0.0;

  /* We start the ROC curve at the lower left corner corresponding to
     the highest threshold possible, which will yield 0 TP and 0 FP */
  old_val = -100000.0;
  i = nallvals-1;
  tp = fp = 0;
  old_tp = old_fp = 0;
  while (i >= 0) {
    if (allvals[i] != old_val) { /* We plot when the threshold actually changes.
				    All ties have been absorbed into a single change.
				    We will see a horizontal segment if only the FP rate changed.
				    We will see a vertical segment if only the TP rate changed.
				    We will see a sloped line if both did. */
      total_area += trap_area(fp,old_fp,tp,old_tp);
      if (prefix != NULL) {
	fprintf(rocfd,"%lf %lf\n",fp/n,tp/p); /* ROC point */
	fprintf(posfd,"%lf %lf\n",allvals[i],1.0-(tp/p)); /* Cumulative histogram of positives */
	fprintf(negfd,"%lf %lf\n",allvals[i],1.0-(fp/n)); /* Cumulative histogram of negatives */
      }

      if ((prefix != NULL) && (nfpthresh > 0)) { /* need to find TP vals for specific FP vals */
	double fpval,oldfprate,fprate,oldtprate,tprate,tpval,w1,w2,gap,percent;
	int k;

	for (k = 0; k < nfpthresh; ++k) {
	  percent = fpthresh[k];
	  fpval = percent/100.0;
	  oldfprate = old_fp/n;
	  fprate = fp/n;
	  oldtprate = old_tp/p;
	  tprate = tp/p;
	  if ((oldfprate < fpval) && (fpval <= fprate)) { /* we've bracketed me */
	    fprintf(rocstatsfd,"TP%%_at_%.2lf_FP%%: ",fpval*100.0);
	    w1 = fpval - oldfprate;
	    w2 = fprate - fpval;
	    if (w1 == 0.0) {
	      tpval = oldtprate;
	    }
	    else if (w2 == 0.0) {
	      tpval = tprate;
	    }
	    else {		/* Linear interpolation */
	      gap = fprate - oldfprate; /* must be > zero if we are here */
	      tpval = oldtprate + (w1/gap)*(tprate-oldtprate);
	    }
	    fprintf(rocstatsfd,"%.3lf\n",100.0*tpval);
	  }
	}
      }

      old_val = allvals[i];
      old_fp = fp;
      old_tp = tp;
    }
    if (labels[i] == 1) {
      ++tp;
    }
    if (labels[i] == 0) {
      ++fp;
    }
    i = i-1;
  }
  total_area += trap_area(fp,old_fp,tp,old_tp);
  total_area = total_area/(tp*fp);
  if (prefix != NULL) {		/* dump last values and close files */
    fprintf(rocfd,"%lf %lf\n",fp/n,tp/p);
    fprintf(posfd,"%lf %lf\n",allvals[0],1.0-(tp/p));
    fprintf(negfd,"%lf %lf\n",allvals[0],1.0-(fp/n));
    fclose(rocfd);
    fclose(posfd);
    fclose(negfd);
  }
  if (prefix != NULL) {
    fprintf(rocstatsfd,"ROC_Area: %.5lf\n",total_area);
    fclose(rocstatsfd);
  }
  return(total_area);
}
