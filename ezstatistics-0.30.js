"use strict";

/* Library version */
var VERSION = "0.30";
/*****************************************************

Functions for statistical tests.

******************************************************/

/**
    R-distribution for Wilcoxon rank-sum, two-tailed.
*/
var Rc005_2t = [
    [ [5,16], [6,18], [6,21], [7,23], [7,26], [8,28], [8,31], [9,33] ],
    [ [6,18], [11,25], [12,28], [12,32], [13,35], [14,38], [15,41], [16,44] ],
    [ [6,21], [12,28], [18,37], [19,41], [20,45], [21,49], [22,53], [24,56] ],
    [ [7,23], [12,32], [19,41], [26,5], [28,56], [29,61], [31,65], [32,70] ],
    [ [7,26], [13,35], [20,45], [28,56], [37,68], [39,73], [41,78], [43,83] ],
    [ [8,28], [14,38], [21,49], [29,61], [39,73], [49,87], [51,93], [54,98] ],
    [ [8,31], [15,41], [22,53], [31,65], [41,78], [51,93], [63,108], [66,114] ],
    [ [9,33], [16,44], [24,56], [32,70], [43,83], [54,98], [66,114], [79,131] ]
];

/**
    R-distribution for Wilcoxon rank-sum, one-tailed.
*/
var Rc005_1t = [
    [ [6,15], [7,17], [7,20], [8,22], [9,24], [9,27], [10,29], [11,31] ],
    [ [7,17], [12,24], [13,27], [14,30], [15,33], [16,36], [17,39], [18,42] ],
    [ [7,20], [13,27], [19,36], [20,40], [22,43], [24,46], [25,50], [26,54] ],
    [ [8,22], [14,30], [20,40], [28,50], [30,54], [32,58], [33,63], [35,67] ],
    [ [9,24], [15,33], [22,43], [30,54], [39,66], [41,71], [43,76], [46,80] ],
    [ [9,27], [16,36], [24,46], [32,58], [41,71], [52,84], [54,90], [57,95] ],
    [ [10,29], [17,39], [25,50], [33,63], [43,76], [54,90], [66,105], [69,111] ],
    [ [11,31], [18,42], [26,54], [35,67], [46,80], [57,95], [69,111], [83,127] ]
];

/**
    W-distribution for Wilcoxon signed-ranks, two-tailed.
    http://users.stat.ufl.edu/~winner/tables/wilcox_signrank.pdf
*/
var Wc005_2t = [-1,0,2,3,5,8,10,13,17,21,25,29,34,40,46,52,58,65,73,81,89,98,107,116,126,137];

/**
    W-distribution for Wilcoxon signed-ranks, one-tailed.
*/
var Wc005_1t = [0,2,3,5,8,10,13,17,21,25,30,35,41,47,53,60,67,75,83,91,100,110,119,130,140,151];

/**
    Calculates a confidence interval for a sample.
    
    Params:
        v1: The sample
        conf: confidence level (0.90, 0.95 or 0.98)
*/
function confidence_interval(v1, conf=0.9) {
    // Check if input params are valid
    if (!Array.isArray(v1)) {
        throw "Parameter is not a valid array";
    }
    if (conf == 0.9 || conf == 0.95 || conf == 0.98) {
        //OK
    }
    else {
        throw "Confidence level must be 0.90, 0.95 or 0.98";
    }
    
    // Sample mean
    let mean = jStat.mean(v1);
    // Sample size
    let n = v1.length;
    // Standard deviations for the samples
    let std = jStat.stdev(v1, true);
    // Standard error
    let stderr = std / Math.sqrt(n);
    // Degrees of freedom
    let DF = n - 1;
    // Alpha
    let alpha = 0.1;
    if (conf == 0.9) alpha = 0.05;
    if (conf == 0.95) alpha = 0.025;
    if (conf == 0.98) alpha = 0.01;
    
    if (n >= 30) {
        // Larg n -> use normal distribution as approximation
        let Z = jStat.normal.inv(1 - alpha, 0, 1);
        // Confidence interval
        let boundary = Z * stderr;
        // Range
        let lower = mean - boundary;
        let upper = mean + boundary;
        return [boundary, lower, upper, alpha];
    }
    else {
        // Small n -> use t-distribution
        let T = jStat.studentt.inv(1 - alpha, DF);
        // Confidence interval
        let boundary = T * stderr;
        // Range
        let lower = mean - boundary;
        let upper = mean + boundary;
        return [boundary, lower, upper, alpha];
    }
}

/**
    Performs a t-test.
    
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
        type: 1 for independent equal variances, 2 for independent unequal variances, 3 for paired
        diff: difference between samples to account for
*/
function ttest(v1, v2, sides=2, tail=1, alpha=0.05, type=1, diff=0) {
    // Check if input params are valid
    if (!Array.isArray(v1) || !Array.isArray(v2)) {
        throw "Parameter is not a valid array";
    }
    if (sides < 1 || sides > 2) {
        throw "Sides must be 1 or 2";
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    if (type < 1 || type > 3) {
        throw "Type must be 1, 2 or 3";
    }   
    
    if (type == 1) {
        // Independent, equal variance
        return ttest_ind(v1, v2, sides, tail, alpha, true, diff);
    }
    else if (type == 2) {
        // Independent, unequal variance
        return ttest_ind(v1, v2, sides, tail, alpha, false, diff);
    }
    else if (type == 3) {
        if (v1.length != v2.length) {
            throw "Sample sizes must be equal";
        }
        // Paired
        return ttest_paired(v1, v2, sides, tail, alpha, diff);
    }
}

/**
    T-test for two independent samples.
  
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
        eq_var: true for equal variances test, false for unequal variances test
        diff: difference between samples to account for
*/
function ttest_ind(v1, v2, sides=2, tail=1, alpha=0.05, eq_var=true, diff=0) {   
    // Sample means
    let mean1 = jStat.mean(v1);
    let mean2 = jStat.mean(v2);
    // Sample sizes
    let n1 = v1.length;
    let n2 = v2.length;
    // Standard deviations for the samples
    let std1 = jStat.stdev(v1, true);
    let std2 = jStat.stdev(v2, true);
    let std1sq = Math.pow(std1, 2);
    let std2sq = Math.pow(std2, 2);
    // Equal variance
    
    let DF = 0;
    let o = 0;
    if (eq_var) {
        // Calculate pooled variance
        let sp = Math.sqrt( ((n1 - 1) * std1sq + (n2 - 1) * std2sq) / (n1 + n2 - 2) );
        // Calculate approximated standard error
        o = sp * Math.sqrt(1 / n1 + 1 / n2);
        // Degrees of freedom
        DF = n1 + n2 - 2;
    }
    // Unequal variance
    else {
        // Calculate approximated standard error
        o = Math.sqrt(std1sq / n1 + std2sq / n2);
        // Degrees of freedom
        let n = Math.pow( std1sq / n1 + std2sq / n2, 2);
        let d = Math.pow(std1sq / n1, 2) / (n1 - 1) + Math.pow(std2sq / n2, 2) / (n2 - 1);
        DF = n / d;
    }
    
    
    // Find critical T-value
    let Tc = 0;
    if (sides == 2) {
        Tc = jStat.studentt.inv(1 - alpha / 2, DF);
    }
    if (sides == 1 && tail == 1) {
        Tc = jStat.studentt.inv(alpha, DF);
    }
    if (sides == 1 && tail == 2) {
        Tc = jStat.studentt.inv(1 - alpha, DF);
    }
    // Difference between sample means
    let Xdiff = mean1 - mean2 - diff;
    // T-score
    let T = Xdiff / o;
    // P-value
    let P = 0;
    if (sides == 2) {
        P = jStat.studentt.cdf(T, DF) * sides;
        if (T > 0) {
            P = (1 - jStat.studentt.cdf(T, DF)) * sides;
        }
    }
    if (sides == 1 && tail == 1) {
        P = jStat.studentt.cdf(T, DF) * sides;
    }
    if (sides == 1 && tail == 2) {
        P = (1 - jStat.studentt.cdf(T, DF)) * sides;
    }
    
    return [P, T, Tc];
}

/**
    T-test for two paired samples.
    
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
        diff: difference between samples to account for
*/
function ttest_paired(v1, v2, sides=2, tail=1, alpha=0.05, diff=0) {
    // Sample size
    let n = v1.length;
    
    // Sample means
    let mean1 = jStat.mean(v1);
    let mean2 = jStat.mean(v2);
    // Calculate differences and squared differences between samples
    let sum_diff = 0;
    let sum_diff_sq = 0;
    for (let i in v1) {
        diff = v2[i] - v1[i];
        sum_diff += diff;
        sum_diff_sq += Math.pow(diff, 2);
    }
    // Calculate standard deviation of the differences
    let s = Math.sqrt( (sum_diff_sq - Math.pow(sum_diff, 2) / n) / (n - 1) );
    // Calculate mean difference
    let dmean = sum_diff / n;
    // T-score
    let T = dmean / (s / Math.sqrt(n)) * -1;
    // Degrees of freedom
    let DF = n - 1;
    // Find critical T-value
    let Tc = 0;
    if (sides == 2) {
        Tc = jStat.studentt.inv(1 - alpha / 2, DF);
    }
    if (sides == 1 && tail == 1) {
        Tc = jStat.studentt.inv(alpha, DF);
    }
    if (sides == 1 && tail == 2) {
        Tc = jStat.studentt.inv(1 - alpha, DF);
    }
    
    // P-value
    let P = 0;
    if (sides == 2) {
        P = jStat.studentt.cdf(T, DF) * sides;
        if (T > 0) {
            P = (1 - jStat.studentt.cdf(T, DF)) * sides;
        }
    }
    if (sides == 1 && tail == 1) {
        P = jStat.studentt.cdf(T, DF) * sides;
    }
    if (sides == 1 && tail == 2) {
        P = (1 - jStat.studentt.cdf(T, DF)) * sides;
    }
    
    return [P, T, Tc];
}

/**
    T-test single sample for the mean.
  
    Params:
        v1: Sample
        mean: Mean to compare with
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
*/
function ttest_single(v1, mean, sides=2, tail=1, alpha=0.05) {
    // Sample mean
    let mean1 = jStat.mean(v1);
    let n = v1.length;
    
    // Standard deviations for the sample
    let std1 = jStat.stdev(v1, true);
    let stderr1 = std1 / Math.sqrt(n);
    
    // Degrees of freedom
    let DF = n - 1;
    
    // Calculate T-score
    let T = (mean1 - mean) / stderr1;
    // Find critical T-value
    let Tc = 0;
    if (sides == 2) {
        Tc = jStat.studentt.inv(1 - alpha / 2, DF);
    }
    if (sides == 1 && tail == 1) {
        Tc = jStat.studentt.inv(alpha, DF);
    }
    if (sides == 1 && tail == 2) {
        Tc = jStat.studentt.inv(1 - alpha, DF);
    }
    
    // P-value
    let P = 0;
    if (sides == 2) {
        P = jStat.studentt.cdf(T, DF) * sides;
        if (T > 0) {
            P = (1 - jStat.studentt.cdf(T, DF)) * sides;
        }
    }
    if (sides == 1 && tail == 1) {
        P = jStat.studentt.cdf(T, DF) * sides;
    }
    if (sides == 1 && tail == 2) {
        P = (1 - jStat.studentt.cdf(T, DF)) * sides;
    }
    
    return [P, T, Tc];
}

/**
    F-test for equal variance.
    
    Params:
        v1: Sample 1
        v2: Sample 2
        alpha: Significance level
*/
function ftest(v1, v2, alpha=0.05) {
    // Check if input params is valid
    if (!Array.isArray(v1) || !Array.isArray(v2)) {
        throw "Parameter is not a valid array";
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    
    // Calculate sample size and standard deviation.
    // Note that sample A must have higher standard deviation than Sample B.
    let n1o = v1.length;
    let n2o = v2.length;
    let std1o = jStat.stdev(v1, true);
    let std2o = jStat.stdev(v2, true);

    let n1 = n2o;
    let n2 = n1o;
    let std1 = std2o;
    let std2 = std1o;
    if (std1o > std2o) {
        let n1 = n1o;
        let n2 = n2o;
        let std1 = std1o;
        let std2 = std2o;
    }

    // Calculate F-score
    let F = Math.pow(std1, 2) / Math.pow(std2, 2);
    // Degrees of freedom
    let DF1 = n1 - 1;
    let DF2 = n2 - 1;
    // Calculate critical F-value
    let Fc = jStat.centralF.inv(1 - alpha / 2, DF1, DF2);
    // P-value
    let P = jStat.centralF.cdf(F, DF1, DF2) * 2;
    if (F > 0) {
        P = (1 - jStat.centralF.cdf(F, DF1, DF2)) * 2;
    }
    
    return [P, F, Fc];
}

/**
  Performs a One-way ANOVA.
  
  Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
*/
function anova(arrs, alpha=0.05) {
    // Check if input params is valid
    if (!Array.isArray(arrs)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    // Number of samples
    let k = arrs.length;
    for (let i = 0; i < k; i++) {
        // Check each sample array
        if (!Array.isArray(arrs[i])) {
            throw "Array " + i + " is not a valid array";
        }
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    if (k < 3) {
        throw "Requires three or more samples";
    }
    
    // Calculate total sum, total sum of squares and total n
    let n_tot = 0;
    let sum_tot = 0;
    let sum_sq_tot = 0;
    for (let s = 0; s < k; s++) {
        let v = arrs[s];
        n_tot += v.length;
        // Sum and squared sum
        sum_tot += jStat.sum(v);
        sum_sq_tot += jStat.sumsqrd(v);
    }
    
    let SST = sum_sq_tot - Math.pow(sum_tot, 2) / n_tot;
    
    // Calculate grand mean
    let grand_mean = 0;
    for (let s = 0; s < k; s++) {
        let v = arrs[s];
        grand_mean += jStat.mean(v) * v.length;
    }
    grand_mean /= n_tot;
    
    // Calculate SSB and SSW
    let SSB = 0;
    for (let s = 0; s < k; s++) {
        let v = arrs[s];
        SSB += v.length * Math.pow(jStat.mean(v) - grand_mean, 2);
    }
    let SSW = SST - SSB;
    
    // Calculate F-score
    let MSB = SSB / (k - 1);
    let MSW = SSW / (n_tot - k);
    let F = MSB / MSW;
    
    // Degrees of freedom
    let DF1 = k - 1;
    let DF2 = n_tot - k;
    
    // Calculate critical F-value
    let Fc = jStat.centralF.inv(1 - alpha, DF1, DF2);
    // P-value
    let P = jStat.centralF.cdf(F, DF1, DF2);
    if (F > 0) {
        P = 1 - jStat.centralF.cdf(F, DF1, DF2);
    }
    
    return [P, MSB, MSW, DF1, DF2, F, Fc];
}

/**
  Performs a Repeated Measures ANOVA.
  See:  
        https://statistics.laerd.com/statistical-guides/repeated-measures-anova-statistical-guide.php
        http://www.real-statistics.com/anova-repeated-measures/repeated-measures-anova-tool/
  
  Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
*/
function rm_anova(arrs, alpha=0.05) {
    // Check if input params is valid
    if (!Array.isArray(arrs)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    // Number of samples
    let k = arrs.length;
    for (let i = 0; i < k; i++) {
        // Check each sample array
        if (!Array.isArray(arrs[i])) {
            throw "Array " + i + " is not a valid array";
        }
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    if (k < 3) {
        throw "Requires three or more samples";
    }
    let n = arrs[0].length;
    for (let i = 0; i < k; i++) {
        let nc = arrs[i].length;
        if (nc != n) {
            throw "Requires equal sample sizes";
        }
    }
    
    // Calculate total sum of squares and total n
    let totN = 0;
    let totSum = 0;
    let totSumSq = 0;
    for (let i = 0; i < k; i++) {
        let vc = arrs[i];
        totN += vc.length;
        for (let s = 0; s < vc.length; s++) {
            totSum += vc[s];
            totSumSq += Math.pow(vc[s], 2);
        }
    }
    let SST = totSumSq - Math.pow(totSum, 2) / totN;
    
    // Calculate grande mean
    let grand_mean = 0;
    for (let i = 0; i < k; i++) {
        let vc = arrs[i];
        grand_mean += jStat.mean(vc) * vc.length;
    }
    grand_mean /= totN;
    
    // Calculate SSB and SSW
    let SStime = 0;
    for (let i = 0; i < k; i++) {
        let vc = arrs[i];
        SStime += vc.length * Math.pow(jStat.mean(vc) - grand_mean, 2);
    }
    let SSW = SST - SStime;
    
    // Calculate subject means
    let sub_means = new Array(n).fill(0);
    for (let s = 0; s < n; s++) {
        for (let i = 0; i < k; i++) {
            let vc = arrs[i];
            sub_means[s] += vc[s];
        }
        sub_means[s] /= k;
    }
    // Calculate SSsubject
    let SSsubj = 0;
    for (let s = 0; s < n; s++) {
        SSsubj += Math.pow(sub_means[s] - grand_mean, 2);
    }
    SSsubj *= k;
    
    // Calculate SSerror
    let SSerror = SSW - SSsubj;
    
    // Calculate F-score
    let MStime = SStime / (k - 1);
    let MSsubj = SSsubj / (n - 1);
    let MSerror = SSerror / ( (n - 1) * (k - 1) );
    let F = MStime / MSerror;
    
    // Degrees of freedome
    let DF1 = k - 1;
    let DF2 = (n - 1) * (k - 1);
    
    // Calculate critical F-value
    let Fc = jStat.centralF.inv(1 - alpha, DF1, DF2);
    // P-value
    let P = jStat.centralF.cdf(F, DF1, DF2);
    if (F > 0) {
        P = 1 - jStat.centralF.cdf(F, DF1, DF2);
    }
    
    return [P, SSW, MSerror, DF1, DF2, F, Fc];
}

/**
    Performs Scheffe's post-test.
    
    Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
        correction: Significance level correction (none or Bonferroni)
        res: ANOVA results
*/
function scheffes_posttest(arrs, alpha=0.05, correction="none", res) {
    // Number of samples
    let k = arrs.length;
    
    // Correction
    if (correction == "bonferroni") {
        let m = k * (k - 1) / 2;
        alpha = alpha / m;
    }
    
    // Get ANOVA results
    let MSW = res[2];
    let DF1 = res[3];
    let DF2 = res[4];
            
    // Pair-wise comparison
    let P_vals = [];
    for (let i = 0; i < k; i++) {
        for (let j = i + 1; j < k; j++) {
            // Samples
            let v1 = arrs[i];
            let v2 = arrs[j];
            let mean1 = jStat.mean(v1);
            let mean2 = jStat.mean(v2);
            // Harmonic mean of sample sizes
            let n = 2 / (1/v1.length + 1/v2.length);
            // Tukey's Q-value
            let Q = Math.abs(mean1 - mean2) / Math.sqrt(MSW / n);
            // Calculate F-value
            let T = Q / Math.sqrt(2);
            let F = Math.pow(T, 2) / DF1;
            // Find critical F
            let Fc = jStat.centralF.inv(1 - alpha, DF1, DF2);
            // P-value
            let P = jStat.centralF.cdf(F, DF1, DF2);
            if (F > 0) {
                P = (1 - jStat.centralF.cdf(F, DF1, DF2));
            }
            // Result
            let postres = [ i + 1, j + 1, P, alpha, F, Fc, "F"];
            P_vals.push(postres);
        }
    }
    
    return P_vals;
}

/**
    Performs Tukey's HSD post-test.
    
    Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
        correction: Significance level correction (none or Bonferroni)
        res: ANOVA results
*/
function tukeys_posttest(arrs, alpha=0.05, correction="none", res) {
    // Number of samples
    let k = arrs.length;
    
    // Correction
    if (correction == "bonferroni") {
        let m = k * (k - 1) / 2;
        alpha = alpha / m;
    }
    
    // Get ANOVA results
    let MSerr = res[2];
    let DF = res[4];
    
    // Pair-wise comparison
    let P_vals = [];
    for (let i = 0; i < k; i++) {
        for (let j = i + 1; j < k; j++) {
            // Samples
            let v1 = arrs[i];
            let v2 = arrs[j];
            let mean1 = jStat.mean(v1);
            let mean2 = jStat.mean(v2);
            // Harmonic mean of sample sizes
            let n = 2 / (1/v1.length + 1/v2.length);
            // Calculate Q-value
            let Q = Math.abs(mean1 - mean2) / Math.sqrt(MSerr / n);
            // Find critical Q
            let Qc = jStat.tukey.inv(1 - alpha, k, DF);
            // P-value
            let P = 1 - jStat.tukey.cdf(Q, k, DF);
            // Result
            let res = [ i + 1, j + 1, P, alpha, Q, Qc, "Q"];
            P_vals.push(res);
        }
    }
    
    return P_vals;
}

/**
    Performs Tukey's HSD post-test.
    
    Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
        correction: Significance level correction (none or Bonferroni)
        res: ANOVA results
*/
function tukeys_rm_posttest(arrs, alpha=0.05, correction="none", res) {
    // Number of samples
    let k = arrs.length;
    
    // Correction
    if (correction == "bonferroni") {
        let m = k * (k - 1) / 2;
        alpha = alpha / m;
    }
    
    // Get ANOVA results
    let MSW = res[2];
    let DF = res[4];
    
    // Pair-wise comparison
    let P_vals = [];
    for (let i = 0; i < k; i++) {
        for (let j = i + 1; j < k; j++) {
            // Samples
            let v1 = arrs[i];
            let v2 = arrs[j];
            let mean1 = jStat.mean(v1);
            let mean2 = jStat.mean(v2);
            // Harmonic mean of sample sizes
            let n = 2 / (1/v1.length + 1/v2.length);
            // Calculate Q-value
            let Q = Math.abs(mean1 - mean2) / Math.sqrt(MSW / n);
            // Find critical Q
            let Qc = jStat.tukey.inv(1 - alpha, k, DF);
            // P-value
            let P = 1 - jStat.tukey.cdf(Q, k, DF);
            // Result
            let res = [ i + 1, j + 1, P, alpha, Q, Qc, "Q"];
            P_vals.push(res);
        }
    }
    
    return P_vals;
}

/**
  Performs a Shapiro-Wilk test to check if a sample is normally distributed.
  See:
        http://www.real-statistics.com/tests-normality-and-symmetry/statistical-tests-normality-symmetry/shapiro-wilk-expanded-test/
        http://www.real-statistics.com/descriptive-statistics/symmetry-skewness-kurtosis/
  
  Params:
        v: Sample array
        alpha: Significance level
*/
function shapiro_wilk(v, alpha=0.05) {
    // Check if input params is valid
    if (!Array.isArray(v)) {
        throw "Parameter is not a valid array";
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    
    // Sample size
    let n = v.length;
    
    // Sort array
    v.sort((a, b) => a - b);
    
    // Calculate mi
    let mi = [];
    for (let i = 0; i < n; i++) {
        let p = ((i + 1) - 0.375) / (n + 0.25);
        mi.push( jStat.normal.inv(p, 0, 1) );
    }
    
    // Calculate m
    let m = 0;
    for (let i = 0; i < n; i++) {
        m += Math.pow(mi[i], 2);
    }
    
    // Calculate ai
    let u = 1 / Math.sqrt(n);
    // Create empty array
    let ai = [];
    for (let i = 0; i < n; i++) {
        ai.push(0);
    }
    ai[n - 1] = -2.70605 * Math.pow(u, 5) + 4.434685 * Math.pow(u, 4) - 2.071190 * Math.pow(u, 3) - 0.147981 * Math.pow(u, 2) + 0.221157 * u + mi[n - 1] * Math.pow(m, -0.5);
    ai[n - 2] = -3.582633 * Math.pow(u, 5) + 5.682633 * Math.pow(u, 4) - 1.752461 * Math.pow(u, 3) - 0.293762 * Math.pow(u, 2) + 0.042981 * u + mi[n - 2] * Math.pow(m, -0.5);
    for (let i = 0; i < n - 2; i++) {
        let ci = i + 1;
        if (ci == 1) {
            ai[i] = -ai[n - 1];
        }
        else if (ci == 2) {
            ai[i] = -ai[n - 2];
        }
        else {
            let eps = (m - 2.0 * Math.pow(mi[n - 1], 2) - 2.0 * Math.pow(mi[n - 2], 2)) / (1.0 - 2.0 * Math.pow(ai[n - 1], 2) - 2.0 * Math.pow(ai[n - 2], 2));
            ai[i] = mi[i] / Math.sqrt(eps);
        }
    }
    
    // Sample mean
    let mean = jStat.mean(v);
    // Sample standard deviation
    let stdev = jStat.stdev(v, true);
    
    // Calculate W
    let sum_A = 0;
    let sum_B = 0;
    for (let i = 0; i < n; i++) {
        sum_A += ai[i] * v[i];
        sum_B += Math.pow(v[i] - mean, 2);
    }
    let W = Math.pow(sum_A, 2) / Math.max(sum_B, 0.0000001); //To avoid division-by-zero
    
    // Calculate mean, stdev and Z
    m = 0.0038915 * Math.pow(Math.log(n), 3) - 0.083751 * Math.pow(Math.log(n), 2) - 0.31082 * Math.log(n) - 1.5861;
    let pwr = 0.0030302 * Math.pow(Math.log(n), 2) - 0.082676 * Math.log(n) - 0.4803;
    let std = Math.pow(Math.E, pwr);
    let Z = (Math.log(1.0 - W) - m) / std;
    let Zc = jStat.normal.inv(alpha, 0, 1) * -1;
    
    // Calculate P
    let P = 1 - jStat.normal.cdf(Z, 0, 1);
    
    // Skew and kurtosis
    let sum2 = 0;
    let sum3 = 0;
    let sum4 = 0;
    for (let i = 0; i < n; i++) {
        sum2 += Math.pow(v[i] - mean, 2);
        sum3 += Math.pow(v[i] - mean, 3);
        sum4 += Math.pow(v[i] - mean, 4);
    }
    
    // Calculate skew
    let S = ( n * sum3 ) / ( (n-1)*(n-2)*Math.pow(stdev, 3) );
    // Calculate kurtosis
    let K = ( n*(n+1)*sum4 ) / ( (n-1)*(n-2)*(n-3)*Math.pow(stdev, 4) ) - ( 3*Math.pow(n-1,2) ) / ( (n-2)*(n-3) );
    
    // Result
    return [P, Z, Zc, S, K, W];
}

/**
  Performs Bartlett's test for equal variances.
  See: 
        https://www.statisticshowto.datasciencecentral.com/bartletts-test/
  
  Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
*/
function bartlett(arrs, alpha=0.05) {
    // Check if input params is valid
    if (!Array.isArray(arrs)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    // Number of samples
    let k = arrs.length;
    let totalN = 0;
    for (let i = 0; i < k; i++) {
        // Check each sample array
        if (!Array.isArray(arrs[i])) {
            throw "Array " + i + " is not a valid array";
        }
        totalN += arrs[i].length;
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    if (k < 3) {
        throw "Requires three or more samples";
    }
    
    // Calculate totals
    let sumA = 0;
    let sumB = 0;
    let sumC = 0;
    for (let i = 0; i < k; i++) {
        let v1 = arrs[i];
        let DF = v1.length - 1;
        sumA += DF * Math.log(jStat.variance(v1, true));
        sumB += 1 / DF - 1 / (totalN - k);
        sumC += (v1.length - 1) * jStat.variance(v1, true);
    }
    
    // Calculate Chi-square
    let sp2 = 1 / (totalN - k) * sumC;
    let num = (totalN - k) * Math.log(sp2) - sumA;
    let den = 1 + 1 / (3 * (k - 1)) * sumB;
    let chi2 = num / den;
    
    // Find critical Chi-square
    let DF = k - 1;
    let chi2c = jStat.chisquare.inv(1 - alpha, DF);
    
    // P-value
    let P = 1 - jStat.chisquare.cdf(chi2, DF);
    
    return [P, chi2, chi2c];
}

/**
    Performs a Wilcoxon test.
    
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
        type: 1 for Wilcoxon Rank-Sum (independent), 2 for Wilcoxon Signed-Ranks (dependent)
*/
function wilcoxon(v1, v2, sides=2, tail=1, alpha=0.05, type=1) {
    // Check if input params are valid
    if (!Array.isArray(v1) || !Array.isArray(v2)) {
        throw "Parameter is not a valid array";
    }
    if (sides < 1 || sides > 2) {
        throw "Sides must be 1 or 2";
    }
    if (type < 1 || type > 2) {
        throw "Type must be 1 or 2";
    }
    
    if (type == 1) {
        // Independent
        return wilcoxon_ranksum(v1, v2, sides, tail, alpha);
    }
    else if (type == 2) {
        if (v1.length != v2.length) {
            throw "Sample sizes must be equal";
        }
        // Dependent
        return wilcoxon_signedranks(v1, v2, sides, tail, alpha);
    }
}

/**
    Sorts a multi-dimensional array by the first column.
*/
function sortByCol(a, b) {
    a = a[0]
    b = b[0]
    return (a === b) ? 0 : (a < b) ? -1 : 1
}

/**
    Creates the ranks table used for independendent tests.
*/
function createRanksTable(arrs) {
    // Check if input params is valid
    if (!Array.isArray(arrs)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    // Number of samples
    let k = arrs.length;
    for (let i = 0; i < k; i++) {
        // Check each sample array
        if (!Array.isArray(arrs[i])) {
            throw "Array " + i + " is not a valid array";
        }
    }
    
    // Fill ranks table with all samples
    let ranks = [];
    for (let s = 0; s < k; s++) {
        let v = arrs[s];
        for (let i = 0; i < v.length; i++) {
            ranks.push( [v[i], 0, s] );
        }
    }
    
    // Sort by value
    ranks.sort(sortByCol);
    // Set ranks
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        re[1] = i+1;
    }
    // Fix rank duplicates
    for (let s = 0; s < ranks.length - 1; s++) {
        // Get first entry
        let re1 = ranks[s];
        let val1 = re1[0];
        let r_sum = re1[1];
        let cnt = 1;
        // Check if following entries have same value
        for (let e = s + 1; e < ranks.length; e++) {
            let re2 = ranks[e];
            let val2 = re2[0];
            if (val1 == val2) {
                r_sum += re2[1];
                cnt += 1;
            }
        }
        // If several entries have same value, update ranks to mean rank
        if (cnt > 1) {
            let r_mean = r_sum / cnt;
            for (let i = s; i < s + cnt; i++) {
                let re = ranks[i];
                re[1] = r_mean;
            }
        }
    }
    
    return ranks;
}

/**
    Creates the ranks table used for dependendent tests.
*/
function createRanksTableDep(v1, v2) {
    // Check if input params is valid
    if (!Array.isArray(v1)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    if (!Array.isArray(v2)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    
    // Fill ranks table with all samples
    let ranks = [];
    for (let i = 0; i < v1.length; i++) {
        let val1 = v1[i];
        let val2 = v2[i];
        let D = val1 - val2;
        let absD = Math.abs(val1 - val2);
        
        ranks.push( [absD, D, val1, val2, 0, 0] );
    }
    
    // Sort by value
    ranks.sort(sortByCol);
    
    // Set ranks
    let r_cnt = 1;
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        if (re[0] > 0) {
            re[4] = r_cnt;
            r_cnt += 1;
        }
    }
    // Fix rank duplicates
    for (let s = 0; s < ranks.length - 1; s++) {
        // Get first entry
        let re1 = ranks[s];
        let absD1 = re1[0];
        let r_sum = re1[4];
        let cnt = 1;
        // Check if following entries have same value
        for (let e = s + 1; e < ranks.length; e++) {
            let re2 = ranks[e];
            let absD2 = re2[0];
            if (absD1 == absD2) {
                r_sum += re2[4];
                cnt += 1;
            }
        }
        // If several entries have same value, update ranks to mean rank
        if (cnt > 1) {
            let r_mean = r_sum / cnt;
            for (let i = s; i < s + cnt; i++) {
                let re = ranks[i];
                re[4] = r_mean;
            }
        }
    }
    // Set signed ranks
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        re[5] = re[4];
        if (re[1] < 0) {
            re[5] = -re[4];
        }
    }
    
    return ranks;
}

/**
    Performs a Wilcoxon Rank-Sum (independent) test.
    
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
*/
function wilcoxon_ranksum(v1, v2, sides=2, tail=1, alpha=0.05) {
    // Create ranks table
    let ranks = createRanksTable([v1, v2]);
    
    // Check which sample has fewest entries
    let n_min = v1.length;
    let n_other = v2.length;
    let s = 0;
    if (v2.length < v1.length) {
        n_min = v2.length;
        n_other = v1.length;
        s = 1;
    }
    // Calculate rank sum for that sample
    let R = 0;
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        if (re[2] == s) {
            R += re[1];
        }
    }
    
    // Calculate approximate P-value using normal approximation (only exact for n1 and n2 >= 10)
    let Wexp = n_min * (n_min + n_other + 1) / 2;
    let Wstdev = Math.sqrt(n_min * n_other * (n_min + n_other + 1) / 12);
    let Z = (R - Wexp) / Wstdev;
    // Find critical Z
    let Zc = 0;
    if (sides == 2) {
        Zc = jStat.normal.inv(1 - alpha / 2, 0, 1);
    }
    if (sides == 1 && tail == 1) {
        Zc = jStat.normal.inv(alpha, 0, 1);
    }
    if (sides == 1 && tail == 2) {
        Zc = jStat.normal.inv(1 - alpha, 0, 1);
    }
    
    // Get P-value
    let P = 0;
    if (sides == 2) {
        P = jStat.normal.cdf(Z, 0, 1) * sides;
        if (Z > 0) {
            P = (1 - jStat.normal.cdf(Z, 0, 1)) * sides;
        }
    }
    if (sides == 1 && tail == 1) {
        P = jStat.normal.cdf(Z, 0, 1) * sides;
    }
    if (sides == 1 && tail == 2) {
        P = (1 - jStat.normal.cdf(Z, 0, 1)) * sides;
    }
    
    // For smaller sample sizes, use Wilcoxon tables
    if (v1.length <= 10 && v2.length <= 10) {
        // Find critical R
        let i1 = n_min - 3;
        let i2 = n_other - 3;
        let Rc = 0;
        if (sides == 2) {
            Rc = Rc005_2t[i1][i2];
        }
        else {
            Rc = Rc005_1t[i1][i2];
        }
        return [P, R, Rc[0], Rc[1], true];
    }
    else {
        return [P, Z, -Zc, Zc, false];
    }
}

/**
    Performs a Wilcoxon Signed-Ranks (dependent) test.
    
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
*/
function wilcoxon_signedranks(v1, v2, sides=2, tail=1, alpha=0.05) {
    // Create ranks table
    let ranks = createRanksTableDep(v1, v2);
    
    // Calculate sum of ranks and actual sample size
    let SR_pos = 0;
    let SR_neg = 0;
    let n = 0;
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        if (re[5] > 0) {
            SR_pos += re[5];
        }
        else {
            SR_neg += re[5];
        }
        // Actual sample size
        if (re[4] > 0) {
            n += 1;
        }
    }
    
    // Calculate W-score
    let W = 0;
    if (sides == 2) {
        W = Math.min(SR_pos, Math.abs(SR_neg));
    }
    if (sides == 1 && tail == 1) {
        W = SR_pos;
    }
    if  (sides == 1 && tail == 2) {
        W = Math.abs(SR_neg);
    }
    
    // Calculate approximate P-value using normal approximation (only exact for n > 25 )
    let r_mean = n * (n + 1) / 4;
    let r_SD = Math.sqrt(n * (n + 1) * (2 * n + 1) / 24);
    let Z = (W - r_mean) / r_SD;
    
    // Find critical Z
    let Zc = jStat.normal.inv(1 - alpha / sides, 0, 1) * -1;
    
    // Get P-value
    let P = jStat.normal.cdf(Z, 0, 1) * sides;
    
    // For smaller sample sizes, use Wilcoxon tables
    if (n >= 6 && n <= 25) {
        // Find critical W
        let Wc = 0;
        if (sides == 2) {
            Wc = Wc005_2t[n - 5];
        }
        else {
            Wc = Wc005_1t[n - 5];
        }
        return [P, W, Wc, true];
    }
    else {
        return [P, Z, Zc, false];
    }
}

/**
    Wilcoxon Signed-Ranks single sample for the mean.
  
    Params:
        v1: Sample
        mean: Mean to compare with
        sides: One- or two-tailed test
        tail: Left (1) or right (2) for one-tailed test
        alpha: Significance level
*/
function wilcoxon_signedranks_single(v1, mean, sides=2, tail=1, alpha=0.05) {
    // Fill ranks table with all samples
    let ranks = [];
    for (let i = 0; i < v1.length; i++) {
        let val1 = v1[i];
        let D = val1 - mean;
        let absD = Math.abs(val1 - mean);
        
        ranks.push( [absD, D, val1, mean, 0, 0] );
    }
    
    // Sort by value
    ranks.sort(sortByCol);
    
    // Set ranks
    let r_cnt = 1;
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        if (re[0] > 0) {
            re[4] = r_cnt;
            r_cnt += 1;
        }
    }
    // Fix rank duplicates
    for (let s = 0; s < ranks.length - 1; s++) {
        // Get first entry
        let re1 = ranks[s];
        let absD1 = re1[0];
        let r_sum = re1[4];
        let cnt = 1;
        // Check if following entries have same value
        for (let e = s + 1; e < ranks.length; e++) {
            let re2 = ranks[e];
            let absD2 = re2[0];
            if (absD1 == absD2) {
                r_sum += re2[4];
                cnt += 1;
            }
        }
        // If several entries have same value, update ranks to mean rank
        if (cnt > 1) {
            let r_mean = r_sum / cnt;
            for (let i = s; i < s + cnt; i++) {
                let re = ranks[i];
                re[4] = r_mean;
            }
        }
    }
    // Set signed ranks
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        re[5] = re[4];
        if (re[1] < 0) {
            re[5] = -re[4];
        }
    }
    
    // Calculate sum of ranks and actual sample size
    let SR_pos = 0;
    let SR_neg = 0;
    let n = 0;
    for (let i = 0; i < ranks.length; i++) {
        let re = ranks[i];
        if (re[5] > 0) {
            SR_pos += re[5];
        }
        else {
            SR_neg += re[5];
        }
        // Actual sample size
        if (re[4] > 0) {
            n += 1;
        }
    }
    
    // Calculate W-score
    let W = 0;
    if (sides == 2) {
        W = Math.min(SR_pos, Math.abs(SR_neg));
    }
    if (sides == 1 && tail == 1) {
        W = SR_pos;
    }
    if  (sides == 1 && tail == 2) {
        W = Math.abs(SR_neg);
    }
    
    // Calculate approximate P-value using normal approximation (only exact for n > 25 )
    let r_mean = n * (n + 1) / 4;
    let r_SD = Math.sqrt(n * (n + 1) * (2 * n + 1) / 24);
    let Z = (W - r_mean) / r_SD;
    
    // Find critical Z
    let Zc = jStat.normal.inv(1 - alpha / sides, 0, 1);
    
    // Get P-value
    let P = jStat.normal.cdf(Z, 0, 1) * sides;
    
    // For smaller sample sizes, use Wilcoxon tables
    if (n >= 6 && n <= 25) {
        // Find critical W
        let Wc = 0;
        if (sides == 2) {
            Wc = Wc005_2t[n - 5];
        }
        else {
            Wc = Wc005_1t[n - 5];
        }
        return [P, W, Wc, true];
    }
    else {
        return [P, Z, Zc, false];
    }
}


/**
    Performs a Kruskal-Wallis test.
  
  Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
*/
function kruskalwallis(arrs, alpha=0.05) {
    // Create ranks table
    let ranks = createRanksTable(arrs);
    
    // Number of samples
    let k = arrs.length;
    let totalN = ranks.length;
    
    // Calculate R-scores (sum of ranks for each sample)
    let R = new Array(k).fill(0);
    for (let i = 0; i < ranks.length; i++) {
        let ri = ranks[i];
        // Sample index
        let si = ri[2];
        // Rank
        R[si] += ri[1];
    }
    
    // Calculate H-score
    let R_sum = 0;
    for (let i = 0; i < k; i++) {
        // Sum of squared rank sums divided by n
        let ns = arrs[i].length;
        R_sum += Math.pow(R[i], 2) / ns;
    }
    let H = 12 / (totalN * (totalN + 1)) * R_sum - 3 * (totalN + 1);
    
    // Find critical H-value
    let DF = k - 1;
    let Hc = jStat.chisquare.inv(1 - alpha, DF);
    // P-value
    let P = 1 - jStat.chisquare.cdf(H, DF);
    
    return [P, H, Hc];
}

/**
    Performs Dunn's post-test.
    
    Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
        correction: Significance level correction (none or Bonferroni)
*/
function dunns_posttest(arrs, alpha=0.05, correction="none") {
    // Check if input params is valid
    if (!Array.isArray(arrs)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    // Number of samples
    let k = arrs.length;
    for (let i = 0; i < k; i++) {
        // Check each sample array
        if (!Array.isArray(arrs[i])) {
            throw "Array " + i + " is not a valid array";
        }
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    if (k < 3) {
        throw "Requires three or more samples";
    }
    
    // Correction
    if (correction == "bonferroni") {
        let m = k * (k - 1) / 2;
        alpha = alpha / m;
    }
    
    // Find critical Z-value
    let Zc = jStat.normal.inv(1 - alpha / 2, 0, 1);
    
    // Create ranks table
    let ranks = createRanksTable(arrs);
    let totalN = ranks.length;
    
    // Calculate sum of ranks and number of ties
    let r_means = new Array(k).fill(0);
    let ties = new Array(k).fill(0);
    let sum_ties = 0;
    for (let i = 0; i < totalN; i++) {
        let ri = ranks[i];
        // Sample index
        let si = ri[2];
        // Rank
        r_means[si] += ri[1];
        // Check if ties (rank is a decimal number)
        if (ri[1] % 1 != 0) {
            ties[si] += 1;
            sum_ties += 1;
        }
    }
    // Mean ranks
    for (let i = 0; i < k; i++) {
        let n = arrs[i].length;
        r_means[i] /= n;
    }
    
    // Pair-wise comparison
    let P_vals = [];
    for (let i = 0; i < k; i++) {
        for (let j = i + 1; j < k; j++) {
            // Samples
            let v1 = arrs[i];
            let v2 = arrs[j];
            
            // Absolute difference between mean ranks
            let absD = Math.abs(r_means[i] - r_means[j]);
            
            // Ties correction
            let sum = 0;
            if (sum_ties > 0) {
                // Calculate sum of ties
                for (let ti = 0; ti < ties.length; ti++) {
                    sum += Math.pow(ties[ti],3) - ties[ti];
                }
                sum /= (totalN - 1);
            }
            
            // Divider
            let div = Math.sqrt(((totalN * (totalN + 1) - sum) / 12) * (1 / v1.length + 1 / v2.length));
            
            // Calculate Z-score
            let Z = absD / div;
            
            // Get P-value
            let P = (1 - jStat.normal.cdf(Z, 0, 1)) * 2;
            
            // Result
            let res = [ i + 1, j + 1, P, alpha, Z, Zc, "Z"];
            P_vals.push(res);
        }
    }
    
    return P_vals;
}

/**
    Performs Wilcoxon pair-wise comparison post-test.
    
    Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
        correction: Significance level correction (none or Bonferroni)
*/
function wilcoxon_posttest(arrs, alpha=0.05, correction="none") {
    // Check if input params is valid
    if (!Array.isArray(arrs)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    // Number of samples
    let k = arrs.length;
    for (let i = 0; i < k; i++) {
        // Check each sample array
        if (!Array.isArray(arrs[i])) {
            throw "Array " + i + " is not a valid array";
        }
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    if (k < 3) {
        throw "Requires three or more samples";
    }
    
    // Correction
    if (correction == "bonferroni") {
        let m = k * (k - 1) / 2;
        alpha = alpha / m;
    }
    
    // Pair-wise comparison
    let P_vals = [];
    for (let i1 = 0; i1 < k; i1++) {
        for (let i2 = i1 + 1; i2 < k; i2++) {
            // Samples
            let v1 = arrs[i1];
            let v2 = arrs[i2];
            
            // Create ranks table
            let ranks = createRanksTable([v1, v2]);

            // Check which sample has fewest entries
            let n_min = v1.length;
            let n_other = v2.length;
            let s = 0;
            if (v2.length < v1.length) {
                n_min = v2.length;
                n_other = v1.length;
                s = 1;
            }
            
            // Calculate rank sum for that sample
            let R = 0;
            for (let i = 0; i < ranks.length; i++) {
                let re = ranks[i];
                if (re[2] == s) {
                    R += re[1];
                }
            }

            let Wexp = n_min * (n_min + n_other + 1) / 2;
            let Wstdev = Math.sqrt(n_min * n_other * (n_min + n_other + 1) / 12);
            let Z = (R - Wexp) / Wstdev;
            let Zc = jStat.normal.inv(1 - alpha / 2, 0, 1);
            // Get P-value
            let P = 0;
            if (Z <= 0) {
                P = jStat.normal.cdf(Z, 0, 1) * 2;
            }
            else {
                P = (1 - jStat.normal.cdf(Z, 0, 1)) * 2;
            }

            // Result
            let res = [ i1 + 1, i2 + 1, P, alpha, Z, Zc, "Z"];
            P_vals.push(res);
        }
    }
    return P_vals;
}

/**
    Performs a Friedman test.
  
  Params:
        arr: Array of sample arrays (3 or more samples)
        alpha: Significance level
*/
function friedman(arrs, alpha=0.05) {
    // Create ranks table
    let k = arrs.length;
    let n = arrs[0].length;
    
    // Fill ranks table with all samples
    let ranks = [];
    for (let i = 0; i < n; i++) {
        let row = [];
        for (let s = 0; s < k; s++) {
            let v = arrs[s];
            row.push( [v[i], 0, s] );
        }
        
        row.sort(sortByCol);
        
        ranks.push( row );
    }
    
    //Fill in ranks
    for (let i = 0; i < n; i++)
    {
        let row = ranks[i];
        for (let s = 0; s < k; s++)
        {
            let e = row[s];
            e[1] = s + 1;
        }
        
        //Fix rank duplicates
        for (let s = 0; s < k - 1; s++)
        {
            let e = row[s];
            let cVal = e[0];
            
            //Check next ranks
            let rankSum = e[1];
            let cnt = 1;
            //Calculate mean ranks
            for (let p = s + 1; p < k; p++)
            {
                let en = row[p];
                let nVal = en[0];
                if (nVal == cVal)
                {
                    rankSum += en[1];
                    cnt++;
                }
            }
            let meanRank = rankSum / cnt;
            //Update ranks for duplicates
            if (cnt > 1)
            {
                for (let p = s; p < k; p++)
                {
                    let en = row[p];
                    let nVal = en[0];
                    if (nVal == cVal)
                    {
                        en[1] = meanRank;
                    }
                }
            }
        }
    }
    //Rankings table done!
    
    // Number of samples
    let totalN = ranks.length;
    
    // Calculate R-scores
    let R = new Array(k).fill(0);
    for (let i = 0; i < n; i++)
    {
        let row = ranks[i];
        for (let s = 0; s < k; s++)
        {
            let e = row[s];
            // Sum ranks for each sample
            R[e[2]] += e[1];
        }
    }
    
    let R_sum = 0;
    for (let s = 0; s < k; s++)
    {
        R_sum += Math.pow(R[s], 2);
    }
    
    //Calculate Q-score
    let Q = 12 / (n * k * (k + 1)) * R_sum - 3 * n * (k + 1);
    
    // Find critical Q-value
    let DF = k - 1;
    let Qc = jStat.chisquare.inv(1 - alpha, DF);
    // P-value
    let P = 1 - jStat.chisquare.cdf(Q, DF);
    
    return [P, Q, Qc];
}

function wilcoxon_sr_posttest(arrs, alpha=0.05, correction="none") {
    // Check if input params is valid
    if (!Array.isArray(arrs)) {
        // Array of sample arrays
        throw "Parameter is not a valid array";
    }
    // Number of samples
    let k = arrs.length;
    for (let i = 0; i < k; i++) {
        // Check each sample array
        if (!Array.isArray(arrs[i])) {
            throw "Array " + i + " is not a valid array";
        }
    }
    if (alpha < 0 || alpha > 1) {
        throw "Alpha must be between 0 and 1";
    }
    if (k < 3) {
        throw "Requires three or more samples";
    }
    
    // Correction
    if (correction == "bonferroni") {
        let m = k * (k - 1) / 2;
        alpha = alpha / m;
    }
    
    // Pair-wise comparison
    let P_vals = [];
    for (let i1 = 0; i1 < k; i1++) {
        for (let i2 = i1 + 1; i2 < k; i2++) {
            // Samples
            let v1 = arrs[i1];
            let v2 = arrs[i2];
            
            // Create ranks table
            let ranks = createRanksTableDep(v1, v2);
            
            // Calculate sum of ranks and actual sample size
            let SR_pos = 0;
            let SR_neg = 0;
            let n = 0;
            for (let i = 0; i < ranks.length; i++) {
                let re = ranks[i];
                if (re[5] > 0) {
                    SR_pos += re[5];
                }
                else {
                    SR_neg += re[5];
                }
                // Actual sample size
                if (re[4] > 0) {
                    n += 1;
                }
            }
            
            // Calculate W-score
            let W = Math.min(SR_pos, Math.abs(SR_neg));
            
            // Calculate approximate P-value using normal approximation (only exact for n > 25 )
            let r_mean = n * (n + 1) / 4;
            let r_SD = Math.sqrt(n * (n + 1) * (2 * n + 1) / 24);
            let Z = (W - r_mean) / r_SD;
            
            // Find critical Z
            let Zc = jStat.normal.inv(1 - alpha / 2, 0, 1) * -1;
            
            // Get P-value
            let P = jStat.normal.cdf(Z, 0, 1) * 2;

            // Result
            let res = [ i1 + 1, i2 + 1, P, alpha, Z, Zc, "Z"];
            P_vals.push(res);
        }
    }
    return P_vals;
}


/**
    Calculates R correlation coefficient between two samples.
  
    Params:
        v1: Sample 1
        v2: Sample 2
        alpha: Significance level
*/
function correlation_R(v1, v2, alpha=0.05) {
    // Sample means
    let mean1 = jStat.mean(v1);
    let mean2 = jStat.mean(v2);
    // Sample sizes
    let n1 = v1.length;
    let n2 = v2.length;
        
    if (n1 != n2) {
            throw("Sample sizes must be equal");
    }
    let n = n1;
    
    // Calculate R coefficient
    let sumX = 0;
    let sumY = 0;
    let sumXY = 0;
    let sumX2 = 0;
    let sumY2 = 0;
    for (let i = 0; i < n; i++) {
        let x = v1[i];
        let y = v2[i];
        
        sumX += x;
        sumY += y;
        sumXY += x * y;
        sumX2 += Math.pow(x, 2);
        sumY2 += Math.pow(y, 2);
    }
    let r = (n*sumXY - sumX*sumY) / Math.sqrt( (n*sumX2 - Math.pow(sumX, 2)) * (n*sumY2 - Math.pow(sumY, 2)) );
    
    // Test for significance
    let T = r / Math.sqrt( (1-Math.pow(r, 2))/(n - 2) );
    // Find critical T-value
    let DF = n - 2;
    let Tc = jStat.studentt.inv(1 - alpha / 2, DF);
    // Get P-value
    let P = jStat.studentt.cdf(T, DF) * 2;
    if (T > 0) {
        P = (1 - jStat.studentt.cdf(T, DF)) * 2;
    }
    
    return [r, P, T, Tc];
}

/**
    Calculates linear regression between two samples.
  
    Params:
        v1: Sample 1
        v2: Sample 2
        alpha: Significance level
*/
function linear_regression(v1, v2, alpha=0.05) {
    // Sample means
    let mean1 = jStat.mean(v1);
    let mean2 = jStat.mean(v2);
    // Sample sizes
    let n1 = v1.length;
    let n2 = v2.length;
        
    if (n1 != n2) {
        throw("Sample sizes must be equal");
    }
    let n = n1;
    
    // Calculate R coefficient
    let sumX = 0;
    let sumY = 0;
    let sumXY = 0;
    let sumX2 = 0;
    let sumY2 = 0;
    for (let i = 0; i < n; i++) {
        let x = v1[i];
        let y = v2[i];
        
        sumX += x;
        sumY += y;
        sumXY += x * y;
        sumX2 += Math.pow(x, 2);
        sumY2 += Math.pow(y, 2);
    }
    
    // Calculate line of best fit
    let b = (n * sumXY - sumX * sumY) / (n * sumX2 - Math.pow(sumX, 2) );
    let a = sumY / n - b * sumX / n;
    
    // Calculate coefficient of determination
    let SST = sumY2 - Math.pow(sumY, 2) / n;
    let SSE = sumY2 - a * sumY - b * sumXY;
    let SSR = SST - SSE;
    let R2 = SSR / SST;
    
    // Calculate F-score
    let F = SSR / (SSE / (n - 2));
        
    // Find critical F
    let DF1 = 1;
    let DF2 = n - 2;
    
    // Calculate critical F-value
    let Fc = jStat.centralF.inv(1 - alpha, DF1, DF2);
    // P-value
    let P = jStat.centralF.cdf(F, DF1, DF2);
    if (F > 0) {
        P = (1 - jStat.centralF.cdf(F, DF1, DF2));
    }
    
    return [a, b, R2, P, F, Fc];
}


/**
    Generalized Extreme Studentized (ESD) test for potential outliers.
  
    Params:
        v1: Sample
        alpha: Significance level
*/
function outliers_esd(v1, alpha=0.05) {
    // Check if input params are valid
    if (!Array.isArray(v1)) {
        throw "Parameter is not a valid array";
    }
    
    // Find number of trials
    let trials = find_trials(v1, alpha);
    
    // Create trials list
    let tlist = [];
    let vc = v1.slice(); //Copy the array
    for (let t = 0; t < trials; t++) {
        // Find outlier
        let index = find_outlier(vc);
        let val = vc[index];
        // Check if outlier is significant
        let sig = check_significance(vc, alpha);
        // Add to list
        tlist.push( [val, index, sig[0], sig[1], sig[2]] );
        // Remove outlier
        vc.splice(index, 1);
    } 
    
    // Iterate over the trials list to remove none-outliers
    // (remove the non-significant entries from the tail of the list)
    let stop = false;
    while (!stop) {
        if (tlist.length == 0) {
            stop = true;
        }
        else {
            let te = tlist[tlist.length-1];
            if (te[2] == false) {
                tlist.splice(tlist.length-1, 1);
            }
            else {
                stop = true;
            }
        }
    }
    
    // Iterate over all values to check for outliers
    let res = [];
    for (let i = 0; i < v1.length; i++) {
        let val = v1[i];
        
        // Check if value is outlier
        let outl = false;
        for (let t = 0; t < tlist.length; t++) {
            let te = tlist[t];
            if (te[0] == val) {
                outl = true;
            }
        }
        
        res.push( [val, outl] );
    }
    
    // Return result
    return [tlist, res];
}

/**
    Finds number of trials for ESD test.
*/
function find_trials(v1, alpha) {
    // Make a copy of the array
    let vc = v1.slice();
    
    let stop = false;
    let no = 0;
    
    while (!stop) {
        // Increase trials
        no += 1;
        
        // Check size of array
        if (vc.length <= 4) {
            stop = true;
        }
        
        // Check if normally distributed
        let res = shapiro_wilk(vc, alpha);
        if (res[0] > alpha) {
            stop = true;
        }
        
        // Find and remove outlier
        let index = find_outlier(vc);
        vc.splice(index, 1);
    }
    
    // Run at least 2 trials
    if (no < 2) {
        no = 2;
    }
    
    return no;
}

/**
    Checks if removing an outlier is significant or not.
*/
function check_significance(vc, alpha) {
    // Calculate G-value
    let min = jStat.min(vc);
    let max = jStat.max(vc);
    let mean = jStat.mean(vc);
    let std = jStat.stdev(vc, true);
    let mean_min = mean - min;
    let max_mean = max - mean;
    let n = vc.length;
    let G = Math.max(mean_min, max_mean) / std;
    
    // Find critical T-value
    let c_alpha = alpha / n;
    let DF = n - 2;
    let Tc = jStat.studentt.inv(1 - c_alpha / 2, DF);
    // Find critical G-value
    let Gc = (n - 1) * Tc / Math.sqrt(n * (DF + Math.pow(Tc,2)));
    
    // Return result
    return [G > Gc, G, Gc];
    
}

/**
    Finds a potential outlier (the value with the largest absolute difference to
    the mean) in a sample.
*/
function find_outlier(vc) {
    // Remove outlier
    let maxD = 0;
    let index = 0;
    let mean = jStat.mean(vc);
    for (let i = 0; i < vc.length; i++) {
        let currentD = Math.abs(vc[i] - mean);
        if (currentD > maxD) {
            maxD = currentD;
            index = i;
        }
    }
    return index;
}
/*****************************************************

Functions for calculating statistical power.

******************************************************/

/**
    Power calculation for three or more samples.
  
    Params:
        arrs: The samples
        alpha: Significance level
        correction: Significance level correction (none or Bonferroni)
        dep: True for paired samples, false for independent samples
        parametric: true for parametric tests, false otherwise
*/
function power_anova(arrs, alpha=0.05, correction="none", dep=false, parametric=false) {   
    // Number of samples
    let no_k = arrs.length;
    
    // Correction
    if (correction == "bonferroni") {
        let m = no_k * (no_k - 1) / 2;
        alpha = alpha / m;
    }
    
    // Find pair with smallest difference in means
    // See discussion here:
    // http://www.3rs-reduction.co.uk/html/6__power_and_sample_size.html
    let i1 = 0;
    let i2 = 0;
    let smallD = 100000;
    for (let i = 0; i < no_k; i++) {
        for (let j = i + 1; j < no_k; j++) {
            // Samples
            let v1 = arrs[i];
            let v2 = arrs[j];
            let mean1 = jStat.mean(v1);
            let mean2 = jStat.mean(v2);
            
            // Calculate difference in means
            let absD = Math.abs(mean1 - mean2);
            if (absD < smallD) {
                let smallD = absD;
                i1 = i;
                i2 = j;
            }
        }
    }
    
    // Calculate power and sample size for the pair with smallest difference in means
    let v1 = arrs[i1];
    let v2 = arrs[i2];
    let res = new Array(8).fill(0);
    
    // Error check - sample means must be different
    let mean1 = jStat.mean(v1);
    let mean2 = jStat.mean(v2);
    if (mean1 == mean2) {
        hide("power");
        return res;
    }
    show("power");
    
    res[0] = calc_min_samplesize_2s(v1, v2, 2, alpha, 0.2, dep, parametric);
    res[1] = calc_min_samplesize_2s(v1, v2, 2, alpha, 0.5, dep, parametric);
    res[2] = calc_min_samplesize_2s(v1, v2, 2, alpha, 0.8, dep, parametric);
    res[3] = calc_power_2s(v1, v2, 2, alpha, dep, parametric);
    
    // Round up to nearest integer
    for (let i = 0; i < 3; i++) {
        res[i] = Math.ceil(res[i]);
    }
    
    return res;
}

/**
    Power calculation for two samples.
  
    Params:
        v1: Sample 1
        v2: Sample 2
        alpha: Significance level
        sides: One- or two-tailed test
        dep: True for paired samples, false for independent samples
        parametric: true for parametric tests, false otherwise
*/
function power_2s(v1, v2, alpha=0.05, sides=2, dep=false, parametric=true) {
    let res = new Array(8).fill(0);
    
    // Error check - sample means must be different
    let mean1 = jStat.mean(v1);
    let mean2 = jStat.mean(v2);
    if (mean1 == mean2) {
        hide("power");
        return res;
    }
    show("power");
    
    res[0] = calc_min_samplesize_2s(v1, v2, sides, alpha, 0.2, dep, parametric);
    res[1] = calc_min_samplesize_2s(v1, v2, sides, alpha, 0.5, dep, parametric);
    res[2] = calc_min_samplesize_2s(v1, v2, sides, alpha, 0.8, dep, parametric);
    res[3] = calc_power_2s(v1, v2, sides, alpha, dep, parametric);
    
    // Round up to nearest integer
    for (let i = 0; i < 3; i++) {
        res[i] = Math.ceil(res[i]);
    }
    
    return res;
}

/**
    Power calculation for two samples.
  
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One or two-tailed
        alpha: Significance level
        dep: True for paired samples, false for independent samples
        parametric: true for parametric tests, false otherwise
*/
function calc_power_2s(v1, v2, sides, alpha, dep, parametric) {
    // Sample sizes
    let n1 = v1.length;
    let n2 = v2.length;
    
    let minD = 100000;
    let pwr = 0;
    
    for (let p = 0; p <= 1000; p++) {
        let cpwr = p * 0.001;
        let n = calc_min_samplesize_2s(v1, v2, sides, alpha, cpwr, dep, parametric);
        
        let cD = Math.abs(n * 2 - (n1 + n2));
        if (cD < minD) {
            minD = cD;
            pwr = cpwr;
        }
    }
    
    return pwr * 100;
}

/**
    Calculates min sample size for two samples.
  
    Params:
        v1: Sample 1
        v2: Sample 2
        sides: One or two-tailed
        alpha: Significance level
        power: Power to calculate sample size for
        dep: True for paired samples, false for independent samples
        parametric: true for parametric tests, false otherwise
*/
function calc_min_samplesize_2s(v1, v2, sides, alpha, power, dep, parametric) {
    // Sample means
    let mean1 = jStat.mean(v1);
    let mean2 = jStat.mean(v2);
    // Sample sizes
    let n1 = v1.length;
    let n2 = v2.length;
    // Standard deviations for the samples
    let std1 = jStat.stdev(v1, true);
    let std2 = jStat.stdev(v2, true);
    
    // Find Z-scores
    let Za = jStat.normal.inv(1 - alpha / sides, 0, 1);
    let Zb = jStat.normal.inv(1 - (1 - power), 0, 1);
    
    // Standard deviation
    let std = 0;
    if (dep == false) {
        // Calculate pooled stdev for 2 independent samples
        std = Math.sqrt( (Math.pow(std1, 2) + Math.pow(std2, 2)) / 2 );
    }
    else {
        // Calculate within pairs stdev for 2 paired samples
        let sumD = 0;
        let sumD2 = 0;
        let s = v1.length;
        for (let i = 0; i < s; i++) {
            sumD += v2[i] - v1[i];
            sumD2 += Math.pow(v2[i] - v1[i], 2);
        }
        std = Math.sqrt( (sumD2 - Math.pow(sumD, 2) / s) / (s-1) );
    }
    
    // Calculate n
    let diff_sq = Math.pow(mean1 - mean2, 2);
    let n = 0;
    if (dep == false) {
        let r = Math.max(n1, n2) / Math.min(n1, n2);
        let k = (r + 1) / r;
        n = k * Math.pow(std, 2) * Math.pow(Zb + Za, 2) / diff_sq;
    }
    else {
        n = 2 * Math.pow(std, 2) * Math.pow(Zb + Za, 2) / diff_sq;
    }
    
    // For non-parametric tests, add 15% according to the discussion at:
    // http://www.graphpad.com/guides/prism/6/statistics/index.htm?stat_sample_size_for_nonparametric_.htm
    if (parametric == false) {
        n = n * 1.15;
    }
    
    return n;
}

/**
    Power calculation for one sample and a mean.
  
    Params:
        v1: Sample
        mean: The mean
        alpha: Significance level
        sides: One- or two-tailed test
        parametric: true for parametric tests, false otherwise
*/
function power_1s(v1, mean, alpha=0.05, sides=2, parametric=true) {   
    let res = new Array(8).fill(0);
    
    // Error check - sample means must be different
    let mean1 = jStat.mean(v1);
    if (mean1 == mean) {
        hide("power");
        return res;
    }
    show("power");
    
    res[0] = calc_min_samplesize_1s(v1, mean, sides, alpha, 0.2, parametric);
    res[1] = calc_min_samplesize_1s(v1, mean, sides, alpha, 0.5, parametric);
    res[2] = calc_min_samplesize_1s(v1, mean, sides, alpha, 0.8, parametric);
    res[4] = calc_power_1s(v1, mean, sides, alpha, parametric);
    
    // Round up to nearest integer
    for (let i = 0; i < 3; i++) {
        res[i] = Math.ceil(res[i]);
    }
    
    return res;
}

/**
    Calculates the power for one sample and a mean.
    
    Params:
        v1: Sample
        mean: The mean
        sides: One or two-tailed
        alpha: Significance level
        parametric: true for parametric tests, false otherwise
*/
function calc_power_1s(v1, mean, sides, alpha, parametric) {
    // Sample size
    let n1 = v1.length;
    
    let minD = 100000;
    let pwr = 0;
    
    for (let p = 0; p <= 1000; p++) {
        let cpwr = p * 0.001;
        let n = calc_min_samplesize_1s(v1, mean, sides, alpha, cpwr, parametric);
        
        let cD = Math.abs(n - n1);
        if (cD < minD) {
            minD = cD;
            pwr = cpwr;
        }
    }
    
    return pwr * 100;
}

/**
    Calculates min sample size for one sample and a mean.
  
    Params:
        v1: Sample
        mean: The mean
        sides: One or two-tailed
        alpha: Significance level
        power: Power to calculate sample size for
        parametric: true for parametric tests, false otherwise
*/
function calc_min_samplesize_1s(v1, mean, sides, alpha, power, parametric) {
    // Sample mean
    let mean1 = jStat.mean(v1);
    // Sample size
    let n1 = v1.length;
    // Standard deviations for the sample
    let std1 = jStat.stdev(v1, true);
    
    // Find Z-scores
    let Za = jStat.normal.inv(1 - alpha / sides, 0, 1);
    let Zb = jStat.normal.inv(1 - (1 - power), 0, 1);
    
    // Calculate n
    let diff_sq = Math.pow(mean - mean1, 2);
    let var1 = jStat.variance(v1, true);
    let n = Math.pow(Za + Zb, 2) * var1 / diff_sq;
    
    // For non-parametric tests, add 15% according to the discussion at:
    // http://www.graphpad.com/guides/prism/6/statistics/index.htm?stat_sample_size_for_nonparametric_.htm
    if (parametric == false) {
        n = n * 1.15;
    }
    
    return n;
}
/*****************************************************

Functions for generating and displaying distributions.

******************************************************/

/**
    Returns the integer value from an input field.
*/
function get_int(id) {
    let str = document.getElementById(id).value;
    str = str.trim();
    let val = parseInt(str);
    if (isNaN(val)) {
        throw("Invalid integer value: " + str);
    }
    return val;
}

/**
    Returns the float value from an input field.
*/
function get_float(id) {
    let str = document.getElementById(id).value;
    str = str.trim();
    let val = parseFloat(str);
    if (isNaN(val)) {
        throw("Invalid integer value: " + str);
    }
    return val;
}

/**
    Generates and visualizes the t-distribution.
*/
function dist_t() {
    try {
        let DF = get_int("df");
        if (DF < 5 || DF > 200) {
            throw("DF must be between 5 and 200");
        }
        
        let vals = generate_dist(-4, 4, 0.1, 0, 0, jStat.studentt.pdf, [DF], 0);
        visualize_dist(vals);
        
        // P-values
        let str = "";
        for (let i = 0; i <= 2; i += 0.2) {
            str += "<tr><td class='border'>" + i.toFixed(1) + "</td>";
            let P1 = (1 - jStat.studentt.cdf(i, DF));
            let P2 = P1 * 2;
            str += "<td class='border'>" + P1.toFixed(4) + "</td>";
            str += "<td class='border'>" + P2.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p1").innerHTML = str;
        
        str = "";
        for (let i = 2.2; i <= 4.1; i += 0.2) {
            str += "<tr><td class='border'>" + i.toFixed(1) + "</td>";
            let P1 = 1 - jStat.studentt.cdf(i, DF);
            let P2 = P1 * 2;
            str += "<td class='border'>" + P1.toFixed(4) + "</td>";
            str += "<td class='border'>" + P2.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p2").innerHTML = str;
        
        show("viz");
        clear_error();
    }
    catch (e) {
        hide("viz");
        show_error(e);
    }
}

/**
    Generates and visualizes the normal distribution.
*/
function dist_norm() {
    try {
        let mean = get_float("mean");
        let stdev = get_float("stdev");
        
        let range = Math.ceil(4 * stdev);
        
        let vals = generate_dist(-range, range, 0.1, 0, 0, jStat.normal.pdf, [0, stdev], mean);
        visualize_dist(vals);
        
        // P-values
        let str = "";
        for (let i = 0; i <= range / 2 + 0.1; i += 0.2) {
            let left = mean - i;
            let right = mean + i;
            str += "<tr><td class='border'>" + left.toFixed(1) + "</td>";
            str += "<td class='border'>" + right.toFixed(1) + "</td>";
            let P1 = 1 - jStat.normal.cdf(i, 0, stdev);
            let P2 = P1 * 2;
            str += "<td class='border'>" + P1.toFixed(4) + "</td>";
            str += "<td class='border'>" + P2.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p1").innerHTML = str;
        
        str = "";
        for (let i = range / 2 + 0.2; i <= range + 0.1; i += 0.2) {
            let left = mean - i;
            let right = mean + i;
            str += "<tr><td class='border'>" + left.toFixed(1) + "</td>";
            str += "<td class='border'>" + right.toFixed(1) + "</td>";
            let P1 = 1 - jStat.normal.cdf(i, 0, stdev);
            let P2 = P1 * 2;
            str += "<td class='border'>" + P1.toFixed(4) + "</td>";
            str += "<td class='border'>" + P2.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p2").innerHTML = str;
        
        show("viz");
        clear_error();
    }
    catch (e) {
        hide("viz");
        show_error(e);
    }
}

/**
    Generates and visualizes the F-distribution.
*/
function dist_f() {
    try {
        let DF1 = get_int("df1");
        if (DF1 < 1 || DF1 > 100) {
            throw("DF<sub>numerator</sub> must be between 1 and 100");
        }
        let DF2 = get_int("df2");
        if (DF2 < 5 || DF2 > 100) {
            throw("DF<sub>denominator</sub> must be between 5 and 100");
        }
        
        let vals = generate_dist(0, 50, 0.2, 0.001, 5, jStat.centralF.pdf, [DF1, DF2], 0);
        visualize_dist(vals);
        
        // Find end values
        let x_vals = vals[0];
        let end = x_vals[x_vals.length - 1];
        let sep = Math.ceil(end / 2);
        
        // P-values
        let str = "";
        for (let i = 0; i <= sep + 0.1; i += 0.2) {
            str += "<tr><td class='border'>" + i.toFixed(1) + "</td>";
            let P = 1 - jStat.centralF.cdf(i, DF1, DF2);
            str += "<td class='border'>" + P.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p1").innerHTML = str;
        
        str = "";
        for (let i = sep + 0.2; i <= end + 0.1; i += 0.2) {
            str += "<tr><td class='border'>" + i.toFixed(1) + "</td>";
            let P = 1 - jStat.centralF.cdf(i, DF1, DF2);
            str += "<td class='border'>" + P.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p2").innerHTML = str;
        
        show("viz");
        clear_error();
    }
    catch (e) {
        hide("viz");
        show_error(e);
    }
}

/**
    Generates and visualizes the Chi-Square distribution.
*/
function dist_chi2() {
    try {
        let DF = get_int("df");
        if (DF < 1 || DF > 30) {
            throw("DF must be between 1 and 30");
        }
        
        let vals = generate_dist(0, 100, 0.5, 0.00001, 10, jStat.chisquare.pdf, [DF], 0);
        visualize_dist(vals);
        
        // Find end values
        let x_vals = vals[0];
        let end = x_vals[x_vals.length - 1];
        let sep = Math.ceil(end / 2);
        
        // P-values
        let str = "";
        for (let i = 0; i <= sep + 0.3; i += 0.5) {
            str += "<tr><td class='border'>" + i.toFixed(1) + "</td>";
            let P = 1 - jStat.chisquare.cdf(i, DF);
            str += "<td class='border'>" + P.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p1").innerHTML = str;
        
        str = "";
        for (let i = sep + 0.5; i <= end + 0.3; i += 0.5) {
            str += "<tr><td class='border'>" + i.toFixed(1) + "</td>";
            let P = 1 - jStat.chisquare.cdf(i, DF);
            str += "<td class='border'>" + P.toFixed(4) + "</td></tr>";
        }
        document.getElementById("p2").innerHTML = str;
        
        show("viz");
        clear_error();
    }
    catch (e) {
        hide("viz");
        show_error(e);
    }
}

/**
    Samples values from a distribution.
    
    Params:
        min_x: First x-value to start sampling from
        max_x: Last x-value to sample from
        inc: Step size
        threshold: Sampling stops when the y-value is below this value (0 for not stopping)
        threshold_x: Threshold stopping is only used after this x-value
        distfunc: Distribution funciton to sample from
        params: Array of parameter values
*/
function generate_dist(min_x, max_x, step, threshold, threshold_x, distfunc, params, x_diff) {
    try {
        let y_vals = [];
        let x_vals = [];
        for (let i = min_x; i <= max_x; i += step) {
            x_vals.push( i + x_diff );
            let y_val = 0;
            if (params.length == 1) {
                // Single DF
                y_val = distfunc(i, params[0]).toFixed(4);
            }
            if (params.length == 2) {
                // Two DF
                y_val = distfunc(i, params[0], params[1]).toFixed(4);
            }
            y_vals.push( y_val );
            
            if (y_val <= threshold && i > threshold_x) break;
        }
        return [x_vals, y_vals];
    }
    catch (e) {
        throw(e);
    }
}

/**
    Visualizes a distribution.
    
    Params:
        vals: First entry is an array of x-values and second entry is an array of sampled y-values
*/
function visualize_dist(vals) {
    try {
        // Plot
        let data = [];
        
        // Data points
        let trace = {
          y: vals[1],
          x: vals[0],
          mode: "line",
          line: {shape: 'spline'},
          name: ""
        };
        data.unshift(trace);
        
        let layout = {
            showlegend: false,
        }
        
        Plotly.newPlot("chart", data, layout);
    }
    catch (e) {
        throw(e);
    }
}
/*****************************************************

Util functions.

******************************************************/

/**
    Modified code from:
    https://cmatskas.com/importing-csv-files-using-jquery-and-html5/
*/
$(document).ready(function() {

    // The event listener for the file upload
    if (document.getElementById('txtFileUpload') != null) {
        document.getElementById('txtFileUpload').addEventListener('change', upload, false);

        // Method that checks that the browser supports the HTML5 File API
        function browserSupportFileUpload() {
            let isCompatible = false;
            if (window.File && window.FileReader && window.FileList && window.Blob) {
                isCompatible = true;
            }
            return isCompatible;
        }

        // Method that reads and processes the selected file
        function upload(evt) {
        if (!browserSupportFileUpload()) {
            alert("The File APIs are not fully supported in this browser!");
            }
            else {
                let data = null;
                let file = evt.target.files[0];
                let reader = new FileReader();
                reader.readAsText(file);
                reader.onload = function(event) {
                    let csvData = event.target.result;
                    data = $.csv.toArrays(csvData);
                    if (data && data.length > 0) {
                        // Import data entries
                        let n = data[0].length;
                        // Fill string array
                        let arrs = [];
                        for (let i = 0; i < n; i++){
                            arrs.push("");        
                        }
                        // Add data values to samples array
                        for (let i in data) {
                            let e = data[i];
                            
                            for (let vi in e) {
                                let v = parseFloat(e[vi]);
                                if (!isNaN(v)) {
                                    arrs[vi] += v + ",";
                                }
                            }
                        }
                        
                        for (let i in arrs) {
                            let str = arrs[i];
                            // Remove last comma
                            str = str.substr(0,str.length - 1);
                            // Set data
                            let key = "samp" + (parseInt(i) + 1);
                            let inp = document.getElementById(key);
                            if (inp != null) {
                                inp.value = str;
                            }
                        }
                    }
                    else {
                        alert("No data to import!");
                    }
                };
                reader.onerror = function() {
                    alert("Unable to read " + file.fileName);
                };
            }
        }
    }
});

/**
    Updates the checked test type radio buttons.
*/
function set_test_type() {
    let sel_type = get_param("type");
    if (sel_type != null) {
        document.forms["testtype"][sel_type].checked=true;
    }
}

/**
    Returns the requested URL parameter, or null if not found.
*/
function get_param(name) {
    let param = null;
    
    // Read all params
    let parts = window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi, function(m,key,value) {
        if (key == name) {
            param = value;
        }
    });
    return param;
}


/**
    Converts a comma-separated string to an array of float values.
*/
function stringToFloatArray(str) {
    let arr = str.split(",");
    for (let i in arr) {
        let val = arr[i];
        val = val.trim();
        arr[i] = parseFloat(val);
        if (isNaN(arr[i])) {
            return NaN;
        }
    }
    return arr;
}

/**
    Shows the div element with the specified id.
*/
function show(id) {
    let e = document.getElementById(id);
    if (e != null) {
        e.style.display = "block";
    }
}

/**
    Hides the div element with the specified id.
*/
function hide(id) {
    let e = document.getElementById(id);
    if (e != null) {
        e.style.display = "none";
    }
}

/**
    Toggles visibility of the div element with the specified id.
*/
function toggle(id) {
    let e = document.getElementById(id);
    if (e.style.display == "none") {
        e.style.display = "block";
    }
    else {
        e.style.display = "none";
    }
}

/**
    Fills the descriptive statistics table with sample arrays.
*/
function fill_descriptive_stats(arrs) {
    // Number of samples
    let k = arrs.length;
    for (let s = 1; s <= k; s++) {
        let v = arrs[s-1];
        // Sample mean
        let mean = jStat.mean(v);
        // Sample size
        let n = v.length;
        // Standard deviations for the sample
        let std = jStat.stdev(v, true);
        // Set values
        document.getElementById("n" + s).innerHTML = n;
        document.getElementById("mean" + s).innerHTML = mean.toFixed(2);
        document.getElementById("stdev" + s).innerHTML = std.toFixed(3);
    }
}

/**
    Clears the data fields.
*/
function clear_fields(no) {
    if (no == -1) {
        // Get number of fields
        let no = get_no_samples();
    }
    
    for (let i = 1; i <= no; i++) {
        document.getElementById("samp" + i).value = "";
    }
    
    if (document.getElementById("spmean") != null) {
        document.getElementById("spmean").value = "";
    }
    if (document.getElementById("alpha") != null) {
        document.getElementById("alpha").value = "0.05";
    }
    hide("test_results");
    clear_error();
}

/**
    Shows an error.
*/
function show_error(error_str) {
    hide("test_results");
    
    document.getElementById("error").innerHTML = "<br/><font color='red'><b>ERROR: " + error_str + "</b></font>";
}

/**
    Clears the error field.
*/
function clear_error() {
    document.getElementById("error").innerHTML = "";
}

/**
    Reads the sample arrays from the html elements.
*/
function get_sample_arrays(k) {
    // Array of sample arrays
    let arrs = [];
    for (let s = 1; s <= k; s++) {
        let str = document.getElementById("samp" + s).value;
        let v = stringToFloatArray(str);
        if (!Array.isArray(v)) {
            let id = String.fromCharCode(65 + s - 1);
            throw("Invalid format in sample " + id);
        }
        arrs.push(v);
    }
    return arrs;
}

/**
    Returns the number of samples.
*/
function get_no_samples() {
    let no_str = document.getElementById("no").value;
    no_str = no_str.trim();
    let no = parseInt(no_str);
    if (isNaN(no)) {
        throw("Invalid number of samples: " + no_str);
    }
    return no;
}

/**
    Returns the selected test type.
*/
function get_type() {
    let v = document.querySelector('input[name="type"]:checked').value;
    let type = parseInt(v);
    if (isNaN(type)) {
        throw("Invalid type: " + v);
    }
    return type;
}

/**
    Returns number of sides (2 for two-tailed, 1 for one-tailed tests).
*/
function get_sides() {
    let v = document.querySelector('input[name="hyp"]:checked').value;
    if (v == "1") {
        return 2;
    }
    return 1;
}

/**
    Returns the tail (1 = left, 2 = right) for one-tailed test.
*/
function get_tail() {
    let v = document.querySelector('input[name="hyp"]:checked').value;
    if (v == "2") {
        return 1;
    }
    if (v == "3") {
        return 2;
    }
    return 0;
}

/**
    Returns the significance level (alpha).
*/
function get_alpha() {
    let strA = document.getElementById("alpha").value;
    let alpha = parseFloat(strA);
    if (isNaN(alpha)) {
        throw("Invalid alpha: " + strA);
    }
    // Set alpha
    let e = document.getElementById("sign_level");
    if (e != null) {
        e.innerHTML = alpha;
    }
    return alpha;
}

/**
    Updates number of samples fields.
*/
function update_no_samples() {
    try {
        // Clear error
        clear_error();
        hide("test_results");

        // Get number of samples
        let no = get_no_samples();

        // Samples table
        let html = "";
        for (let i = 1; i <= no; i++) {
            // Convert number to char
            let s = String.fromCharCode(65 + i - 1);

            let samp = document.getElementById("samp" + i);
            let old_vals = "";
            if (samp != null) {
                old_vals = samp.value;
            }

            html += "<tr><td width=110'>Sample " + s + ":</td>";
            html += "<td><input class='sample' name='samp" + s + "' id='samp" + i + "' value='" + old_vals + "'></td></tr>";
        }

        // Set samples
        let tcont = document.getElementById("samples");
        tcont.innerHTML = html;

        // Summary table
        html = "";
        for (let i = 1; i <= no; i++) {
            // Convert number to char
            let s = String.fromCharCode(65 + i - 1);

            html += "<tr><th class='dark'>" + s + "</th>";
            html += "<td class='border' id='n" + i + "'>&nbsp;</td>";
            html += "<td class='border' id='mean" + i + "'>&nbsp;</td>";
            html += "<td class='border' id='stdev" + i + "'>&nbsp;</td>";
            html += "</tr>";
        }

        // Set summary
        tcont = document.getElementById("summary");
        tcont.innerHTML = html;
    }
    catch (e) {
        show_error(e);
    }
}

/**
    Visualizes the samples in a box-and-whiskers plot.
*/
function visualize(arrs, names) {
    let data = [];
    
    for (let i = 0; i < arrs.length; i++) {
        let name = "Sample " + String.fromCharCode(65 + i);
        if (names != null) {
            name = names[i];
        }
        
        let trace = {
          x: arrs[i],
          type: "box",
          name: name
        };
        data.unshift(trace); //push unshift
    }
    
    let layout = {
        legend: {
            traceorder: "reversed",
        }
    }
    
    Plotly.newPlot("chart", data, layout);
}

/**
    Visualizes the samples in a histogram.
*/
function visualize_histogram(arrs) {
    let data = [];
    
    for (let i = 0; i < arrs.length; i++) {
        let name = "Sample " + String.fromCharCode(65 + i);
        
        let trace = {
          x: arrs[i],
          type: "histogram",
          name: name
        };
        data.unshift(trace); //push unshift
    } 
    
    let layout = {
        bargap: 0.05,
        yaxis: {
            title: "y Axis",
            titlefont: {
              family: "Courier New, monospace",
              size: 18,
              color: "#7f7f7f"
            }
        }
    }
    
    Plotly.newPlot("hist", data, layout);
}

/**
    Visualizes linear regression.
*/
function visualize_regression(x_vals, y_vals, a, b) {
    let data = [];
    
    // Calculate line points
    let line_y = [];
    let line_x = [];
    let min_x = Math.floor(jStat.min(x_vals));
    let max_x = Math.ceil(jStat.max(x_vals));
    for (let i = min_x; i <= max_x; i++) {
        let y = a + b * i;
        line_x.push(i);
        line_y.push(y.toFixed(2));
    }
    
    // Data points
    let trace = {
      y: y_vals,
      x: x_vals,
      mode: "markers",
      name: "Data points"
    };
    data.unshift(trace);
    // Line of best fit
    trace = {
      y: line_y,
      x: line_x,
      mode: "lines",
      name: "Line of best fit"
    };
    data.unshift(trace);
    
    let layout = {
        showlegend: true,
        legend: {
            traceorder: "reversed",
        },
        xaxis: {
            title: {
                text: "Sample X"
            }
        },
        yaxis: {
            title: {
                text: "Sample Y"
            }
        }
    }
    
    Plotly.newPlot("chart", data, layout);
}

/**
    Visualizes correlation between two samples.
*/
function visualize_correlation(x_vals, y_vals) {
    let data = [];
    
    // Data points
    let trace = {
      y: y_vals,
      x: x_vals,
      mode: "markers",
      name: ""
    };
    data.unshift(trace);
    
    let layout = {
        showlegend: false,
        xaxis: {
            title: {
                text: "Sample A"
            }
        },
        yaxis: {
            title: {
                text: "Sample B"
            }
        }
    }
    
    Plotly.newPlot("chart", data, layout);
}
/*****************************************************

Main functions for running statistical tests.

******************************************************/

// 2 samples, 2 sided
var diff_2p_2s = "<font color='green'>The means of the samples are different</font>";
var nodiff_2p_2s = "<font color='red'>There is no difference between the means of the samples</font>";
// 2 samples, 1 sided
var diff_2p_1s_AltB = "<font color='green'>The mean of sample A is less than the mean of sample B</font>";
var diff_2p_1s_AgtB = "<font color='green'>The mean of sample A is greater than the mean of sample B</font>";
var nodiff_2p_1s_AltB = "<font color='red'>The mean of sample A is not less than the mean of sample B</font>";
var nodiff_2p_1s_AgtB = "<font color='red'>The mean of sample A is not greater than the mean of sample B</font>";
var nodiff_2p_1s_eq = "<font color='red'>There is no difference between the means of the samples</font>";

// 3 or more samples
var diff_3s = "<font color='green'>There is a difference between the means of the samples</font>";
var nodiff_3s = "<font color='red'>There is no difference between the means of the samples</font>";
// Post-test
var post_diff = "<font color='green'>Difference</font>";
var post_nodiff = "<font color='red'>No difference</font>";

// 1 sample, 2 sided
var diff_2p_1s = "<font color='green'>The mean of the sample is different than the specified mean</font>";
var nodiff_2p_1s = "<font color='red'>There is no difference between the mean of the sample and the specified mean</font>";

// 2 sample, 1 sided
var diff_1p_1s_AltB = "<font color='green'>The mean of the sample is less than the specified mean</font>";
var diff_1p_1s_AgtB = "<font color='green'>The mean of the sample is greater than the specified mean</font>";
var nodiff_1p_1s_AltB = "<font color='red'>The mean of the sample is not less than the specified mean</font>";
var nodiff_1p_1s_AgtB = "<font color='red'>The mean of the sample is not greater than the specified mean</font>";
var nodiff_1p_1s_eq = "<font color='red'>There is no difference between the mean of the sample and the specified mean</font>";

// Equal variances
var eq_var = "<font color='green'>The samples have equal variances</font>";
var uneq_var = "<font color='red'>The sample variances are unequal</font>";

// Normally distributed samples
var norm_dist = "<font color='green'>The sample is normally distributed</font>";
var not_norm_dist = "<font color='red'>The sample is not normally distributed</font>";

// Correlation coefficient
var corr_sign = "<font color='green'>There is a relationship between the samples</font>";
var corr_not_sign = "<font color='red'>There is no relationship between the samples</font>";

function calc_confidence_intervals() {
    // Clear error
    clear_error();
    
    try {
        // Get sample array
        let arrs = get_sample_arrays(1);
        let v = arrs[0];
        
        // Sample mean
        let mean = jStat.mean(v);
        // Sample size
        let n = v.length;
        // Standard deviations for the sample
        let std = jStat.stdev(v, true);
        // Set values
        document.getElementById("n").innerHTML = n;
        document.getElementById("mean").innerHTML = mean.toFixed(2);
        document.getElementById("stdev").innerHTML = std.toFixed(3);
        
        let carrs = [];
        let cnames = [];
        
        let res = confidence_interval(v, 0.90);
        document.getElementById("ci90a").innerHTML = "&nbsp;±" + res[0].toFixed(3) + "&nbsp;";
        document.getElementById("ci90b").innerHTML = "&nbsp;" + res[1].toFixed(3) +" to " + res[2].toFixed(3) + "&nbsp;";
        document.getElementById("ci90c").innerHTML = "&nbsp;<font color='green'>90% certain that the true population mean is between " + res[1].toFixed(3) +" and " + res[2].toFixed(3) + "&nbsp;</font>";
        carrs.push( [jStat.min(v), res[1].toFixed(2), res[1].toFixed(2), jStat.median(v), res[2].toFixed(2), res[2].toFixed(2), jStat.max(v)] );
        cnames.push( "90%" );
        
        res = confidence_interval(v, 0.95);
        document.getElementById("ci95a").innerHTML = "&nbsp;±" + res[0].toFixed(3) + "&nbsp;";
        document.getElementById("ci95b").innerHTML = "&nbsp;" + res[1].toFixed(3) +" to " + res[2].toFixed(3) + "&nbsp;";
        document.getElementById("ci95c").innerHTML = "&nbsp;<font color='green'>95% certain that the true population mean is between " + res[1].toFixed(3) +" and " + res[2].toFixed(3) + "&nbsp;</font>";
        carrs.push( [jStat.min(v), res[1].toFixed(2), res[1].toFixed(2), jStat.median(v), res[2].toFixed(2), res[2].toFixed(2), jStat.max(v)] );
        cnames.push( "95%" );

        res = confidence_interval(v, 0.98);
        document.getElementById("ci98a").innerHTML = "&nbsp;±" + res[0].toFixed(3) + "&nbsp;";
        document.getElementById("ci98b").innerHTML = "&nbsp;" + res[1].toFixed(3) +" to " + res[2].toFixed(3) + "&nbsp;";
        document.getElementById("ci98c").innerHTML = "&nbsp;<font color='green'>98% certain that the true population mean is between " + res[1].toFixed(3) +" and " + res[2].toFixed(3) + "&nbsp;</font>";
        carrs.push( [jStat.min(v), res[1].toFixed(2), res[1].toFixed(2), jStat.median(v), res[2].toFixed(2), res[2].toFixed(2), jStat.max(v)] );
        cnames.push( "98%" );
        
        // Data visualization
        visualize( carrs, cnames );
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function calc_correlation() {
    // Clear error
    clear_error();
    
    try {
        // Get sample arrays
        let arrs = get_sample_arrays(2);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        let res = correlation_R(arrs[0], arrs[1], alpha);
        
        document.getElementById("corr").innerHTML = res[0].toFixed(3);
        document.getElementById("p").innerHTML = res[1].toFixed(5);
        if (res[1] <= alpha) {
            document.getElementById("res").innerHTML = corr_sign;
        }
        else {
            document.getElementById("res").innerHTML = corr_not_sign;
        }
        // Score
        if (res[2] < 0) {
            document.getElementById("score").innerHTML = res[2].toFixed(2) + "<br>Correlation coefficient is significant if T < -" + res[3].toFixed(2);
        }
        else {
            document.getElementById("score").innerHTML = res[2].toFixed(2) + "<br>Correlation coefficient is significant if T > " + res[3].toFixed(2);
        }
        
        // Normally distributed samples
        // Run Shapiro-Wilk test
        if (arrs[0].length < 10) {
            document.getElementById("sw_p_1").innerHTML = "";
            document.getElementById("sw_res_1").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[0], alpha);
            // Set results
            document.getElementById("sw_p_1").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_1").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_1").innerHTML = norm_dist;
            }
        }
        if (arrs[1].length < 10) {
            document.getElementById("sw_p_2").innerHTML = "";
            document.getElementById("sw_res_2").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[1], alpha);
            // Set results
            document.getElementById("sw_p_2").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_2").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_2").innerHTML = norm_dist;
            }
        }
        
        // Data visualization
        //visualize(arrs);
        visualize_correlation(arrs[0], arrs[1]);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function calc_linearregression() {
    // Clear error
    clear_error();
    
    try {
        // Get sample arrays
        let arrs = get_sample_arrays(2);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        let res = linear_regression(arrs[0], arrs[1], alpha);
        
        if (res[1] >= 0) {
            document.getElementById("line").innerHTML = "y = " + res[0].toFixed(3) + " + " + res[1].toFixed(3) + "x";
        }
        else {
            document.getElementById("line").innerHTML = "y = " + res[0].toFixed(3) + " - " + Math.abs(res[1].toFixed(3)) + "x";
        }
        let R2 = res[2];
        let R2p = R2 * 100;
        document.getElementById("r2").innerHTML = R2.toFixed(3) + "<br>" + R2p.toFixed(1) + "% of the variation in B is explained by A";
        document.getElementById("p").innerHTML = res[3].toFixed(5);
        if (res[3] <= alpha) {
            document.getElementById("res").innerHTML = corr_sign;
        }
        else {
            document.getElementById("res").innerHTML = corr_not_sign;
        }
        // Score
        if (res[4] < 0) {
            document.getElementById("score").innerHTML = res[4].toFixed(2) + "<br>Correlation coefficient is significant if T < -" + res[5].toFixed(2);
        }
        else {
            document.getElementById("score").innerHTML = res[4].toFixed(2) + "<br>Correlation coefficient is significant if T > " + res[5].toFixed(2);
        }
        
        // Data visualization
        visualize_regression(arrs[0], arrs[1], res[0], res[1]);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function find_outliers() {
    // Clear error
    clear_error();
    
    try {
        // Get sample array
        let arrs = get_sample_arrays(1);
        let v = arrs[0];
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
    
        let res = outliers_esd(v, alpha);
        let vals = res[1];
        
        // Create result HTML
        let new_v = [];
        let html = "";
        for (let i = 0; i < vals.length; i++) {
            let val = vals[i];
            if (val[1] == false) {
                html += "<font color='green'>" + val[0] + "</font>";
                new_v.push(val[0]);
            }
            else {
                html += "<font color='red'><b>" + val[0] + "</b></font>";
            }
            if (i < vals.length - 1) {
                html += ", ";
            }
            
            if ((i+1) % 20 == 0) {
                html += "<br>";
            }
        }
        document.getElementById("vals").innerHTML = html;
        
        // Create trials HTML
        let tlist = res[0];
        if (tlist.length > 0) {
            let html = "";
            let tlist = res[0];
            for (let i = 0; i < tlist.length; i++) {
                let te = tlist[i];
                html += "<tr>";
                html += "<td class='border'>" + (i+1) + "</td>";
                html += "<td class='border'>" + te[0] + "</td>";
                html += "<td class='border'>" + te[3].toFixed(3) + "</td>";
                html += "<td class='border'>" + te[4].toFixed(3) + "</td>";
                if (te[2] == false) {
                    html += "<td class='border'><font color='green'>No</font></td>";
                }
                else {
                    html += "<td class='border'><font color='red'>Yes</font></td>";
                }
                html += "</tr>";
            }
            document.getElementById("trials").innerHTML = html;
            document.getElementById("restext").innerHTML = "The sample values marked as <font color='red'><b>red</b></font> are potential outliers.";
            
            show("trials_tab");
        }
        else {
            // No outliers found
            document.getElementById("restext").innerHTML = "No outliers found in the sample.";
            hide("trials_tab");
        }
        
        // Data visualization
        visualize( [v, new_v], ["Original sample", "Outliers removed"] );
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_ttest() {
    // Clear error
    clear_error();
    
    try {
        // Get sample arrays
        let arrs = get_sample_arrays(2);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Test type
        let type = get_type();
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Sides (one-tailed or two-tailed) and tail (left or right)
        let sides = get_sides();
        let tail = get_tail();
        
        // Set mean diff
        let mean1 = jStat.mean(arrs[0]);
        let mean2 = jStat.mean(arrs[1]);
        document.getElementById("mean").innerHTML = (mean1 - mean2).toFixed(3);
        
        // Run test
        let res = ttest(arrs[0], arrs[1], sides, tail, alpha, type, 0);
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (sides == 2) {
            if (res[0] <= alpha) {
                document.getElementById("res").innerHTML = diff_2p_2s;
            }
            else {
                document.getElementById("res").innerHTML = nodiff_2p_2s;
            }
            
        }
        else {
            if (res[0] <= alpha) {
                let txt = diff_2p_1s_AltB;
                if (tail == 2) {
                    txt = diff_2p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
            else {
                let txt = nodiff_2p_1s_AltB;
                if (tail == 2) {
                    txt = nodiff_2p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
        }
        
        // Show T-score
        ttest_score(res, sides, tail, "t");
        
        // Independent or paired?
        let dep = false;
        if (type == 3) {
            // Paired
            dep = true;
        }
        
        // Power analysis
        let minN = power_2s(arrs[0], arrs[1], alpha, sides, dep, true);
        document.getElementById("n_low").innerHTML = minN[0];
        document.getElementById("n_medium").innerHTML = minN[1];
        document.getElementById("n_high").innerHTML = minN[2];
        document.getElementById("pwr").innerHTML = minN[3].toFixed(2) + "%";
        
        // Equal variances
        if (type != 3) {
            show("vartest");
            // Run F-test
            let res = ftest(arrs[0], arrs[1], alpha);
            // Set results
            document.getElementById("f_p").innerHTML = res[0].toFixed(5) + "<br>Variances are equal if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("f_res").innerHTML = uneq_var;
            }
            else {
                document.getElementById("f_res").innerHTML = eq_var;
            }
        }
        else {
            hide("vartest");
        }
        
        // Normally distributed samples
        // Run Shapiro-Wilk test
        if (arrs[0].length < 10) {
            document.getElementById("sw_p_1").innerHTML = "";
            document.getElementById("sw_res_1").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[0], alpha);
            // Set results
            document.getElementById("sw_p_1").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_1").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_1").innerHTML = norm_dist;
            }
        }
        if (arrs[1].length < 10) {
            document.getElementById("sw_p_2").innerHTML = "";
            document.getElementById("sw_res_2").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[1], alpha);
            // Set results
            document.getElementById("sw_p_2").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_2").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_2").innerHTML = norm_dist;
            }
        }
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function ttest_score(res, sides, tail, id) {
    if (sides == 2) {
        document.getElementById(id).innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if T is outside ±" + res[2].toFixed(2);
    }
    if (sides == 1 && tail == 1) {
        document.getElementById(id).innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if T < " + res[2].toFixed(2);
    }    
    if (sides == 1 && tail == 2) {
        document.getElementById(id).innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if T > " + res[2].toFixed(2);
    }
}

function run_wilcoxon() {
    // Clear error
    clear_error();
    
    try {
        // Get sample arrays
        let arrs = get_sample_arrays(2);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Test type
        let type = get_type();
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Sides (one-tailed or two-tailed) and tail (left or right)
        let sides = get_sides();
        let tail = get_tail();
        
        // Set mean diff
        let mean1 = jStat.mean(arrs[0]);
        let mean2 = jStat.mean(arrs[1]);
        document.getElementById("mean").innerHTML = (mean1 - mean2).toFixed(3);
        
        // Run test
        let res = wilcoxon(arrs[0], arrs[1], sides, tail, alpha, type, 0);
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (sides == 2) {
            if (res[0] <= alpha) {
                document.getElementById("res").innerHTML = diff_2p_2s;
            }
            else {
                document.getElementById("res").innerHTML = nodiff_2p_2s;
            }
            
        }
        else {
            if (res[0] <= alpha) {
                let txt = diff_2p_1s_AltB;
                if (tail == 2) {
                    txt = diff_2p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
            else {
                let txt = nodiff_2p_1s_AltB;
                if (tail == 2) {
                    txt = nodiff_2p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
        }
        
        // Critical score
        wilcoxon_score(res, sides, tail, type, "r");
        
        // Independent or paired?
        let dep = false;
        if (type == 2) {
            // Paired
            dep = true;
        }
        
        // Power analysis
        let minN = power_2s(arrs[0], arrs[1], alpha, sides, dep, false);
        document.getElementById("n_low").innerHTML = minN[0];
        document.getElementById("n_medium").innerHTML = minN[1];
        document.getElementById("n_high").innerHTML = minN[2];
        document.getElementById("pwr").innerHTML = minN[3].toFixed(2) + "%";
        
        // Normally distributed samples
        // Run Shapiro-Wilk test
        if (arrs[0].length < 10) {
            document.getElementById("sw_p_1").innerHTML = "";
            document.getElementById("sw_res_1").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[0], alpha);
            // Set results
            document.getElementById("sw_p_1").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_1").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_1").innerHTML = norm_dist;
            }
        }
        if (arrs[1].length < 10) {
            document.getElementById("sw_p_2").innerHTML = "";
            document.getElementById("sw_res_2").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[1], alpha);
            // Set results
            document.getElementById("sw_p_2").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_2").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_2").innerHTML = norm_dist;
            }
        }
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function wilcoxon_score(res, sides, tail, type, id) {
    // Rank-sum
    if (type == 1) {
        if (res[4] == true) {
            // Critical R-score
            if (sides == 2) {
                document.getElementById(id).innerHTML = "R: " + res[1] + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if R is outside " + res[2] + " ≤ R ≤ " + res[3];
            }
            if (sides == 1 && tail == 1) {
                document.getElementById(id).innerHTML = "R: " + res[1] + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if R < " + res[2];
            }    
            if (sides == 1 && tail == 2) {
                document.getElementById(id).innerHTML = "R: " + res[1] + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if R > " + res[3];
            }
        }
        else {
            // Critical Z-score
            if (sides == 2) {
                document.getElementById(id).innerHTML = "Z: " + res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Z is outside ±" + res[3].toFixed(2);
            }
            if (sides == 1 && tail == 1) {
                document.getElementById(id).innerHTML = "Z: " + res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Z < " + res[3].toFixed(2);
            }    
            if (sides == 1 && tail == 2) {
                document.getElementById(id).innerHTML = "Z: " + res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Z > " + res[3].toFixed(2);
            }
        }
    }
    // Signed-ranks
    else {
        if (res[3] == true) {
            // Critical W-score
            document.getElementById(id).innerHTML = "W: " + res[1] + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if W ≤ " + res[2];
        }
        else {
            // Critical Z-score
            document.getElementById(id).innerHTML = "Z: " + res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Z ≤ " + res[2].toFixed(2);
        }
    }
} 

function run_anova() {
    try {
        // Clear error
        clear_error();
        
        // Returns the number of samples
        let no = get_no_samples();
        
        // Get sample arrays
        let arrs = get_sample_arrays(no);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Run ANOVA
        let res = anova(arrs, alpha);
        
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (res[0] <= alpha) {
            document.getElementById("res").innerHTML = diff_3s;
        }
        else {
            document.getElementById("res").innerHTML = nodiff_3s;
        }
        // F-value
        document.getElementById("f").innerHTML = res[5].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if F > " + res[6].toFixed(2);
        
        // Run post-test if significant differences
        if (res[0] <= alpha) {
            // Post-test to do
            let e = document.getElementById("posttest");
            let type = e.options[e.selectedIndex].value;
            let name = e.options[e.selectedIndex].text;
            
            // Alpha correction
            e = document.getElementById("correction");
            let correction = e.options[e.selectedIndex].value;
            
            // Run post-test
            let P_vals = null;
            if (type == "scheffes") {
                P_vals = scheffes_posttest(arrs, alpha, correction, res);
            }
            else if (type == "tukeys") {
                P_vals = tukeys_posttest(arrs, alpha, correction, res);
            }
            
            let html = "<h3 class='f14'>" + name + "</h3>";
            html += "<table class='border'>";
            html += "<thead><tr><th class='dark' width='70'>Pair</th><th class='dark' width='250'>Score</th><th class='dark' width='80'>P-value</th><th class='dark' width='80'>Alpha</th><th class='dark' width='150'>Result</th></tr></thead><tbody>";
            
            for (let i = 0; i < P_vals.length; i++) {
                let res = P_vals[i];
                html += "<tr>";
                
                // Convert number to char
                let s1 = String.fromCharCode(65 + res[0] - 1);
                let s2 = String.fromCharCode(65 + res[1] - 1);
                
                // Pair and P-value
                html += "<th class='dark'>" + s1 + " - " + s2 + "</th>";
                html += "<td class='border'>" + res[6] + "-score: " + res[4].toFixed(3) + ", difference if " + res[6] + " > " + res[5].toFixed(3) + "</td>";
                html += "<td class='border'>" + res[2].toFixed(5) + "</td>";
                html += "<td class='border'>" + res[3].toFixed(3) + "</td>";
                
                // Result
                let res_str = post_nodiff;
                if (res[2] <= res[3]) {
                    res_str = post_diff;
                }
                html += "<td class='border'>" + res_str + "</td>";
                html += "</tr>";
            }
            html += "</tbody></table>";
            
            // Set result table
            document.getElementById("post").innerHTML = html;
            
            show("postres");
        }
        else {
            hide("postres");
        }
        
        // Power analysis
        let minN = power_anova(arrs, alpha, correction, false, true);
        document.getElementById("n_low").innerHTML = minN[0];
        document.getElementById("n_medium").innerHTML = minN[1]
        document.getElementById("n_high").innerHTML = minN[2];
        document.getElementById("pwr").innerHTML = minN[3].toFixed(2) + "%";
        
        // Equal variances
        // Run Bartlett's test
        res = bartlett(arrs, alpha);
        // Set results
        document.getElementById("v_p").innerHTML = res[0].toFixed(5) + "<br>Variances are equal if P > " + alpha;
        if (res[0] <= alpha) {
            document.getElementById("v_res").innerHTML = uneq_var;
        }
        else {
            document.getElementById("v_res").innerHTML = eq_var;
        }
        
        // Normally distributed samples
        // Run Shapiro-Wilk test
        let html = "";
        for (let si = 0; si < no; si++) {
            let res = shapiro_wilk(arrs[si], alpha);
            html += "<tr><td class='dark' colspan='2'><center><b>Sample " + String.fromCharCode(65 + si); + "</b></center></td></tr>";
            html += "<tr><td class='dark' width='100'><b>Result:</b></td><td class='border'>";
            
            if (arrs[si].length < 10) {
                html += "<font color='blue'>Requires a sample size of at least 10</font>";
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
            }
            else {
                html += res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
                if (res[0] <= alpha) {
                    html += not_norm_dist;
                }
                else {
                    html += norm_dist;
                }
            }
            html += "</td></tr>";
        }
        document.getElementById("normtest_cont").innerHTML = html;
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_rm_anova() {
    try {
        // Clear error
        clear_error();
        
        // Returns the number of samples
        let no = get_no_samples();
        
        // Get sample arrays
        let arrs = get_sample_arrays(no);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Run Repeated Measures ANOVA
        let res = rm_anova(arrs, alpha);
        
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (res[0] <= alpha) {
            document.getElementById("res").innerHTML = diff_3s;
        }
        else {
            document.getElementById("res").innerHTML = nodiff_3s;
        }
        // F-value
        document.getElementById("f").innerHTML = res[5].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if F > " + res[6].toFixed(2);
        
        // Run post-test if significant differences
        if (res[0] <= alpha) {
            // Post-test to do
            let e = document.getElementById("posttest");
            let type = e.options[e.selectedIndex].value;
            let name = e.options[e.selectedIndex].text;
            
            // Alpha correction
            e = document.getElementById("correction");
            let correction = e.options[e.selectedIndex].value;
            
            // Run post-test
            let P_vals = null;
            if (type == "tukeys") {
                P_vals = tukeys_rm_posttest(arrs, alpha, correction, res);
            }
            
            let html = "<h3 class='f14'>" + name + "</h3>";
            html += "<table class='border'>";
            html += "<thead><tr><th class='dark' width='70'>Pair</th><th class='dark' width='250'>Score</th><th class='dark' width='80'>P-value</th><th class='dark' width='80'>Alpha</th><th class='dark' width='150'>Result</th></tr></thead><tbody>";
            
            for (let i = 0; i < P_vals.length; i++) {
                let res = P_vals[i];
                html += "<tr>";
                
                // Convert number to char
                let s1 = String.fromCharCode(65 + res[0] - 1);
                let s2 = String.fromCharCode(65 + res[1] - 1);
                
                // Pair and P-value
                html += "<th class='dark'>" + s1 + " - " + s2 + "</th>";
                html += "<td class='border'>" + res[6] + "-score: " + res[4].toFixed(3) + ", difference if " + res[6] + " > " + res[5].toFixed(3) + "</td>";
                html += "<td class='border'>" + res[2].toFixed(5) + "</td>";
                html += "<td class='border'>" + res[3].toFixed(3) + "</td>";
                
                // Result
                let res_str = post_nodiff;
                if (res[2] <= res[3]) {
                    res_str = post_diff;
                }
                html += "<td class='border'>" + res_str + "</td>";
                html += "</tr>";
            }
            html += "</tbody></table>";
            
            // Set result table
            document.getElementById("post").innerHTML = html;
            
            show("postres");
        }
        else {
            hide("postres");
        }
        
        // Normally distributed samples
        // Run Shapiro-Wilk test
        let html = "";
        for (let si = 0; si < no; si++) {
            let res = shapiro_wilk(arrs[si], alpha);
            html += "<tr><td class='dark' colspan='2'><center><b>Sample " + String.fromCharCode(65 + si); + "</b></center></td></tr>";
            html += "<tr><td class='dark' width='100'><b>Result:</b></td><td class='border'>";
            
            if (arrs[si].length < 10) {
                html += "<font color='blue'>Requires a sample size of at least 10</font>";
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
            }
            else {
                html += res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
                if (res[0] <= alpha) {
                    html += not_norm_dist;
                }
                else {
                    html += norm_dist;
                }
            }
            html += "</td></tr>";
        }
        document.getElementById("normtest_cont").innerHTML = html;
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_kruskalwallis() {
    try {
        // Clear error
        clear_error();
        
        // Returns the number of samples
        let no = get_no_samples();
        
        // Get sample arrays
        let arrs = get_sample_arrays(no);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Run Kruskal-Wallis
        let res = kruskalwallis(arrs, alpha);
        
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (res[0] <= alpha) {
            document.getElementById("res").innerHTML = diff_3s;
        }
        else {
            document.getElementById("res").innerHTML = nodiff_3s;
        }
        // H-value
        document.getElementById("h").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if H > " + res[2].toFixed(2);
        
        // Run post-test if significant differences
        if (res[0] <= alpha) {
            // Post-test to do
            let e = document.getElementById("posttest");
            let type = e.options[e.selectedIndex].value;
            let name = e.options[e.selectedIndex].text;
            
            // Alpha correction
            e = document.getElementById("correction");
            let correction = e.options[e.selectedIndex].value;
            
            // Run post-test
            let P_vals = null;
            if (type == "dunn") {
                P_vals = dunns_posttest(arrs, alpha, correction);
            }
            // Run post-test
            if (type == "wilcoxon") {
                P_vals = wilcoxon_posttest(arrs, alpha, correction);
            }
            
            let html = "<h3 class='f14'>" + name + "</h3>";
            html += "<table class='border'>";
            html += "<thead><tr><th class='dark' width='70'>Pair</th><th class='dark' width='250'>Score</th><th class='dark' width='80'>P-value</th><th class='dark' width='80'>Alpha</th><th class='dark' width='150'>Result</th></tr></thead><tbody>";
            
            for (let i = 0; i < P_vals.length; i++) {
                let res = P_vals[i];
                html += "<tr>";
                
                // Convert number to char
                let s1 = String.fromCharCode(65 + res[0] - 1);
                let s2 = String.fromCharCode(65 + res[1] - 1);
                
                // Pair, P-value and alpha
                html += "<th class='dark'>" + s1 + " - " + s2 + "</th>";
                html += "<td class='border'>" + res[6] + "-score: " + res[4].toFixed(3) + ", difference if " + res[6] + " > " + res[5].toFixed(3) + "</td>";
                html += "<td class='border'>" + res[2].toFixed(5) + "</td>";
                html += "<td class='border'>" + res[3].toFixed(3) + "</td>";
                
                // Result
                let res_str = post_nodiff;
                if (res[2] <= res[3]) {
                    res_str = post_diff;
                }
                html += "<td class='border'>" + res_str + "</td>";
                
                html += "</tr>";
            }
            
            html += "</tbody></table>";
            
            // Set result table
            document.getElementById("post").innerHTML = html;
            
            show("postres");
        }
        else {
            hide("postres");
        }
        
        // Power analysis
        // Power analysis
        let minN = power_anova(arrs, alpha, correction, false, false);
        document.getElementById("n_low").innerHTML = minN[0];
        document.getElementById("n_medium").innerHTML = minN[1]
        document.getElementById("n_high").innerHTML = minN[2];
        document.getElementById("pwr").innerHTML = minN[3].toFixed(2) + "%";
        
        // Normally distributed samples
        // Run Shapiro-Wilk test
        let html = "";
        for (let si = 0; si < no; si++) {
            let res = shapiro_wilk(arrs[si], alpha);
            html += "<tr><td class='dark' colspan='2'><center><b>Sample " + String.fromCharCode(65 + si); + "</b></center></td></tr>";
            html += "<tr><td class='dark' width='100'><b>Result:</b></td><td class='border'>";
            
            if (arrs[si].length < 10) {
                html += "<font color='blue'>Requires a sample size of at least 10</font>";
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
            }
            else {
                html += res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
                if (res[0] <= alpha) {
                    html += not_norm_dist;
                }
                else {
                    html += norm_dist;
                }
            }
            html += "</td></tr>";
        }
        document.getElementById("normtest_cont").innerHTML = html;
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_friedman() {
    try {
        // Clear error
        clear_error();
        
        // Returns the number of samples
        let no = get_no_samples();
        
        // Get sample arrays
        let arrs = get_sample_arrays(no);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Run Kruskal-Wallis
        let res = friedman(arrs, alpha);
        
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (res[0] <= alpha) {
            document.getElementById("res").innerHTML = diff_3s;
        }
        else {
            document.getElementById("res").innerHTML = nodiff_3s;
        }
        // H-value
        document.getElementById("q").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Q > " + res[2].toFixed(2);
        
        // Run post-test if significant differences
        if (res[0] <= alpha) {
            // Post-test to do
            let e = document.getElementById("posttest");
            let type = e.options[e.selectedIndex].value;
            let name = e.options[e.selectedIndex].text;
            
            // Alpha correction
            e = document.getElementById("correction");
            let correction = e.options[e.selectedIndex].value;
            
            // Run post-test
            let P_vals = null;
            if (type == "wilcoxon") {
                P_vals = wilcoxon_sr_posttest(arrs, alpha, correction);
            }
            
            let html = "<h3 class='f14'>" + name + "</h3>";
            html += "<table class='border'>";
            html += "<thead><tr><th class='dark' width='70'>Pair</th><th class='dark' width='250'>Score</th><th class='dark' width='80'>P-value</th><th class='dark' width='80'>Alpha</th><th class='dark' width='150'>Result</th></tr></thead><tbody>";
            
            for (let i = 0; i < P_vals.length; i++) {
                let res = P_vals[i];
                html += "<tr>";
                
                // Convert number to char
                let s1 = String.fromCharCode(65 + res[0] - 1);
                let s2 = String.fromCharCode(65 + res[1] - 1);
                
                // Pair, P-value and alpha
                html += "<th class='dark'>" + s1 + " - " + s2 + "</th>";
                html += "<td class='border'>" + res[6] + "-score: " + res[4].toFixed(3) + ", difference if " + res[6] + " < " + res[5].toFixed(3) + "</td>";
                html += "<td class='border'>" + res[2].toFixed(5) + "</td>";
                html += "<td class='border'>" + res[3].toFixed(3) + "</td>";
                
                // Result
                let res_str = post_nodiff;
                if (res[2] <= res[3]) {
                    res_str = post_diff;
                }
                html += "<td class='border'>" + res_str + "</td>";
                
                html += "</tr>";
            }
            
            html += "</tbody></table>";
            
            // Set result table
            document.getElementById("post").innerHTML = html;
            
            show("postres");
        }
        else {
            hide("postres");
        }
        
        // Normally distributed samples
        // Run Shapiro-Wilk test
        let html = "";
        for (let si = 0; si < no; si++) {
            
            let res = shapiro_wilk(arrs[si], alpha);
            html += "<tr><td class='dark' colspan='2'><center><b>Sample " + String.fromCharCode(65 + si); + "</b></center></td></tr>";
            html += "<tr><td class='dark' width='100'><b>Result:</b></td><td class='border'>";
            
            if (arrs[si].length < 10) {
                html += "<font color='blue'>Requires a sample size of at least 10</font>";
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
            }
            else {
                html += res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
                html += "</td></tr><tr><td class='dark'><b>P-value:</b></td><td class='border'>";
                if (res[0] <= alpha) {
                    html += not_norm_dist;
                }
                else {
                    html += norm_dist;
                }
            }
            html += "</td></tr>";
        }
        document.getElementById("normtest_cont").innerHTML = html;
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_ttest_single() {
    // Clear error
    clear_error();
    
    try {
        // Get sample array
        let arrs = get_sample_arrays(1);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Sides (one-tailed or two-tailed) and tail (left or right)
        let sides = get_sides();
        let tail = get_tail();
        
        // Get specified mean
        let mean_str = document.getElementById("spmean").value;
        let mean = parseFloat(mean_str);
        if (isNaN(mean)) {
            throw("Invalid mean: " + mean_str);
        }
        
        // Set mean diff
        let mean1 = jStat.mean(arrs[0]);
        document.getElementById("mean").innerHTML = (mean1 - mean).toFixed(3);
        
        // Run test
        let res = ttest_single(arrs[0], mean, sides, tail, alpha);
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (sides == 2) {
            if (res[0] <= alpha) {
                document.getElementById("res").innerHTML = diff_2p_1s;
            }
            else {
                document.getElementById("res").innerHTML = nodiff_2p_1s;
            }
            
        }
        else {
            if (res[0] <= alpha) {
                let txt = diff_1p_1s_AltB;
                if (tail == 2) {
                    txt = diff_1p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
            else {
                let txt = nodiff_1p_1s_AltB;
                if (tail == 2) {
                    txt = nodiff_1p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
        }
        
        // Show T-score
        ttest_score(res, sides, tail, "t");
        
        // Power analysis
        let minN = power_1s(arrs[0], mean, alpha, sides, true);
        document.getElementById("n_low").innerHTML = minN[0];
        document.getElementById("n_medium").innerHTML = minN[1];
        document.getElementById("n_high").innerHTML = minN[2];
        document.getElementById("pwr").innerHTML = minN[4].toFixed(2) + "%";
        
        // Normally distributed sample
        // Run Shapiro-Wilk test
        if (arrs[0].length < 10) {
            document.getElementById("sw_p_1").innerHTML = "";
            document.getElementById("sw_res_1").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[0], alpha);
            // Set results
            document.getElementById("sw_p_1").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_1").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_1").innerHTML = norm_dist;
            }
        }
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_wilcoxon_single() {
    // Clear error
    clear_error();
    
    try {
        // Get sample array
        let arrs = get_sample_arrays(1);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Sides (one-tailed or two-tailed) and tail (left or right)
        let sides = get_sides();
        let tail = get_tail();
        
        // Get specified mean
        let mean_str = document.getElementById("spmean").value;
        let mean = parseFloat(mean_str);
        if (isNaN(mean)) {
            throw("Invalid mean: " + mean_str);
        }
        
        // Set mean diff
        let mean1 = jStat.mean(arrs[0]);
        document.getElementById("mean").innerHTML = (mean1 - mean).toFixed(3);
        
        // Run test
        let res = wilcoxon_signedranks_single(arrs[0], mean, sides, tail, alpha);
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (sides == 2) {
            if (res[0] <= alpha) {
                document.getElementById("res").innerHTML = diff_2p_1s;
            }
            else {
                document.getElementById("res").innerHTML = nodiff_2p_1s;
            }
            
        }
        else {
            if (res[0] <= alpha) {
                let txt = diff_1p_1s_AltB;
                if (tail == 2) {
                    txt = diff_1p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
            else {
                let txt = nodiff_1p_1s_AltB;
                if (tail == 2) {
                    txt = nodiff_1p_1s_AgtB;
                }
                document.getElementById("res").innerHTML = txt;
            }
        }
        
        // Critical score
        wilcoxon_score(res, sides, tail, 2, "r");
        
        // Power analysis
        let minN = power_1s(arrs[0], mean, alpha, sides, false);
        document.getElementById("n_low").innerHTML = minN[0];
        document.getElementById("n_medium").innerHTML = minN[1];
        document.getElementById("n_high").innerHTML = minN[2];
        document.getElementById("pwr").innerHTML = minN[4].toFixed(2) + "%";
        
        // Normally distributed sample
        // Run Shapiro-Wilk test
        if (arrs[0].length < 10) {
            document.getElementById("sw_p_1").innerHTML = "";
            document.getElementById("sw_res_1").innerHTML = "<font color='blue'>Requires a sample size of at least 10</font>";
        }
        else {
            let res = shapiro_wilk(arrs[0], alpha);
            // Set results
            document.getElementById("sw_p_1").innerHTML = res[0].toFixed(5) + "<br>Normally distributed if P > " + alpha;
            if (res[0] <= alpha) {
                document.getElementById("sw_res_1").innerHTML = not_norm_dist;
            }
            else {
                document.getElementById("sw_res_1").innerHTML = norm_dist;
            }
        }
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_ftest() {
    // Clear error
    clear_error();
    
    try {
        // Get sample arrays
        let arrs = get_sample_arrays(2);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Run F-test
        let res = ftest(arrs[0], arrs[1], alpha);
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (res[0] <= alpha) {
            document.getElementById("res").innerHTML = uneq_var;
        }
        else {
            document.getElementById("res").innerHTML = eq_var;
        }
        
        // Critical F
        if (res[1] < 0) {
            document.getElementById("f").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if F < -" + res[2].toFixed(2);
        }
        else {
            document.getElementById("f").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if F > " + res[2].toFixed(2);
        }
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

function run_bartlett() {
    try {
        // Clear error
        clear_error();
        
        // Returns the number of samples
        let no = get_no_samples();
        
        // Get sample arrays
        let arrs = get_sample_arrays(no);
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Run Bartlett's test for equal variances
        let res = bartlett(arrs, alpha);
        
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (res[0] <= alpha) {
            document.getElementById("res").innerHTML = uneq_var;
        }
        else {
            document.getElementById("res").innerHTML = eq_var;
        }
        
        // Critical F
        if (res[1] < 0) {
            document.getElementById("chi2").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Chi-square < -" + res[2].toFixed(2);
        }
        else {
            document.getElementById("chi2").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Chi-square > " + res[2].toFixed(2);
        }
        
        // Data visualization
        visualize(arrs);
        
        show("test_results");
    }
    catch (e) {
        show_error(e);
    }
}

/**
    See explanation and example here:
    http://www.real-statistics.com/tests-normality-and-symmetry/statistical-tests-normality-symmetry/shapiro-wilk-expanded-test/
*/
function run_shapirowilk() {
    // Clear error
    clear_error();
    
    try {
        // Get sample arrays
        let arrs = get_sample_arrays(1);
        
        // Copy array
        let vcopy = [];
        let v1 = arrs[0];
        for (let i = 0; i < v1.length; i++) {
            vcopy.push(v1[i]);
        }
        
        
        // Fill in descriptice statistics
        fill_descriptive_stats(arrs);
        
        // Significance level, alpha
        let alpha = get_alpha();
        
        // Check size
        if (arrs[0].length < 3) {
            show_error("Sample size should be at least 10");
            return;
        }
        
        // Run Shapiro-Wilk test
        let res = shapiro_wilk(arrs[0], alpha);
        
        // Set results
        document.getElementById("p").innerHTML = res[0].toFixed(5);
        if (res[0] <= alpha) {
            document.getElementById("res").innerHTML = not_norm_dist;
        }
        else {
            document.getElementById("res").innerHTML = norm_dist;
        }
        // Set skew and kurtosis
        document.getElementById("skew").innerHTML = res[3].toFixed(5);
        document.getElementById("kurtosis").innerHTML = res[4].toFixed(5);
        document.getElementById("w").innerHTML = res[5].toFixed(5);
        
        // Critical Z
        if (res[1] < 0) {
            document.getElementById("z").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Z < -" + res[2].toFixed(2);
        }
        else {
            document.getElementById("z").innerHTML = res[1].toFixed(2) + "<br>Reject H<sub>0</sub> (accept H<sub>1</sub>) if Z > " + res[2].toFixed(2);
        }
        
        // Data visualization
        visualize(arrs);
        visualize_histogram(arrs);
        
        show("test_results");
        
    }
    catch (e) {
        show_error(e);
    }
}
/*****************************************************

Demonstration of central limit theorem.

******************************************************/

// Frequency of drawn sums
var freq = new Array(25).fill(0);
// Number of drawn samples
var no = 0;
// X-labels
var x_labels = [-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12];

/**
    Draws 1000 random samples and updates the bar chart.
*/
function iterate() {
    // Draw 1000 samples
    for (let i = 0; i < 1000; i++) {
        let sum = 0;
        // Each sample is the sum of five dice rolls
        for (let d = 0; d < 6; d++) {
            sum += Math.floor(Math.random() * 5) + 1;
        }
        // Update the frequency table for the drawn sum
        freq[sum - 6] += 1;
        
    }
    
    // Increase number of drawn samples
    no += 1000;
    document.getElementById("no").innerHTML = no;
    
    update_chart();
}

/**
    Resets the demo.
*/
function reset() {
    freq = new Array(25).fill(0);
    no = 0;
    document.getElementById("no").innerHTML = no;
    
    update_chart();
}

/**
    Updates the bar chart.
*/
function update_chart() {
    // Generate the bar chart
    let data = [];
    
    // Data points
    let trace = {
      y: freq,
      x: x_labels,
      mode: "line",
      line: {shape: 'spline'},
      name: ""
    };
    data.unshift(trace);
    
    let layout = {
        showlegend: false,
    }
    
    Plotly.newPlot("chart", data, layout);
}
