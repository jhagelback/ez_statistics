
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
