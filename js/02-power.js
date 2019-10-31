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
