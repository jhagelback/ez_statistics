
/*****************************************************

Main functions for running statistical tests.

******************************************************/

/**
    Global status texts.
*/

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

/**
    Runs the calculate confidence intervals test.
    See: confintervals.html
*/
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

/**
    Runs the correlation test.
    See: correlation.html
*/
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

/**
    Runs the linear regression test.
    See: linearregression.html
*/
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

/**
    Runs the ESD test to find outliers.
    See: outliers.html
*/
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

/**
    Runs the two-sample t-tests.
    See: ttest.html
*/
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

/**
    Helper function for showing t-test scores.
*/
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

/**
    Runs the two-sample Wilcoxon tests.
    See: wilcoxon.html
*/
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

/**
    Helper function for showing the Wilcoxon scores.
*/
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

/**
    Runs the One-way ANOVA test.
    See: anova.html
*/
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

/**
    Runs the Repeated Measures ANOVA test.
    See: anova_rm.html
*/
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

/**
    Runs the Kruskal-Wallis test.
    See: kruskalwallis.html
*/
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

/**
    Runs the Friedman test.
    See: friedman.html
*/
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

/**
    Runs the single-sample t-test.
    See: ttest_single.html
*/
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

/**
    Runs the single-sample Wilcoxon test.
    See: wilcoxon_single.html
*/
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

/**
    Runs the F-test for equal variances.
    See: ftest.html
*/
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

/**
    Runs the Bartlett's test for equal variances.
    See: bartlett.html
*/
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
    Runs the Shapiro-Wilk Expanded test for normality check.
    See: shapiro_wilk.html

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
