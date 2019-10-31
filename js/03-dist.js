
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
