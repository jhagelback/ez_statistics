
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
