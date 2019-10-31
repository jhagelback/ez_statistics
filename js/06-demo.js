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
