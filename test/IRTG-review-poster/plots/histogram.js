// set the dimensions and margins of the graph
const margin = {top: 10, right: 30, bottom: 50, left: 70},
    width = 460 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;

// append the svg object to the body of the page
const svg = d3.select(".canvas")
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          `translate(${margin.left},${margin.top})`);

// get the data
//d3.csv("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/1_OneNum.csv")
d3.csv("exec_times.csv").then( function(data) {
    console.log(data)

  // X axis: scale and draw:
  const x = d3.scaleLinear()
      .domain([0, d3.max(data, (d) => {
          return d.time
      })])     // can use this instead of 1000 to have the max of data: d3.max(data, function(d) { return +d.price })
      .range([0, width]);
  svg.append("g")
      .attr("transform", `translate(0, ${height})`)
      .call(d3.axisBottom(x));

  // set the parameters for the histogram
  const histogram = d3.histogram()
      .value(function(d) { return d.time; })   // I need to give the vector of value
      .domain(x.domain())  // then the domain of the graphic
      .thresholds(x.ticks(10)); // then the numbers of bins

  // And apply this function to data to get the bins
  const bins = histogram(data);
  console.log(bins)

  // Y axis: scale and draw:
  const y = d3.scaleLinear()
      .range([height, 0]);
      y.domain([0, d3.max(bins, function(d) { return d.length; })]);   // d3.hist has to be called before the Y axis obviously
  svg.append("g")
      .call(d3.axisLeft(y));

  // append the bar rectangles to the svg element
  svg.selectAll("rect")
      .data(bins)
      .join("rect")
        .attr("x", 1)
    .attr("transform", function(d) { return `translate(${x(d.x0)} , ${y(d.length)})`})
        .attr("width", function(d) { return x(d.x1) - x(d.x0) -1})
        .attr("height", function(d) { return height - y(d.length); })
        .style("fill", "#59BAE7")

    
    //x-axis lable
    svg.append("text")
        .attr("text-anchor", "end")
        .attr("x", width/2. + margin.left)
        .attr("y", height + margin.top + 40)
        .text("Execution time (ms)")
        .style('font-size', '18px') 

    // Y axis label:
    svg.append("text")
        .attr("text-anchor", "end")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left+20)
        .attr("x", -margin.top -height/3.)
        .style('font-size', '18px') 
        .text("No. of measurements")



    svg.selectAll('rect')
        .on('mouseover', (e, d) => {
        //tip.show(e,d);
        //console.log("hi")
        handleMouseOver(e,d);
    })
        .on('mouseout', (e,d) => {
        //tip.hide();
        handleMouseOut(e,d);
    })
    .on('click', handleClick)
});

const colorPallet = ['#feb236', '#d64161', '#ff7b25','#6b5b95' ]
var currentColorId =  0

const handleMouseOver = (e,d) => {
    //console.log(e.currentTarget);

    d3.select(e.currentTarget)
        .transition('changeSliceFill').duration(300) // name the transition so that they do not interfere with other transitioons
            .style('fill', colorPallet[currentColorId])

}

const handleMouseOut = (e,d) => {

    //reset the color when mouse hovers out
    d3.select(e.currentTarget)
        .transition('changeSliceFill').duration(100) // same name as as hover as they belong to same type
            .style('fill', '#59BAE7')

}

const handleClick = (e,d) => {
    //go to database and delete the slice
    console.log(d.map(d => d.time));
    const algs = ["A", "A", "B"];
    console.log([...new Set(algs)]);

    d3.select(e.currentTarget)
    .transition('changeSliceFill').duration(300) // name the transition so that they do not interfere with other transitioons
        .style('fill', colorPallet[currentColorId])

    d3.select(e.currentTarget)
        .on('mouseout', null)
        .on('mouseover', null)


    if(currentColorId <  3)
        currentColorId += 1

}


function saveSvg(svgEl, name) {
    svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    var svgData = svgEl.outerHTML;
    var preface = '<?xml version="1.0" standalone="no"?>\r\n';
    var svgBlob = new Blob([preface, svgData], {type:"image/svg+xml;charset=utf-8"});
    var svgUrl = URL.createObjectURL(svgBlob);
    var downloadLink = document.createElement("a");
    downloadLink.href = svgUrl;
    downloadLink.download = name;
    document.body.appendChild(downloadLink);
    downloadLink.click();
    document.body.removeChild(downloadLink);
}