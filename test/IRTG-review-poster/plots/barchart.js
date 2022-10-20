const svg = d3.select('.canvas')
    .append('svg')
    .attr('width',600)
    .attr('height',600)

    //create a group => creates margins from svg container
    const margin = {top:20, right:20, bottom:100, left:100}; // values are in pixel
    const graphWidth = 600 - margin.left - margin.right;
    const graphHeight = 600 - margin.top - margin.bottom;

    // create margins from the container
    const graph = svg.append('g')
        .attr('width',graphWidth)
        .attr('height', graphHeight)
        .attr('transform', `translate(${margin.left},${margin.top})`);
        //.attr('transform','rotate(90)');

    // group for x axis and y axis elements
    const xAxisGroup = graph.append('g')
        .attr('transform', `translate(0, ${graphHeight})`); // should be in bottom
    const yAxisGroup = graph.append('g')


    d3.json('menu.json').then(data => {

        const y = d3.scaleLinear()
            .domain([0,d3.max(data, d=>d.orders)]) // min and max val of input data
            .range([graphHeight,0]) // 0 to max height of graph. Scale is reversed.
    

       // console.log(y(900));

        // tells where the bars start
        const x = d3.scaleBand()
            .domain(data.map(item => item.name)) 
            .range([0,500])
            .paddingInner(0.2) // space between the bars (between 0 and 1)
            .paddingOuter(0.3); // space from edges of the graph area
    
    
        // join the data to rects
        const rects = graph.selectAll('rect')
            .data(data)
    
        rects.attr('width',x.bandwidth) //bandwidth() wiill invoke the function. But we are just passing the function
            .attr('height', d=> graphHeight - y(d.orders))
            .attr('fill','orange')
            .attr('x', d=>x(d.name)) // not the best way of doing it. use band scale - next lecture
            .attr('y', d=>y(d.orders)); // top edge of the rect
    
        rects.enter()
            .append('rect')
            .attr('width',x.bandwidth)
            .attr('height', d=> graphHeight - y(d.orders))
            .attr('fill','orange')
            .attr('x', (d) => x(d.name))
            .attr('y', d=>y(d.orders));

        // create axis
        const xAxis = d3.axisBottom(x); // ticks come from bottoom (ie elements "veg soup" are at bottom of the axis)
        const yAxis = d3.axisLeft(y)
            .ticks(3) // number of vals in y axis. this is optional arg. Not always exact
            .tickFormat(d=>d+ ' orders') // some formatting
        
        xAxisGroup.call(xAxis);
        yAxisGroup.call(yAxis);

        xAxisGroup.selectAll('text')
            .attr('transform', 'rotate(-40)')
            .attr('text-anchor', 'end') //start, middle of end. Take the end of text and rotate it.
            .attr('fill', 'orange');

        //saveSvg(document.querySelector('svg'),"barchart.svg")
    
    })


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


