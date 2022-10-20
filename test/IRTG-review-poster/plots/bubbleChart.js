const data = [
    { name: 'news', parent: '' },

    { name: 'Variants {19,17,16}', parent: 'news', amount:17 },
    { name: 'Variants {18,15,14}', parent: 'news', amount:20 },
    { name: 'Variant {2}', parent: 'news', amount:2 },

    { name: 'Variant {1}', parent: 'news', amount: 7 },
    { name: 'Variant {3}', parent: 'news', amount: 5 },
    { name: 'Variant {4,5}', parent: 'news', amount: 4 },
    { name: 'Variant {6}', parent: 'news', amount: 6 },
    { name: 'Variant {7}', parent: 'news', amount: 8 },
    { name: 'Variant {8}', parent: 'news', amount: 8 },

    { name: 'Variant {9,10}', parent: 'news', amount: 6.5 },
    { name: 'Variant {11}', parent: 'news', amount: 9 },
    { name: 'Variant {12}', parent: 'news', amount: 5 },
  ];


  draw_bubble_chart(data);

function draw_bubble_chart(data){

    const svg = d3.select('.canvas')
        .append('svg')
        .attr('width', 1060)
        .attr('height',800)

    const graph = svg.append('g')
        .attr('transform', 'translate(50,50)');

    const stratify = d3.stratify()
    .id(d=>d.name) // for each object what is the id
    .parentId(d=>d.parent); // for each object which is the parent

    console.log(stratify(data)); 
    console.log(stratify(data).data); //root node
    console.log(stratify(data).children[0]); 

    //gets the root node
    const rootNode = stratify(data);

    //create a value for each data. which property should be used to create a value.
    // we need to sum up the amount. sum appends a value for each node, which will be needed by the visualizationos
    rootNode.sum(d=>d.amount);
    //console.log(rootNode.value);

    const pack = d3.pack()
        .size([960,700])
        .padding(5)

    //appends the dimensions x and y needed for the circle. the vales depends on the size given to the pack generator
    console.log(pack(rootNode));

    //flattens the array
    console.log(pack(rootNode).descendants());

    const bubbleData = pack(rootNode).descendants()

    //create ordinal scales
    //const colour = d3.scaleOrdinal(d3['schemeSet1']);
    //colour.domain(bubbleData.map(d => d.id));

    //use the depth value: 0,1,2 in data
    const colour = d3.scaleOrdinal(['#d1c4e9', '#b39ddb', '#9575cd']); //different shades of purple

    //join data to  group .
    const nodes = graph.selectAll('g')
        .data(bubbleData)
        

    //why var? we need to append circles and text separately within a group and not attach text inside circle    
    const nodesEnter =  nodes.enter()
        .append('g')
        .attr('transform', d =>  `translate(${d.x}, ${d.y})` ) //for each d in data, position each group according to x and y
        
    nodesEnter.append('circle')
        .attr('r', d => d.r)
        .attr('stroke', 'white')
        .attr('stroke-width', 1)
        //.attr('fill', d => colour(d.id))
        .attr('fill', d => colour(d.depth))

    nodesEnter.filter(d=> !d.children)
        .selectAll('circle')
        .attr('class', 'childnode') 

    //console.log(nodes)
    //set name only to noodes that dont have children
    nodesEnter.filter(d=> !d.children) //d.children is False if there are no children   
        .append('text')
        .attr('text-anchor', 'middle')
        .attr('dy', '0.3em') //offset from default y position
        .attr('fill', 'white')
        .style('font-size', '18px') //larger the value, larger the font size
        .text(d => d.data.name);


    //tool tip setup
    const tip = d3.tip()
        .attr('class', 'd3-tip card') // card is materialize class to add style
        .html((e,d) => {
        // let is when a variables scope should be local
        //let content = `<div class="name">{variant 2, variant 3}</div>`;
        let content = ``;
        content += `<div class="cost">Avg. Exec. time: ${d.value} ms </div>`;
        return content
    })

    graph.call(tip);

    d3.select('.d3-tip')
        .style('background', '#333')
        .style('color', '#fff');

    graph.selectAll('.childnode')
    .on('mouseover', (e, d) => {
        tip.show(e,d);
        //console.log("hi")
        handleMouseOver(e,d);
    })
    .on('mouseout', (e,d) => {
        tip.hide();
        handleMouseOut(e,d);
    })
    //.on('click', handleClick)
}

const handleMouseOver = (e,d) => {
    //console.log(e.currentTarget);

    // on hover change color to white
    //name the transitions. otherwise, when the page is refreshed and you hoover over before the chart draw transitiono is complete, the chart will break.
    d3.select(e.currentTarget)
        .transition('changeSliceFill').duration(300) // name the transition so that they do not interfere with other transitioons
            .attr('fill', '#9575cd')

}

const handleMouseOut = (e,d) => {

    //reset the color when mouse hovers out
    d3.select(e.currentTarget)
        .transition('changeSliceFill').duration(100) // same name as as hover as they belong to same type
            .attr('fill', '#b39ddb')

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