var ContinuousInteractionHandler = function(width, height, spaceWidth, spaceHeight, ctx){

    // Find scaled size:
    const scaledWidth = Math.floor(width / spaceWidth);
    const scaledHeight = Math.floor(height / spaceHeight);
  
    const lineHeight = 10;
  
      // list of standard rendering features to ignore (and key-values in the portrayal will be added )
    const ignoredFeatures = [
        'Shape',
        'Filled',
        'Color',
        'r',
        'w',
        'h',
        'width',
        'height',
        'heading_x',
        'heading_y',
        'stroke_color',
        'text_color'
    ];
  
    // Set a variable to hold the lookup table and make it accessible to draw scripts
    var mouseoverLookupTable = this.mouseoverLookupTable = buildLookupTable(spaceWidth, spaceHeight);
    function buildLookupTable(spaceWidth, spaceHeight){
      var lookupTable;
      this.init = function(){
         lookupTable = {};
         lookupTable["spaceWidth"] = spaceWidth;
         lookupTable["spaceHeight"] = spaceHeight;
      }
  
      this.set = function(x, y, value){
              let point = "(" + x + "," + y + ")";
              if (point in Object.keys(lookupTable)){
                      lookupTable[point].push(value);
              }
              else{
                      lookupTable[point] = [value];
                }
      }
  
      this.get = function(x, y){
        let point = "(" + x + "," + y + ")";
        if (point in Object.keys(lookupTable)){
                return lookupTable[point];
        }
        return [];
      }
  
      return this;
    }
  
    var coordinateMapper = function(event){
      return {
        x: Math.floor(event.offsetX/scaledWidth),
        y: Math.floor(event.offsetY/scaledHeight)
      };
    };
  /*  this.setCoordinateMapper = function(){  
       default coordinate mapper for grids

    };
  
    this.setCoordinateMapper();*/
  
  
    // wrap the rect styling in a function
    function drawTooltipBox(ctx, x, y, width, height){
      ctx.fillStyle = "#F0F0F0";
      ctx.beginPath();
      ctx.shadowOffsetX = -3;
      ctx.shadowOffsetY = 2;
      ctx.shadowBlur = 6;
      ctx.shadowColor = "#33333377";
      ctx.rect(x, y, width, height);
      ctx.fill();
      ctx.shadowColor = "transparent";
    }
  
    var listener; var tmp
    this.updateMouseListeners = function(portrayalLayer){tmp = portrayalLayer
  
        // Remove the prior event listener to avoid creating a new one every step
        ctx.canvas.removeEventListener("mousemove", listener);
  
        // define the event litser for this step
        listener = function(event){
                // clear the previous interaction
          ctx.clearRect(0, 0, width, height);
          // map the event to x,y coordinates
          const position = coordinateMapper(event);
         // const yPosition = Math.floor(event.offsetY/scaledHeight);
         //       const xPosition = Math.floor(event.offsetX/scaledWidth);
  
          // look up the portrayal items the coordinates refer to and draw a tooltip
          mouseoverLookupTable.get(position.x, position.y).forEach((portrayalIndex, nthAgent) => {
            const agent = portrayalLayer[portrayalIndex];
            const features = Object.keys(agent).filter(k => ignoredFeatures.indexOf(k) < 0);
            const textWidth = Math.max.apply(null, features.map(k => ctx.measureText(`${k}: ${agent[k]}`).width));
            const textHeight = features.length * lineHeight;
            const y = Math.max(lineHeight * 2, Math.min(height - textHeight, event.offsetY - textHeight/2));
            const rectMargin = 2 * lineHeight;
            var x = 0;
            var rectX = 0;
            if(event.offsetX < width/2){
              x = event.offsetX + rectMargin + nthAgent * (textWidth + rectMargin);
              ctx.textAlign = "left";
              rectX = x - rectMargin/2;
            } else {
              x = event.offsetX - rectMargin - nthAgent * (textWidth + rectMargin + lineHeight );
              ctx.textAlign = "right";
              rectX = x - textWidth - rectMargin/2;
            }
  
              // draw a background box
            drawTooltipBox(ctx, rectX, y - rectMargin, textWidth + rectMargin, textHeight + rectMargin);
  
              // set the color and draw the text
            ctx.fillStyle = "black";
            features.forEach((k,i) => {
              ctx.fillText(`${k}: ${agent[k]}`, x, y + i * lineHeight)
            })
          })
  
        };
      ctx.canvas.addEventListener("mousemove", listener);
    };
  
    return this;
  }