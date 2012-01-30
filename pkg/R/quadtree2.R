setClass("SearchTree2", representation = list(ref = "externalptr", numNodes = "integer", dataNodes = "integer", maxDepth = "integer", maxBucket = "integer", totalData = "integer", "VIRTUAL"))
setClass("QuadTree2", contains="SearchTree2")

setGeneric("rectLookup",
           function(tree, ptOne, ptTwo)
           standardGeneric("rectLookup")
           )

setMethod("rectLookup", "QuadTree2",
          function(tree, ptOne, ptTwo)
          {
            ptOne = as.numeric(sort(ptOne))
            ptTwo = as.numeric(sort(ptTwo))
            .Call("R_Rectangle_Lookup", tree, ptOne, ptTwo)
          }
          )


createTree = function(data, type = "quad", columns = c(1, 2), ...)
  {
    if (tolower(type) == "quad")
      {
        if(length(columns) != 2)
          stop("wrong number of columns for this index type.")
        x = data[,columns[1]] 
        y = data[,columns[2]]
        quadTree2(x,y, ... )
      }   
  }
        
quadTree2 = function(x, y, maxDepth = 7, minNodeArea, ...)
  { 
    xrange = range(x)
    yrange = range(y)
    
    if (!missing(minNodeArea))
      {
        totArea = (xrange[2] - xrange[1]) * (yrange[2] - yrange[1])
        areas = which(totArea / (4^(1:10)) <= minNodeArea)
        
        if(!length(areas))
          {
            warning("The minNodeArea selected lead to a maximum depth > 10, which is very memory intensive for negligable benefit. Using maximum depth of 10.")
           
            maxDepth = 10
          } else {
            maxDepth = areas[1]
          }

      }
    x = as.numeric(x)
    y = as.numeric(y)
    .Call("R_Build_Quadtree_Pt", x, y, max(x), min(x), max(y), min(y), as.integer(maxDepth))

  }
