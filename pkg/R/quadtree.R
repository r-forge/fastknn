setClass("SearchTree", representation = list(ref = "externalptr", numNodes = "integer", dataNodes = "integer", maxDepth = "integer", maxBucket = "integer", points = "integer", "VIRTUAL"))
setClass("QuadTree", contains="SearchTree")

setGeneric("getPointsInRect",
           function(tree, ptOne, ptTwo, data, columns = 1:2)
           standardGeneric("getPointsInRect")
           )
setMethod("getPointsInRect", "QuadTree",
          function(tree, ptOne, ptTwo, data, columns)
          {
            left = min(ptOne[1], ptTwo[1])
            right = max(ptOne[1], ptTwo[1])
            down = min(ptOne[2], ptTwo[2])
            up = max(ptOne[2], ptTwo[2])
            if(is(data, "matrix"))
               mode(data) = "numeric"
            else
              {
                for (i in columns)
                  data[[i]] = as.numeric(data[[i]])
              }
            getPointsInBox(tree, left, right, down, up, data=data, cols = columns)
          }
          )

setGeneric("findKNN",
           function(tree, newdat, fulldat, newcols=1:2, fullcols = 1:2,  k = 5)
           standardGeneric("findKNN")
           )
setMethod("findKNN", "QuadTree",
          function(tree, newdat, fulldat, newcols, fullcols, k)
          {
            if (is(newdat, "matrix"))
              {
                newtype = 1L
                if (mode(newdat) != "double")
                  mode( newdat ) = "double"
              } else if (is(newdat, "data.frame")) {
                newtype = 2L
                
              } else {
                stop("newdat must be either a matrix or a data.frame")
              }
           if (is(fulldat, "matrix"))
              {
                fulltype = 1L
                if (mode(fulldat) != "double")
                  mode(fulldat) = "double"
              } else if (is(fulldat, "data.frame")) {
                fulltype = 2L
                
              } else {
                stop("fulldat must be either a matrix or a data.frame")
              }
            newcols = as.integer(newcols)
            fullcols = as.integer(fullcols)
            if (length(newcols) != 2)
              stop("Incorrect number of columns specified for new data")
            if (length(fullcols) != 2 )
              stop("Incorrect number of columns specified for full data")
            k = as.integer(k)
            
            inds = .Call("R_Find_KNN", tree, newdat, fulldat, newcols, fullcols, newtype, fulltype, k, nrow(newdat), nrow(fulldat))
            matrix(inds, byrow = TRUE, ncol = k)
          }
          )

createIndex = function(data, type = "quad", columns = c(1, 2), ...)
  {
    if (tolower(type) == "quad")
      {
        if(length(columns) != 2)
          stop("wrong number of columns for this index type.")
        x = data[,columns[1]] 
        y = data[,columns[2]]
        quadTree(x,y, ... )
      }   
  }
        
quadTree = function(x, y, maxDepth = 7, minNodeArea, ...)
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
    .Call("R_Build_Quadtree", x, y, max(x), min(x), max(y), min(y), as.integer(maxDepth))

  }

getTestOrder = function(tree, len)
  {
    .Call("R_Get_Top_Left_Order", tree, as.integer(len))
  }

getMaxDepth = function(tree)
  { 
    .Call("R_Find_Max_Depth", tree)
  }

getPointsInBox =  function(tree, left, right, down, up, x,y, data = NULL, cols = 1:2)
          {
            cols = as.integer(cols)
            if (!is.null(data))
              {
                if(is(data, "matrix"))
                  dattype = 1L
                else if (is(data, "data.frame"))
                  dattype = 2L
                else
                  stop(paste("Unrecognized data class:", class(data)))
                datframe = TRUE
              } else {
                dattype = 0L
                datframe = FALSE
                
                if (length(x) != length(y))
                  stop("Data vectors of unequal length detected.")
                if (length(x) != tree@points)
                  stop("Length of data vectors does not match number of points in tree.")
               }
            if(datframe)
              .Call("R_Get_Points_In_Box", tree, left, right, down, up, data, nrow(data), dattype, cols)
            else
              .Call("R_Get_Points_In_Box", tree, left, right, down, up, x, y, dattype, cols)
          }
          
