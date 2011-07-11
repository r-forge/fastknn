setClass("SearchTree", representation = list(ref = "externalptr", numNodes = "integer", dataNodes = "integer", maxDepth = "integer", maxBucketSize = "integer", points = "integer"))
setClass("QuadTree", contains="SearchTree")

setGeneric("getPointsInRect",
           function(tree, left, right, down, up, data)
           standardGeneric("getPointsInRect")
           )
setMethod("getPointsInRect", "QuadTree",
          function(tree, left, right, down, up, data)
          {
            getPointsInBox(tree, left, right, down, up, data=data)
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
            
            .Call("R_Find_KNN", tree, newdat, fulldat, newcols, fullcols, newtype, fulltype, k, nrow(newdat), nrow(fulldat))
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

getKNNIndices = function(tree, newx, newy = NULL, allx, ally, k = 5)
  {
    #stop("KNN has not been modified to deal with the new buckets yet!")
    if (length(allx) != length(ally))
    stop("Old data vectors of unequal length detected.")
  if (length(allx) != tree@points)
    stop("Length of old data vectors does not match number of points in tree.")
    if (is.null(newy))
      {
        if (dim(newx) && dim(newx)[2] > 1)
          {
            newy = newx[,2]
            newx = newx[,1]
          } else {
            stop("Not able to find new y values.")
          }
      }
        
    res = .Call("R_Find_KNN", tree, newx, newy, allx, ally, as.integer(k));
    if (length(newx) > 1)
      matrix(res, ncol = k, byrow = TRUE)
    else
      res
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
          
