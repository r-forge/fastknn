setClass("KDTree", representation = list(ref = "externalptr", k = "integer", maxDepth = "integer", npoints = "integer"));

kdTree = function(data)
  {
    if(is.null(dim(data)))
      stop(paste("Data does not appear to be of the right class. is of class", class(data)))

    storage.mode(data) <- "double"
    ords = apply(data, 2, order)
    k = ncol(data)
    n = nrow(data) 
    res = .Call("R_Build_kdTree", data, ords, k, n)
    res@npoints = n
    res@k = k
    res
  }

getKDTestOrder = function(tree)
  {
    .Call("R_Get_Left_Order", tree)
  }
