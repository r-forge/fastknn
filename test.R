library(SearchTrees)
dotest = function()
  {
    x = rnorm(100)
    y = rnorm(100)
    dat = cbind(x,y)
    newdat = cbind(c(.3, -.4), c(-.4, .3))
    tree = createIndex(dat)
    inds = findKNN(tree, newdat, dat)
    ds = as.matrix(dist(rbind(newdat, dat)))
    true1 =  order(ds[1,][-(1:2)])[1:5]
    true2 =  order(ds[2,][-(1:2)])[1:5]
    trueknn = rbind(true1, true2)
    toret = all(inds == trueknn)
    if (!toret)
      browser()
    toret
  }
