library(SearchTrees)
dotest = function()
  {
    x = rnorm(100)
    y = rnorm(100)
    dat = cbind(x,y)
    newdat = cbind(c(.3, -.4), c(-.4, .3))
    tree = createTree(dat)
    inds = knnLookup(tree, newdat =  newdat)
    ds = as.matrix(dist(rbind(newdat, dat)))
    true1 =  order(ds[1,][-(1:2)])[1:5]
    true2 =  order(ds[2,][-(1:2)])[1:5]
    trueknn = rbind(true1, true2)
    toret = all(inds == trueknn)
    if (!toret)
      browser()
    toret
  }



dat = cbind(x,y)
tree = createTree(dat)
thing = rectLookup(tree, c(0,0), c(1,1))

plot(x,y)
abline(h=c(0, 1), v=c(0,1))
points(dat[thing,], pch="*", col = "blue")


x2 = x + runif(100, .5, 2)
y2 = y + runif(100, .5, 2)
dat2 = cbind(x, y, x2, y2)
tree2 = createTree(dat2, dataType="rect", columns= 1:4)


thing2 = rectLookup(tree2, xlim = c(0,1), ylim=c(0, 1))


