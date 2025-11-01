
storeNetwork <- function(plotmat, filename=NULL, 
                         grps=NULL, predset=NULL, 
                         voronoi=FALSE,
                         width=2500, height=2000, pointsize=10,
                         repulsion=1.1, ...){
  args<-list(...)
  
  if(!is.null(grps)){
    plotgrps <- factor(grps$name[match(colnames(plotmat), grps$state)])
    plotcols <- grps$color[match(levels(plotgrps), grps$name)]
  }else{
    plotgrps <- NULL
    plotcols <- "white"
  }
  
  if(!is.null(predset)){
    plotpcts <- predset$pct[match(colnames(plotmat), predset$state)]
    pieBorder <- 0.3
  }else{
    plotpcts <- NULL
    pieBorder <- 0
  }
  
  std.args<-list(input=plotmat, repulsion=repulsion, layout="spring",
                 label.cex=5, label.scale.equal=T, curveDefault=.25,
                 labels=colnames(plotmat),
                 border.color="grey", normalize=T,
                 groups=plotgrps, legend.cex=2, asize=2,
                 color=plotcols, posCol="blue", negCol="red",
                 pie=plotpcts, pieBorder=pieBorder,
                 mar=rep(2,4))
  std.args[names(args)] <- args
  
  if(!is.null(filename)){
    png(filename, width=width, height=height, pointsize=pointsize)
  }
  
  do.call(qgraph, std.args)
  
  if(voronoi > 0){
    qgr <- do.call(qgraph, c(std.args, DoNotPlot=T))
    tesselation <- deldir::deldir(qgr$layout[,1], qgr$layout[,2], sort=F)
    tiles <- deldir::tile.list(tesselation)
    if(!is.null(grps)){
      brighten <- .25
      cols <- plotcols[match(plotgrps, levels(plotgrps))] %>%
        col2rgb() %>% t() 
      cols <- floor(cols*brighten+255-255*brighten) %>%
        as.hexmode() %>% format(width=2, upper.case=T) %>% matrix(ncol=3) %>%
        apply(1, function(x) {paste(x, collapse="")})
      cols <- paste0("#", cols)
    }else{
      cols <- "#FFFFFF"
    }
    if(voronoi > 1){
      s <- seq(0, 2 * pi, length.out = 3000)
      circle <- list(x = 1.1 * (cos(s)),
                     y = 1.1 * (sin(s)))
      plot(tiles, pch = 19, axes=F, xlab="", ylab="", close=T, add=T, 
           showpoints=F, fillcol=cols, clipp=circle)
    }else{
      plot(tiles, pch = 19, axes=F, xlab="", ylab="", close=T, add=T, 
           showpoints=F, fillcol=cols)
    }
    do.call(qgraph,c(std.args, plot=F, legend=F))
  }
  
  if(!is.null(filename)){
    dev.off()
  }
}
