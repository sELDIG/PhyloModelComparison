modelMetricDivergent <- function(mat, mat2, xCol = c("darkblue","blue", "lightblue", "gray94", "pink","red", "darkred"), 
                         range = c(-1,1), x.labels = dimnames(mat)[[1]],      
                         y.labels = dimnames(mat)[[2]], alpha = .01, cex.metrics = 0.5, 
                         cex.models = 0.7, cex.processes = 1,
                         cex.asterisk = 1,
                         boxes = FALSE, #draw a highlight box around consistent responses
                         boxProcesses, #vector of processes to highlight (of same length as boxMetrics, and processes can be repeated to highlight diff metrics)
                         boxMetrics, #vector of metrics to highlight (of same length as boxProcesses, and metrics can be repeated to highlight diff processes)
                         lwd = 4,
                         srt = 60, 
                         metric.order = 1:dim(mat)[1] #default is to list metrics in the order they appear in the columns of simuData, but can optionally pass an alternative order, e.g., that derived from the metric-correlation clustering
)
{ 
  
  processPosition = data.frame(process = c('env', 'dis', 'nic', 'mut', 'com'),
                               xPos = c(0.6, 10.5, 20.5, 30.5, 40.5))
  
  metricPosition = data.frame(metrics = rownames(mat)[metric.order], 
                              yPos = seq(0.5, nrow(mat)-0.5, by = 1))
  
  newMat = matrix(nrow = dim(mat)[1], ncol = 50)
  row.index = 0
  for(i in 1:dim(mat)[1]){   
    row.index = row.index + 1
    for(k in 1:dim(mat)[2]) {
      newMat[i,(1+(k-1)*10):((k-1)*10 + dim(mat)[3])] = mat[dend.order[row.index],k,]
    }
  }
  
  newMat <- t(newMat)
  xpos = rep(0:(dim(mat)[2]-1), each = dim(mat)[3])*10 + 1:dim(mat)[3] 
  ypos = 1:dim(mat)[1]
  fields::image.plot(x = 1:50, y = 1:dim(mat)[1], z = newMat, axes = FALSE, zlim = range, 
                     col = colorRampPalette(xCol)(30), xlab = "", ylab = "")
  
  abline(v = c(0,10,20,30,40,50) + 0.5)
  abline(v = c(0,10,20,30,40,50) - 1.5)
  abline(h = c(1:(dim(mat)[1] +1)) + 0.5)
  axis(2, at = ypos, labels = x.labels[metric.order], las = 2, cex.axis = cex.metrics)
  axis(3, at = seq(5, 45, by = 10), labels = y.labels, cex.axis = cex.processes, las = 1)
  axis(1, at = xpos, labels = F, tck = -0.01)
  box() 
  
  text(x = xpos, y = 0, srt = srt, adj = 1, xpd = TRUE, labels = rep(dimnames(mat)[[3]], 5), cex = cex.models)
  
  # Add *s where correlation's p is below alpha
  for(i in 1:dim(mat2)[1]){     
    for(k in 1:dim(mat2)[2]){
      for(j in 1:dim(mat2)[3])
        if(!is.na(mat2[i,k,j]) & mat2[i,k,j] < alpha){
          text(x = 10*(k-1) +j,  y = (1:dim(mat2)[1])[metricPosition$metrics == x.labels[i]], labels = "*", cex = cex.asterisk )
        }
    }
  } # end * loop
  
  if(boxes) {
    for (b in 1:length(boxProcesses)) {
      rect(processPosition$xPos[processPosition$process == boxProcesses[b]],
           metricPosition$yPos[metricPosition$metric == boxMetrics[b]],
           processPosition$xPos[processPosition$process == boxProcesses[b]] + 8,
           metricPosition$yPos[metricPosition$metric == boxMetrics[b]] + .95,
           lwd = lwd)
      
    } # end box loop
  } # end box if
  
} #end function





modelMetricContinous <- function(mat, xCol = c( "gray94", "lightblue", "blue", "darkblue"), 
                             range = c(0,1), x.labels = dimnames(mat)[[1]],      
                             y.labels = dimnames(mat)[[2]], alpha = .01, cex.metrics = 0.5, 
                             cex.models = 0.7, cex.processes = 1,
                             lwd = 4,
                             srt = 60, 
                             metric.order = 1:dim(mat)[1] #default is to list metrics in the order they appear in the columns of simuData, but can optionally pass an alternative order, e.g., that derived from the metric-correlation clustering
)
{ 
  
  processPosition = data.frame(process = c('env', 'dis', 'nic', 'mut', 'com'),
                               xPos = c(0.6, 10.5, 20.5, 30.5, 40.5))
  
  newMat = matrix(nrow = dim(mat)[1], ncol = 50)
  row.index = 0
  for(i in 1:dim(mat)[1]){   
    row.index = row.index + 1
    for(k in 1:dim(mat)[2]) {
      newMat[i,(1+(k-1)*10):((k-1)*10 + dim(mat)[3])] = mat[dend.order[row.index],k,]
    }
  }
  
  newMat <- t(newMat)
  xpos = rep(0:(dim(mat)[2]-1), each = dim(mat)[3])*10 + 1:dim(mat)[3] 
  ypos = 1:dim(mat)[1]
  fields::image.plot(x = 1:50, y = 1:dim(mat)[1], z = newMat, axes = FALSE, 
                     col = colorRampPalette(xCol)(30), xlab = "", ylab = "")
  
  abline(v = c(0,10,20,30,40,50) + 0.5)
  abline(v = c(0,10,20,30,40,50) - 1.5)
  abline(h = c(1:(dim(mat)[1] +1)) + 0.5)
  axis(2, at = ypos, labels = x.labels[metric.order], las = 2, cex.axis = cex.metrics)
  axis(3, at = seq(5, 45, by = 10), labels = y.labels, cex.axis = cex.processes, las = 1)
  axis(1, at = xpos, labels = F, tck = -0.01)
  box() 
  
  text(x = xpos, y = 0, srt = srt, adj = 1, xpd = TRUE, labels = rep(dimnames(mat)[[3]], 5), cex = cex.models)
  
  
} #end function



image.real <- function(mat, xCol = c("blue", "white", "white", "red"), 
                       range = c(-1,1), x.labels = rownames(mat), 
                       y.labels = colnames(mat), cex.axis = 1) { 
  mat <- t(mat)[,nrow(mat):1]
  fields::image.plot(mat, axes = FALSE, zlim = range, 
                     col = colorRampPalette(xCol)(30))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = x.labels, cex.axis = cex.axis)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = y.labels, las = 2, cex.axis = cex.axis)
  box() 
}


crossPredictabilityPlot = function(CpredictabilityR2,
                                   cexValues = 1, cexMain = 2, 
                                   cexAxis = 1.5, cexLab = 1.8,
                                   mainText, cexLabSub = 1.2,
                                   range = c(-1,1),
                                   xCol = c("darkblue","blue",
                                            "lightblue", "white",
                                            "pink","red", "darkred"))
{
  
  image.real(CpredictabilityR2, range = range, xCol = xCol, cex.axis = cexAxis)
  title(main = mainText, xlab = "predicted to trees from", 
        cex.lab = cexLab, cex.main = cexMain)
  mtext("Model trained on\n", 2, line = 2.7, cex = cexLab)
  mtext("(modelers' assumption of reality)", 2, line = 3.4, cex = cexLabSub)
  mtext("(what if reality works like this model?)", 1, line = 4.5, cex = cexLabSub)
  
  # add values to cells
  for(i in 1:nrow(CpredictabilityR2)){
    for(j in 1:nrow(CpredictabilityR2)){
      text((j-1)/(nrow(CpredictabilityR2) - 1),
           1-(i-1)/(nrow(CpredictabilityR2) - 1), 
           labels = format(round(CpredictabilityR2[i,j], digits = 2), nsmall = 2),
           cex = cexValues)
      
    }
  }
}


