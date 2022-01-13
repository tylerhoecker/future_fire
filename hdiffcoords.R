

h1 <- eg_hex(ref_df$ndvi, ref_df$cbi)[['hex']]
h2 <- eg_hex(fut_df$ndvi, fut_df$cbi)[['hex']]

hdiffcoords <- function(bin1, bin2) {
  
  ## Checks modified from: https://github.com/edzer/hexbin/blob/master/R/hdiffplot.R
  if(is.null(bin1) | is.null(bin1)) {
    stop("Need 2 hex bin objects")
  } else {
    if(bin1@shape != bin2@shape)
      stop("Bin objects must have same shape parameter")
    if(all(bin1@xbnds == bin2@xbnds) & all(bin1@ybnds == bin2@ybnds))
      equal.bounds <- TRUE
    else stop("Bin objects need the same xbnds and ybnds")
    if(bin1@xbins != bin2@xbins)
      stop("Bin objects need the same number of bins")
  }
  
  ## Find overlapping and unique bins
  
  hd1 <- data.frame(hcell2xy(bin1), count_bin1 = bin1@count, cell_bin1 = bin1@cell)
  hd2 <- data.frame(hcell2xy(bin2), count_bin2 = bin2@count, cell_bin2 = bin2@cell)
  
  overlapping_hd1 <- apply(hd1, 1, function(r, A){ sum(A$x==r[1] & A$y==r[2]) }, hd2)
  overlapping_hd2 <- apply(hd2, 1, function(r, A){ sum(A$x==r[1] & A$y==r[2]) }, hd1)
  
  overlaps <- merge(hd1[as.logical(overlapping_hd1),], hd2[as.logical(overlapping_hd2),])
  
  unique_hd1 <- hd1[!as.logical(overlapping_hd1),]
  unique_hd2 <- hd2[!as.logical(overlapping_hd2),]
  
  ## Return list of data.frames
  
  list(unique_bin1 = unique_hd1, unique_bin2 = unique_hd2, overlapping = overlaps)
  
}


hdiff_df <- hdiffcoords(h1,h2)

difference <- 
  ggplot() + 
  geom_point(data = hdiff_df$unique_bin1, aes(x = x, y = y), color = "red", size = 10) + 
  geom_point(data = hdiff_df$unique_bin2, aes(x = x, y = y), color = "cyan", size = 10) +
  geom_point(data = hdiff_df$overlapping, aes(x = x, y = y, color = count_bin2-count_bin1), size = 10) + theme_bw()




hdiffplot(bin1 = h1, bin2 = h2, xbnds = c(0,1), ybnds = c(0,3),
          dif.cols = 'yellow')


