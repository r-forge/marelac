
# Reads the original file from the code by David Jackett
# nx=91,ny=44,nz=45
# longs, lats: latitiudes, longitudes for which salinity correction is known
# longs[91], lats[44]
# ndepth: depths of salinity correction ndepth[long,lat]
# del_sa:salinity correction del_sal[nz,ny,nx]

dsal <- read.table("gsw_data.dat")
dsal[dsal=="NaN"]<- -99

sw_sfac <- list()
sw_sfac$longs <- dsal[1:91,1]
sw_sfac$lats  <- dsal[92:(92+43),1]
is <- 92+44
sw_sfac$p  <- dsal[is:(is+44),]
is <- 92+44+45
sw_sfac$ndepth <- matrix(nc=91,nr=44,dsal[is:(is+(91*44)-1),1])
ii <- is+(91*44)
sw_sfac$del_sa <- array(dim=c(45,44,91),data=dsal[ii:(ii+45*44*91-1),1])
rm(dsal,is,ii)
