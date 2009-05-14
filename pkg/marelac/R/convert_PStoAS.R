## -----------------------------------------------------------------------------
## Convert practical salinity to absolute salinity and vice versa
## -----------------------------------------------------------------------------

# internal functions

sw_delta_SA <- function (p0=0,longs0=0,lats0=0) {

  RES <-.Fortran("gsw_delta_sa",as.double(p0),as.double(longs0),as.double(lats0),
           sw_sfac$longs,sw_sfac$lats,sw_sfac$p,sw_sfac$ndepth,sw_sfac$del_sa,
           delta=as.double(0.))
  RES$delta
}

sw_adjust_baltic <- function (SA,PS,long,lat,flag) {
  S <- 1
  RES <- .Fortran("adjust_baltic",as.double(SA), as.double(PS),
     as.double(long),as.double(lat),as.integer(flag), S=as.double(S))
  RES$S
}

# Practical -> absolute salinity
convert_PStoAS <- function(S=35, p=max(0,P-1.013253), P=1.013253,
                           lat=NULL, lon=NULL) {

  SA <- 35.16504 /35.0 *S
  if (! (is.null(lat)) | ! (is.null(lon))) {
    if (is.null(lat))
      stop (" latitude should be given if longitude is")
    if (is.null(lon))
      stop (" longitude should be given if latitude is")
    if (lat > 90 | lat < -90)
      stop (" latitude, 'lat' should be within -90,+90")

    if (lon > 360 | lon < 0)
      stop (" longitude, 'lon' should be within 0,360")

    SA <- SA + sw_delta_SA(p*10,lon,lat)

    flag <- 1

    SA <- sw_adjust_baltic(SA, S, lon, lat, flag)
  }
  SA
}

# Absolute -> practical salinity

convert_AStoPS <- function(S=35, p=max(0,P-1.013253), P=1.013253,
                           lat=NULL, lon=NULL) {
  PS <- 35.0/35.16504 *S
  if (! (is.null(lat)) | ! (is.null(lon))) {
    if (is.null(lat))
      stop (" latitude should be given if longitude is")
    if (is.null(lon))
      stop (" longitude should be given if latitude is")
    if (lat > 90 | lat < -90)
      stop (" latitude, 'lat' should be within -90,+90")

    if (lon > 360 | lon < 0)
      stop (" longitude, 'lon' should be within 0,360")

    PS <- PS - 35.0/35.16504 *sw_delta_SA(p*10,lon,lat)
    flag <- 2

    PS <- sw_adjust_baltic(S,PS,lon,lat,flag)
  }
  PS
}

