
! ******************************************************************
! This is the gsw Absolute Salinity algorithm implemented in Fortran
! ******************************************************************
!
! as defined in 
!
!   McDougall, T.J., Jackett, D.R. and Millero, F.J., 2009: An algorithm for estimating  
!                       Absolute Salinity in the global ocean, Ocean Science, submitted.  
!     
!   
! David Jackett
! January 2009
!
! Karline Soetaert - May 2009: made subroutines instead of functions
! isnan -> -99
!
!******************************************************************************

!******************************************************************************

        SUBROUTINE gsw_delta_sa(p0,longs0,lats0,longs,lats,p,ndepth,           &
     &                         del_sa,delta)

! calculate the Absolute Salinity anomaly
!
! SP                  : Practical Salinity                 [psu]
! p0                  : sea (gauge) pressure               [dbar]
! longs0              : longitude                          [deg E]     
! lats0               : latitude                           [deg N]
!
! result              : Absolute Salinity anomaly          [g/kg]

        implicit none

        integer nx, ny, nz
        parameter(nx=91,ny=44,nz=45)

        integer indx0, indy, indz, i, icalled
        integer k, deli(4), delj(4), nmean
        real*8 p0, p0_original, longs0, lats0, sa_upper
        real*8 sa_lower, delta
        real*8 longs(nx), lats(ny), p(nz), del_sa(nz,ny,nx)
	      real*8 dlongs, dlats
        real*8 r1, s1, t1, dsa_mean, dsa(4), ndepth(ny,nx), ndepth_max
        logical ISNAN

        data deli/0,1,1,0/, delj/0,0,1,1/

        dlongs = longs(2)-longs(1)
		    dlats = lats(2)-lats(1)

        indx0 = floor(1 + (nx-1)*(longs0-longs(1))/(longs(nx)-longs(1)))
        if(indx0.eq.nx) then
         indx0 = nx-1
        end if

        indy = floor(1 + (ny-1)*(lats0-lats(1))/(lats(ny)-lats(1)))
        if(indy.eq.ny) then
         indy = ny-1
        end if

        ndepth_max = -1.d0
        do k = 1,4
         if(ndepth(indy+delj(k),indx0+deli(k)).gt.0.d0)                              &
     &    ndepth_max =max(ndepth_max,ndepth(indy+delj(k),indx0+deli(k)))
        end do

        if(ndepth_max.eq.-1.d0) then
         delta = 0.d0 !/0.
         return
        end if    

        p0_original = p0

        if(p0.gt.p(int(ndepth_max))) p0 = p(int(ndepth_max))

        call indx(p,nz,p0,indz)
    
        r1 = (longs0-longs(indx0))/(longs(indx0+1)-longs(indx0))
        s1 = (lats0-lats(indy))/(lats(indy+1)-lats(indy))
        t1 = (p0-p(indz))/(p(indz+1)-p(indz))

        isnan = .FALSE.
        do k = 1,4
            dsa(k) = del_sa(indz,indy+delj(k),indx0+deli(k))
	      if (dsa(k) .LE. -90) isnan=.TRUE. 
        end do

        if(260.d0.le.longs0.and.longs0.le.291.999d0.and.                       &
     &     3.4d0.le.lats0.and.lats0.le.19.55d0) then
           CALL dsa_add_barrier(dsa,longs0,lats0,longs(indx0),                &
     &           lats(indy),dlongs,dlats)

        elseif(isnan) then 
          call dsa_add_mean(dsa,longs0,lats0)
        end if

        sa_upper = (1.d0-s1)*(dsa(1) + r1*(dsa(2)-dsa(1))) + s1*(dsa(4)        &
     &      +  r1*(dsa(3)-dsa(4)))

        isnan = .FALSE.
        do k = 1,4
            dsa(k) = del_sa(indz+1,indy+delj(k),indx0+deli(k))
	      if (dsa(k) .LE. -90) isnan=.TRUE.
        end do

        if(260.d0.le.longs0.and.longs0.le.291.999d0.and.3.4d0.le.lats0         &
     &      .and.lats0.le.19.55d0) then
            CALL dsa_add_barrier(dsa,longs0,lats0,longs(indx0),               &
     &     lats(indy),dlongs,dlats)
        elseif(isnan) then 
            isnan = .FALSE.
            call dsa_add_mean(dsa,longs0,lats0)
            do k=1,4
		    if (dsa(k) .LE. -90) isnan=.TRUE.
            enddo 
        end if

        sa_lower = (1.d0-s1)*(dsa(1) + r1*(dsa(2)-dsa(1))) + s1*(dsa(4)        &
     &    + r1*(dsa(3)-dsa(4)))

        if(isnan) sa_lower = sa_upper
                            
        delta = sa_upper + t1*(sa_lower-sa_upper)

!        if(isnan(delta)) delta = 0.d0

!        p0 = p0_original
  
        return
        end

!******************************************************************************

        subroutine adjust_baltic(SA,SP,longs,lats,flag,S)

! for the Baltic Sea, overwrite Absolute Salinity with a value
! computed analytically from Practical Salinity, or vice versa
!
! SA                  : Absolute Salinity                  [g/kg]
! SP                  : Practical Salinity                 [psu]
! longs               : longitude                          [deg E]     
! lats                : latitude                           [deg N]
! flag                : flag - 1 or 2 
!
! adjust_Baltic       : Absolute Salinity                  [g/kg]
!                         when flag = 1
!                     : Practical Salinity                 [psu]
!                         when flag = 2

       implicit none
       integer n2, n3, flag
       real*8 SA, SP, S, longs, lats, xinterp1
       real*8 xb_left(3), xb_right(2), yb_left(3), yb_right(2)                 
       real*8 xx_left, xx_right

       data xb_left/12.6d0, 7.d0, 26.d0/, yb_left/50.d0, 59.d0, 69.d0/
       data xb_right/45.d0, 26.d0/, yb_right/50.d0, 69.d0/

       external xinterp1

       n2 = 2
       n3 = 3
       if(flag.eq.1) then
         S = SA
       elseif(flag.eq.2) then
         S = SP
       end if

       if(xb_left(2).lt.longs .and. longs.lt.xb_right(1) .and.                 &
     &    yb_left(1).lt.lats .and. lats.lt.yb_left(3)) then
  
         xx_left = xinterp1(yb_left, xb_left, n3, lats)
    
         xx_right = xinterp1(yb_right, xb_right, n2, lats)
    
         if(xx_left.le.longs .and. longs.le.xx_right) then
           if(flag.eq.1) then
              S = (35.16504d0/35.d0)*SP + 0.124d0*(1.d0-SP/35.d0)
           elseif(flag.eq.2) then
              S = (35.d0/35.04104d0)*(SA- 0.124d0)
           end if
         end if
     
       end if

       return
       end

!******************************************************************************

       subroutine dsa_add_barrier(dsa,longs0,lats0,longs,lats,dlongs,dlats)

! add a barrier through Central America (Panama) and then average
! over the appropriate side of the barrier
! 
! dsa                 : absolute salinity anomaly                [g/kg]
! longs0              : longitudes of data                       [deg E]
! lats0               : latitudes of data                        [deg N]
! longs               : longitudes of regular grid               [deg E]
! lats                : latitudes of regular grid                [deg N]
! dlongs              : longitude difference of regular grid     [deg longitude]
! dlats               : latitudes difference of regular grid     [deg latitude]
!
! result              : absolute salinity anomaly of data        [g/kg]

       implicit none

       integer k, nmean, above_line0, above_line(4)

       real*8 dsa_add(4), dsa(4), dsa_mean
       real*8 longs_pan(6), lats_pan(6), longs0, lats0, r, lats_line
       real*8 longs, lats, dlongs, dlats

       data longs_pan/260.0000, 272.5900, 276.5000, 278.6500,                  &
     &      280.7300, 292.000/ 
       data  lats_pan/ 19.5500,  13.9700,   9.6000,   8.1000,                  &
     &      9.3300,   3.400/ 


        call indx(longs_pan,6,longs0,k)                            !   the long0/lat0 point
        r = (longs0-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
        lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))
        if(lats_line.le.lats0) then
            above_line0 = 1
        else
            above_line0 = 0
        end if

!print *, 'above_line = ', above_line0

        call indx(longs_pan,6,longs,k)                                     !  the 1 and 4 long/lat points 
        r = (longs-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
        lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))
        if(lats_line.le.lats) then
            above_line(1) = 1
        else
            above_line(1) = 0
        end if
        if(lats_line.le.lats+dlats) then
            above_line(4) = 1
        else
            above_line(4) = 0
        end if

        call indx(longs_pan,6,longs+dlongs,k)                              !  the 2 and 3 long/lat points 
        r = (longs+dlongs-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
        lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))
        if(lats_line.le.lats) then
            above_line(2) = 1
        else
            above_line(2) = 0
        end if
        if(lats_line.le.lats+dlats) then
            above_line(3) = 1
        else
            above_line(3) = 0
        end if

        nmean = 0
        dsa_mean = 0.d0
        dsa_add = dsa
        do k = 1,4
            if (dsa(k).GT.-90.and.above_line0.eq.above_line(k)) then
                nmean = nmean+1 
				dsa_mean = dsa_mean+dsa(k)
            end if
        end do

        dsa_mean = dsa_mean/nmean

        do k = 1,4
            if(dsa(k).LE. -90 .or.above_line0.ne.above_line(k)) then
                 dsa_add(k) = dsa_mean
            end if
        end do

        dsa = dsa_add
        return
        end

!******************************************************************************

        subroutine dsa_add_mean(dsa,longs0,lats0)

! replace nans with namean   
!
! dsa           : absolute salinity anomaly                [g/kg]
!                 of the 4 adjacent neighbours     
!
! result        : nanmean of the 4 adjacent neighbours     [g/kg]

        implicit none

        integer k, nmean
        real*8 dsa(4)
        real*8 dsa_mean, dsa_add(4)
        real*8 longs0, lats0

        nmean = 0
        dsa_mean = 0.d0
        dsa_add = dsa

        do k = 1,4
            if (dsa(k).GT. -90) then
                nmean = nmean+1 
				dsa_mean = dsa_mean+dsa(k)
            end if
        end do

        dsa_mean = dsa_mean/nmean

        do k = 1,4
            if(dsa(k).LE. -90) then
                 dsa_add(k) = dsa_mean
            end if
        end do

        dsa = dsa_mean
        return
        end

!******************************************************************************

        subroutine indx(x,n,z,k)

!   DESCRIPTION:    Find the index of a real number in a
!                   monotonically increasing real array
!
!   INPUT:          x       array of increasing values
!                   n       length of array
!                   z       real number
!
!   OUTPUT:         k       index k - if x(k) <= z < x(k+1), or
!                   n-1             - if z = x(n)
!
!   CREATED:        June 1993
!

       implicit none
	   character(len=80) msg

       integer n,k,ku,km,kl
	   real*8 x(n),z

       if(x(1).lt.z.and.z.lt.x(n)) then

         kl=1
         ku=n

         do while (ku-kl.gt.1)
           km=(ku+kl)/2
           if(z.gt.x(km))then
             kl=km
           else
             ku=km
           endif
         end do

         k=kl

         if(z.eq.x(k+1)) k = k+1

       else

        if(z.eq.x(1)) then
          k = 1
        elseif(z.eq.x(n)) then
          k = n-1
        else
          call rwarn("ERROR in indx.f : out of range")
          write (msg,*) z,n,x
 	    call rexit(msg)
        end if

      end if

      return
      end

!******************************************************************************

        function xinterp1(x,y,n,x0)

! linearly interpolate a real array   
!
! x                   : x array (monotonic)               
! y                   : y array     
! n                   : length of x and y arrays
! x0                  : x value to be interpolated
!
! result              : interpolated value

        implicit none
        integer n, k
        real*8 x(n), y(n), x0, r, xinterp1

        call indx(x,n,x0,k)

        r = (x0-x(k))/(x(k+1)-x(k))

        xinterp1 = y(k) + r*(y(k+1)-y(k))

        return
        end

!******************************************************************************