cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pinput2( maxnlay,maxnzone,maxnr,
     &                     re,ratc,ratl,
     &tlen,np,omegai,imin,imax,
     &                     nlayer,
     &nzone,vrmin,vrmax,rho,
     &                     vpv,vph,vsv,vsh,eta,qmu,qkappa,
     &r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Parameter Input
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer maxnlay,maxnzone,maxnr
        integer np
        integer imin,imax,nlayer(*)
        integer nzone,nr
        real*8 tlen,omegai,re,ratc,ratl
        real*8 vrmin(*),vrmax(*),rho(4,*)
        real*8 vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
        real*8 qmu(*),qkappa(*)
        real*8 r0,mt(3,3),theta(*),phi(*),lat(*),lon(*)
        real*8 eqlat,eqlon,stlat,stlon,eqlattmp
        character*80 output(*)
        integer i
        character*80 dummy,tmpfile
c
        data tmpfile / 'work' /
c
c temporary file open
        open( unit=11, file=tmpfile, status='unknown' )
c writing to the temporary file
  100   continue
        read(5,110) dummy
  110   format(a80)
        if ( dummy(1:1).eq.'c' ) goto 100
        if ( dummy(1:3).eq.'end' ) goto 120
        write(11,110) dummy
        goto 100
  120   continue
c temporary file close
        close(11)
c
c temporary file open
        open( unit=11, file=tmpfile, status='unknown' )
c reading the parameter
        read(11,*) tlen,np
        read(11,*) re  ! relative error (vertical grid)
        read(11,*) ratc  ! ampratio (vertical grid cut-off)
        read(11,*) ratl  ! ampratio (for l-cutoff)
        read(11,*) omegai  ! artificial damping
        omegai = - dlog(omegai) / tlen
        read(11,*) imin,imax
        read(11,*) nzone
        if ( nzone.gt.maxnzone )
     &      stop 'nzone is too large. (pinput)'
        do 140 i=1,nzone
          read(11,*) vrmin(i),vrmax(i),
     &    rho(1,i),rho(2,i),rho(3,i),rho(4,i),
     &    vpv(1,i), vpv(2,i), vpv(3,i), vpv(4,i),
     &    vph(1,i), vph(2,i), vph(3,i), vph(4,i),
     &    vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i),
     &    vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i),
     &    eta(1,i), eta(2,i), eta(3,i), eta(4,i),
     &    qmu(i),qkappa(i)
  140   continue
        read(11,*) r0,eqlat,eqlon
        eqlattmp = eqlat
        call translat(eqlattmp,eqlattmp)
        read(11,*) mt(1,1),mt(1,2),mt(1,3),mt(2,2),mt(2,3),mt(3,3)
        read(11,*) nr
        if ( nr.gt.maxnr )
     &  stop 'nr is too large. (pinput)'
        do 150 i=1,nr
          read(11,*) lat(i),lon(i)
          stlat = lat(i)
          stlon = lon(i)
          call translat(stlat,stlat)
          call calthetaphi(eqlattmp,eqlon,stlat,stlon,theta(i),phi(i))
  150   continue
        do 160 i=1,nr
          read(11,110) output(i)
  160   continue
c temporary file close
        close(11)
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        real*8 pi
        parameter ( pi = 3.1415926535897932d0 )
c
        real*8 ievla,ievlo,istla,istlo
        real*8 evla,evlo,stla,stlo
        real*8 theta,phi
        real*8 gcarc,az
        real*8 tc,ts
c
c transformation to spherical coordinates
c
        evla = 90.d0 - ievla
        stla = 90.d0 - istla
c
        evla = evla / 1.8d2 * pi
        evlo = ievlo / 1.8d2 * pi
        stla = stla / 1.8d2 * pi
        stlo = istlo / 1.8d2 * pi
c
        gcarc = dacos( dcos(evla) * dcos(stla)
     &       + dsin(evla) * dsin(stla) * dcos(evlo - stlo) )
c
        tc = ( dcos(stla) * dsin(evla)
     &       - dsin(stla) * dcos(evla) * dcos(stlo - evlo) )
     &       / dsin(gcarc)
        ts = dsin(stla) * dsin(stlo - evlo) / dsin(gcarc)
c
        az = dacos(tc)
        if( ts .lt. 0.d0 ) az = -1.d0 * az
c
        az = az * 1.8d2 / pi

        gcarc = gcarc * 1.8d2 / pi
c
        theta = gcarc
        phi   = 180.d0 - az
        return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine translat(geodetic,geocentric)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        real*8 flattening,pi
c      parameter ( flattening = 1.d0 / 297.d0 )
        parameter ( flattening = 1.d0 / 298.25d0)
        parameter ( pi = 3.1415926535897932d0 )
        real*8 geocentric, geodetic

        real*8 tmp
        integer flag
c      read(5,*) geodetic
        flag = 0
        if(geodetic .gt. 90.d0) then
          geodetic = 1.8d2 - geodetic
          flag = 1
        endif
c
        geodetic = geodetic / 1.8d2 * pi
        geocentric = datan( (1.d0 - flattening) * (1.d0 - flattening)
     &       * dtan(geodetic) )
        geocentric = geocentric * 1.8d2 / pi
c      if(geocentric .lt. 0.d0 ) geocentric = 1.8d2 + geocentric
        if(flag .eq. 1) then
          geocentric = 1.8d2 - geocentric
        endif
c      write(6,*) 'geocentric latitude', geocentric
        return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calnl( nzone,vs,iphase,nsl,nll )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c counting of nsl and nll.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer nzone,iphase(*),nsl,nll
        real*8 vs(4,*)
        integer i
c
        nsl = 0
        nll = 0
        do 100 i=1,nzone
          if ( ( vs(1,i).eq.0.d0 ).and.( vs(2,i).eq.0.d0 ).and.
     &    ( vs(3,i).eq.0.d0 ).and.( vs(4,i).eq.0.d0 ) ) then
            nll = nll + 1
            iphase(i) = 2
          else
            nsl = nsl + 1
            iphase(i) = 1
          endif
  100   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calgrid( nzone,vrmin,vrmax,vp,vs,rmin,rmax,
     &imax,lmin,tlen,vmin,gridpar,dzpar )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real*8 pi
        parameter ( pi=3.1415926535897932d0 )
c
        integer nzone,imax,lmin
        real*8 vrmin(*),vrmax(*),vp(4,*),vs(4,*)
        real*8 rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
        integer izone,i,j
        real*8 coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp
c
        do 130 izone=1,nzone
c computing the S-velocity at each zone
          if ( vs(1,izone).eq.0.d0 ) then
            do 100 i=1,4
              v(i) = vp(i,izone)
  100       continue
          else
            do 110 i=1,4
              v(i) = vs(i,izone)
  110       continue
          endif
          vs1 = 0.d0
          vs2 = 0.d0
          do 120 j=1,4
            if ( j.eq.1 ) then
              coef1 = 1.d0
            else
              coef1 = coef1 * ( vrmin(izone) / rmax )
            endif
            if ( j.eq.1 ) then
              coef2 = 1.d0
            else
              coef2 = coef2 * ( vrmax(izone) / rmax )
            endif
            vs1 = vs1 + v(j) * coef1
            vs2 = vs2 + v(j) * coef2
  120     continue
c computing rh
          rh = vrmax(izone) - vrmin(izone)
c computing omega,amax
          omega = 2.d0 * pi * dble(imax) / tlen
          if ( vs1.ge.vs2 ) then
            vmin(izone) = vs2
          else
            vmin(izone) = vs1
          endif
          amax = vrmax(izone)
          gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) )
     &    - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) )
     &    / ( amax * amax )
          if ( gtmp.gt.0.d0 ) then
            dzpar(izone)   = dsqrt( 1.d0/gtmp )
            gridpar(izone) = rh / dzpar(izone)
          else
            dzpar(izone)   = 0.d0
            gridpar(izone) = 0.d0
          endif
  130   continue
c rearangement of gridpar
        gtmp = 0.d0
        do 140 izone=1,nzone
          gtmp = gtmp + gridpar(izone)
  140   continue
        do 150 izone=1,nzone
          if ( gridpar(izone).gt.0.d0 ) then
            gridpar(izone) = gridpar(izone) / gtmp
          else
            rh = vrmax(izone) - vrmin(izone)
            gridpar(izone) = rh / ( rmax - rmin ) * 0.1d0
          endif
  150   continue
c re-rearangement of gridpar
        gtmp = 0.d0
        do 160 izone=1,nzone
          gtmp = gtmp + gridpar(izone)
  160   continue
        do 170 izone=1,nzone
          gridpar(izone) = gridpar(izone) / gtmp
  170   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calra( maxnlay,maxnslay,maxnllay,maxnzone,
     &nlayer,inlayer,jnlayer,jnslay,jnllay,
     &gridpar,dzpar,nzone,vrmin,vrmax,iphase,
     &rmin,rmax,r0,nslay,nllay,nnl,ra,re )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the number and the location of grid points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        real*8 pi
        parameter ( pi=3.1415926535897932d0 )
c
        integer maxnlay,maxnslay,maxnllay,maxnzone
        integer nlayer,inlayer,jnlayer,jnslay,jnllay
        integer nzone,iphase(*),nslay,nllay,nnl(maxnzone)
        real*8 gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,r0
        real*8 ra(maxnlay+maxnzone+1)
        integer izone,itmp,i,ntmp
        real*8 rh,re
c
c Initializing the data
        nslay = 0
        nllay = 0
        inlayer = 0
        do 100 i=1,maxnlay+maxnzone+1
          ra(i) = 0.d0
  100   continue
        do 110 izone=1,nzone
          nnl(izone) = 0
  110   continue
        jnlayer = 0
        jnslay = 0
        jnllay = 0
c
c	tnlayer = nlayer / (2**(idr-1))
c computing the number and the location of the grid points
        ra(1) = rmin
        itmp = 1
        do 140 izone=1,nzone
          rh = vrmax(izone) - vrmin(izone)
          if(dzpar(izone).eq.0.d0) then
            ntmp = 1
          else
            ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone)
     &                   / 2.d0 / pi  / 7.d-1 + 1 )
          endif
c                             ! ntmp (see Geller & Takeuchi 1995 6.2)
c	  nnl(izone) = dint( dble(nlayer) * gridpar(izone) ) + 1
          nnl(izone) = ntmp
          if ( nnl(izone).lt.5 ) nnl(izone)=5
          if ( iphase(izone).eq.1 )
     &    nslay = nslay + nnl(izone)
          if ( nslay.gt.maxnslay )
     &         stop  'nslay is too large. (calra)'
          if ( iphase(izone).eq.2 )
     &    nllay = nllay + nnl(izone)
          if ( nllay.gt.maxnllay )
     &         stop  'nllay is too large. (calra)'
          do 130 i=1,nnl(izone)
            itmp = itmp + 1
            if ( itmp.gt.maxnlay )
     &          stop  'nlay is too large. (calra)'
            ra(itmp) = vrmin(izone)
     &      + rh * dble(i) / dble( nnl(izone) )
  130     continue
  140   continue
c
c recouting the total number of grid points
        inlayer = 0
        do 150 izone=1,nzone
          inlayer = inlayer + nnl(izone)
  150   continue
        jnlayer = jnlayer + inlayer
        jnslay  = jnslay  + nslay
        jnllay  = jnllay  + nllay
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calsp( maxnzone,ndc,nsl,nll,
     &iphase,nlayer,nslay,nllay,
     &isp,jsp,ksp,issp,ilsp,lsp,jssp,
     &isdr,jsdr,ildr,jdr,kdr )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the stack points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer maxnzone
        integer ndc,nsl,nll,iphase(*),nlayer(maxnzone)
        integer nslay,nllay
        integer isp(maxnzone),jsp(maxnzone),ksp(maxnzone)
        integer issp(maxnzone),ilsp(maxnzone)
        integer lsp(maxnzone),jssp(maxnzone)
        integer isdr,jsdr,ildr,jdr,kdr
        integer i,isl,ill
c
c Initialization of the data
        do 100 i=1,maxnzone
          isp(i)  = 0
          jsp(i)  = 0
          ksp(i)  = 0
          issp(i) = 0
          ilsp(i) = 0
          lsp(i)  = 0
          jssp(i) = 0
  100   continue
        isdr = 0
        jsdr = 0
        ildr = 0
        jdr = 0
        kdr = 0
c computation of isp,jsp,ksp,issp,ilsp,lsp
        isp(1)  = 1
        jsp(1)  = 1
        ksp(1)  = 1
        issp(1) = 1
        ilsp(1) = 1
        lsp(1)  = 1
        jssp(1) = 1
        isl = 0
        ill = 0
        do 120 i=1,ndc
          isp(i+1) = isp(i) + nlayer(i)
          if ( iphase(i).eq.1 ) then
            jsp(i+1) = jsp(i) + 16 * nlayer(i)
            ksp(i+1) = ksp(i) + 2 * ( nlayer(i) + 1 )
            lsp(i+1) = lsp(i) + 4 * nlayer(i)
            isl = isl + 1
            if ( isl.ne.nsl ) then
              issp(isl+1) = issp(isl) + 4 * nlayer(i)
              jssp(isl+1) = jssp(isl) + nlayer(i) + 1
            endif
          else
            jsp(i+1) = jsp(i) + 4 * nlayer(i)
            ksp(i+1) = ksp(i) + ( nlayer(i) + 1 )
            lsp(i+1) = lsp(i) + 2 * nlayer(i)
            ill = ill + 1
            if ( ill.ne.nll )
     &      ilsp(ill+1) = ilsp(ill) + 4 * nlayer(i)
          endif
  120   continue
        isdr = 0
        jsdr = 0
        ildr = 0
        jdr = 0
        isdr = isdr + issp(nsl)-1 + 4 * nlayer(ndc+1)
        jsdr = jsdr + jssp(nsl)-1 + nlayer(ndc+1) + 1
        ildr = ildr + 4 * nllay
        jdr =  jdr  + jsp(ndc+1)-1 + 16 * nlayer(ndc+1)
        kdr =  kdr + ksp(ndc+1)-1 + 2 * ( nlayer(ndc+1)+1 )
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calspo( maxnlay,maxnzone,
     &ndc,rdc,iphase,inlayer,
     &r0,rmin,rmax,ra,isp,spo,spn )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the source location.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer maxnlay,maxnzone,ndc,iphase(*)
        integer inlayer,isp(maxnzone),spn
        real*8 rdc(*),r0,rmin,rmax,ra(maxnlay+maxnzone+1),spo
        integer itmp,idr
c
c checking the parameter
        if ( (r0.lt.rmin).or.(r0.gt.rmax) )
     &  stop 'The source location is improper.(calspo)'
        spo = 0
c computing 'spo'
        if ( r0.eq.rmax ) then
          spo = dble( inlayer ) - 0.01d0
          r0 = ra(inlayer)
     &    + (spo-dble(inlayer-1)) * ( ra(inlayer+1) -ra(inlayer) )
        else
          itmp = 2
  110     continue
          if ( r0.lt.ra(itmp) ) then
            continue
          else
            itmp = itmp + 1
            goto 110
          endif
          spo = dble(itmp-2)
     &    + ( r0-ra(itmp-1) )   / ( ra(itmp)-ra(itmp-1) )
c temporal handling
          if ( (spo-dble(itmp-2)).lt.0.01d0 ) then
            spo = dble(itmp-2) + 0.01d0
            r0 = ra(itmp-1)
     &      + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
          endif
          if ( (spo-dble(itmp-2)).gt.0.99d0 ) then
            spo = dble(itmp-2) + 0.99d0
            r0 = ra(itmp-1)
     &      + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
          endif
        endif
c computing 'spn'
        spn = 0
        itmp = 1
  130   continue
        if ( iphase(itmp).eq.1 ) then
          spn = spn + 1
          if ( r0.le.rdc(itmp) ) then
            continue
          else
            itmp = itmp + 1
            goto 130
          endif
        else
          spn = spn + 1
          if ( r0.le.rdc(itmp) )
     &    stop 'The source is in the liquid layer.(calspo)'
          itmp = itmp + 1
          goto 130
        endif
c changing 'spo'
        spo = spo - dble( isp(spn) - 1 )
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine xmatinit( maxnslay,maxnllay,maxnzone,
     &a,t,h1x,h1y,h1z,h2L,h2N,
     &                       h3ax,h3ay,h3az,h4aL,h4aN,
     &                       h5ax,h5ay,h5az,h6aL,h6aN,
     &h3x,h3y,h3z,h4L,h4N,
     &                       h5x,h5y,h5z,h6L,h6N,
     &h7x,h7y,h7z,h8L,h8N,
     &h3mx,h3my,h3mz,h5mx,h5my,h5mz,
     &                       h4m1L,h4m1N,h4m2L,h4m2N,
     &                       h6m1L,h6m1N,h6m2L,h6m2N,
     &p1,p2,p3 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer maxnslay,maxnllay,maxnzone
        complex*16 a( 4, 2*(maxnslay+1) + (maxnllay+1) )
        real*8  t(8*maxnslay)
        real*8 h1x( 8*maxnslay ),h1y( 8*maxnslay ),h1z( 8*maxnslay )
        real*8 h2L( 8*maxnslay ),h2N( 8*maxnslay )
        real*8 h3ax( 8*maxnslay )
        real*8 h3ay( 8*maxnslay ),h3az( 8*maxnslay )
        real*8 h4aL( 8*maxnslay ),h4aN( 8*maxnslay )
        real*8 h5ax( 8*maxnslay ),h5ay( 8*maxnslay ),h5az( 8*maxnslay )
        real*8 h6aL( 8*maxnslay ),h6aN( 8*maxnslay )
        real*8 h3x( 8*maxnslay ),h3y( 8*maxnslay ),h3z( 8*maxnslay )
        real*8 h4L( 8*maxnslay ),h4N( 8*maxnslay )
        real*8 h5x( 8*maxnslay ),h5y( 8*maxnslay ),h5z( 8*maxnslay )
        real*8 h6L( 8*maxnslay ),h6N( 8*maxnslay )
        real*8 h7x( 8*maxnslay ),h7y( 8*maxnslay ),h7z( 8*maxnslay )
        real*8 h8L( 8*maxnslay ),h8N( 8*maxnslay )
        real*8  h3mx( -2:1,2*(maxnslay+maxnzone) )
        real*8  h3my( -2:1,2*(maxnslay+maxnzone) )
        real*8  h3mz( -2:1,2*(maxnslay+maxnzone) )
        real*8  h5mx( -1:2,2*(maxnslay+maxnzone) )
        real*8  h5my( -1:2,2*(maxnslay+maxnzone) )
        real*8  h5mz( -1:2,2*(maxnslay+maxnzone) )
        real*8 h4m1L( -1:2,2*(maxnslay+maxnzone) )
        real*8 h4m1N( -1:2,2*(maxnslay+maxnzone) )
        real*8 h4m2L( -2:1,2*(maxnslay+maxnzone) )
        real*8 h4m2N( -2:1,2*(maxnslay+maxnzone) )
        real*8 h6m1L( -1:2,2*(maxnslay+maxnzone) )
        real*8 h6m1N( -1:2,2*(maxnslay+maxnzone) )
        real*8 h6m2L( -2:1,2*(maxnslay+maxnzone) )
        real*8 h6m2N( -2:1,2*(maxnslay+maxnzone) )
        real*8 p1(8*maxnllay),p2(8*maxnllay),p3(8*maxnllay)
c
        call cmatinit( 4,2*(maxnslay+1)+(maxnllay+1),a )
        call vecinit( 8*maxnslay,t )
        call vecinit( 8*maxnslay,h1x )
        call vecinit( 8*maxnslay,h1y )
        call vecinit( 8*maxnslay,h1z )
        call vecinit( 8*maxnslay,h2L )
        call vecinit( 8*maxnslay,h2N )
        call vecinit( 8*maxnslay,h3ax )
        call vecinit( 8*maxnslay,h3ay )
        call vecinit( 8*maxnslay,h3az )
        call vecinit( 8*maxnslay,h4aL )
        call vecinit( 8*maxnslay,h4aN )
        call vecinit( 8*maxnslay,h5ax )
        call vecinit( 8*maxnslay,h5ay )
        call vecinit( 8*maxnslay,h5az )
        call vecinit( 8*maxnslay,h6aL )
        call vecinit( 8*maxnslay,h6aN )
        call vecinit( 8*maxnslay,h3x )
        call vecinit( 8*maxnslay,h3y )
        call vecinit( 8*maxnslay,h3z )
        call vecinit( 8*maxnslay,h4L )
        call vecinit( 8*maxnslay,h4N )
        call vecinit( 8*maxnslay,h5x )
        call vecinit( 8*maxnslay,h5y )
        call vecinit( 8*maxnslay,h5z )
        call vecinit( 8*maxnslay,h6L )
        call vecinit( 8*maxnslay,h6N )
        call vecinit( 8*maxnslay,h7x )
        call vecinit( 8*maxnslay,h7y )
        call vecinit( 8*maxnslay,h7z )
        call vecinit( 8*maxnslay,h8L )
        call vecinit( 8*maxnslay,h8N )
        call matinit( 4,2*(maxnslay+maxnzone),h3mx(-2,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h3my(-2,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h3mz(-2,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h4m2L(-2,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h4m2N(-2,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h6m2L(-2,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h6m2N(-2,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h5mx(-1,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h5my(-1,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h5mz(-1,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h4m1L(-1,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h4m1N(-1,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h6m1L(-1,1) )
        call matinit( 4,2*(maxnslay+maxnzone),h6m1N(-1,1) )
        call vecinit( 8*maxnllay,p1 )
        call vecinit( 8*maxnllay,p2 )
        call vecinit( 8*maxnllay,p3 )
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calstg( maxnlay,maxnzone,
     &nzone,iphase,rrho,
     &                     vpv,vph,vsv,vsh,eta,nnl,ra,rmax,
     &vnp,vra,rho,kappa,ecKx,ecKy,ecKz,
     &                     mu,ecL,ecN,
     &                     r0,spn,ecC0,ecF0,ecL0 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the structure grid points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer maxnlay,maxnzone,nzone,iphase(*),nnl(*),vnp,spn
        real*8 rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
        real*8 ra(*),rmax
        real*8 vra(*),rho(*),kappa(*),ecKx(*),ecKy(*),ecKz(*)
        real*8 mu(*),ecL(*),ecN(*)
        real*8 ecA,ecC,ecF
        real*8 r0,ecA0,ecC0,ecF0,ecL0
        real*8 trho,tvpv,tvph,tvsv,tvsh,teta,coef
        integer izone,i,j,itmp,jtmp
c
c initializing the data
        call vecinit( maxnlay+2*maxnzone+1,vra )
        call vecinit( maxnlay+2*maxnzone+1,rho )
        call vecinit( maxnlay+2*maxnzone+1,kappa )
        call vecinit( maxnlay+2*maxnzone+1,ecKx )
        call vecinit( maxnlay+2*maxnzone+1,ecKy )
        call vecinit( maxnlay+2*maxnzone+1,ecKz )
        call vecinit( maxnlay+2*maxnzone+1,mu )
        call vecinit( maxnlay+2*maxnzone+1,ecL )
        call vecinit( maxnlay+2*maxnzone+1,ecN )
c computing the structure grid points
        itmp = 0
        jtmp = 0
        do 130 izone=1,nzone
          do 120 i=1,nnl(izone)+1
            itmp = itmp + 1
            jtmp = jtmp + 1
            vra(itmp) = ra(jtmp)
c --- evaluating the density and elastic constants at this point
            trho = 0.d0
            tvpv = 0.d0
            tvph = 0.d0
            tvsv = 0.d0
            tvsh = 0.d0
            teta = 0.d0
            do 110 j=1,4
              if ( j.eq.1 ) then
                coef = 1.d0
              else
                coef = coef * ( vra(itmp) / rmax )
              endif
              trho  = trho  + rrho(j,izone)  * coef
              tvpv  = tvpv  + vpv(j,izone)   * coef
              tvph  = tvph  + vph(j,izone)   * coef
              tvsv  = tvsv  + vsv(j,izone)   * coef
              tvsh  = tvsh  + vsh(j,izone)   * coef
              teta  = teta  + eta(j,izone)   * coef
  110       continue
            rho(itmp) = trho
            ecL(itmp)  = rho(itmp) * tvsv * tvsv
            ecN(itmp)  = rho(itmp) * tvsh * tvsh
            ecA = trho * tvph * tvph
            ecC = trho * tvpv * tvpv
            ecF = teta * ( ecA - 2.d0 * ecL(itmp) )
            kappa(itmp) = ( 4.d0 * ecA + ecC
     &                 + 4.d0 * ecF - 4.d0 * ecN(itmp) ) / 9.d0
            ecKx(itmp) = ecA - 4.d0 / 3.d0 * ecN(itmp)
            ecKy(itmp) = ecF + 2.d0 / 3.d0 * ecN(itmp)
            ecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
  120     continue
          jtmp = jtmp - 1
  130   continue
        vnp = itmp
c
        trho = 0.d0
        tvpv = 0.d0
        tvph = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        teta = 0.d0
        do 140 j=1,4
          if ( j.eq.1 ) then
            coef = 1.d0
          else
            coef = coef * ( r0 / rmax )
          endif
          trho  = trho  + rrho(j,spn)  * coef
          tvpv  = tvpv  + vpv(j,spn)   * coef
          tvph  = tvph  + vph(j,spn)   * coef
          tvsv  = tvsv  + vsv(j,spn)   * coef
          tvsh  = tvsh  + vsh(j,spn)   * coef
          teta  = teta  + eta(j,spn)   * coef
  140   continue
        ecL0 = trho * tvsv * tvsv
        ecA0 = trho * tvph * tvph
        ecC0 = trho * tvpv * tvpv
        ecF0 = teta * ( ecA0 - 2.d0 * ecL0 )
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine caltstg( maxnlay,maxnzone,
     &nzone,rrho,vpv,vph,vsv,vsh,eta,
     &                      nnl,ra,rmax,
     &tvra,tkappa,tecKx,tecKy,tecKz,
     &                      tmu,tecL,tecN)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the structure grid points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        integer maxnlay,maxnzone,nzone,nnl(*)
        real*8 rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
        real*8 ra(*),rmax
        real*8 tvra(*),tkappa(*),tmu(*)
        real*8 tecKx(*),tecKy(*),tecKz(*),tecL(*),tecN(*)
        real*8 trho,tvpv,tvph,tvsv,tvsh,teta,coef
        real*8 ecA,ecC,ecF
        integer izone,i,j,itmp,jtmp
c
        call vecinit( maxnlay+2*maxnzone+1,tvra )
        call vecinit( maxnlay+2*maxnzone+1,tkappa )
        call vecinit( maxnlay+2*maxnzone+1,tecKx )
        call vecinit( maxnlay+2*maxnzone+1,tecKy )
        call vecinit( maxnlay+2*maxnzone+1,tecKz )
        call vecinit( maxnlay+2*maxnzone+1,tmu )
        call vecinit( maxnlay+2*maxnzone+1,tecL )
        call vecinit( maxnlay+2*maxnzone+1,tecN )
c computing the structure grid points
        itmp = 0
        jtmp = 0
        do 130 izone=1,nzone
          do 120 i=1,nnl(izone)+1
            itmp = itmp + 1
            jtmp = jtmp + 1
            tvra(itmp) = ra(jtmp)
c --- evaluating the density and elastic constants at this point
            trho = 0.d0
            tvpv = 0.d0
            tvph = 0.d0
            tvsv = 0.d0
            tvsh = 0.d0
            teta = 0.d0
            do 110 j=1,4
              if ( j.eq.1 ) then
                coef = 1.d0
              else
                coef = coef * ( tvra(itmp) / rmax )
              endif
              trho = trho + rrho(j,izone) * coef
              tvpv  = tvpv  + vpv(j,izone)   * coef
              tvph  = tvph  + vph(j,izone)   * coef
              tvsv  = tvsv  + vsv(j,izone)   * coef
              tvsh  = tvsh  + vsh(j,izone)   * coef
              teta  = teta  + eta(j,izone)   * coef
  110       continue
            tecL(itmp)  = trho * tvsv * tvsv
            tecN(itmp)  = trho * tvsh * tvsh
            ecA = trho * tvph * tvph
            ecC = trho * tvpv * tvpv
            ecF = teta * ( ecA - 2.d0 * tecL(itmp) )
            tkappa(itmp) = ( 4.d0 * ecA + ecC
     &               + 4.d0 * ecF - 4.d0 * tecN(itmp) )/ 9.d0
            tecKx(itmp) = ecA - 4.d0 / 3.d0 * tecN(itmp)
            tecKy(itmp) = ecF + 2.d0 / 3.d0 * tecN(itmp)
            tecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
  120     continue
          jtmp = jtmp - 1
  130   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calinv(vnp,rho,kappa,rhoinv,kappainv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the inverse of density and elastic constant.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer vnp,i
        real*8 rho(*),kappa(*),rhoinv(*),kappainv(*)
c
        do 100 i=1,vnp
          rhoinv(i)   = 1.d0 / rho(i)
          kappainv(i) = 1.d0 / kappa(i)
  100   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine submat( nlayer,ha,hb,h )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subtracting matrix `hb' from matrix `ha'.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer nlayer
        real*8 ha(*),hb(*),h(*)
        integer i
c
        do 100 i=1,4*nlayer
          h(i) = ha(i) - hb(i)
  100   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calspdr( maxnzone,nzone,iphase,nlayer,
     &jjdr,kkdr )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer maxnzone,nzone,iphase(*)
        integer nlayer(maxnzone),jjdr(*),kkdr(*)
        integer izone
c
        jjdr(1) = 1
        kkdr(1) = 1
        do 100 izone=1,nzone-1
          if ( iphase(izone).eq.1 ) then
            jjdr(izone+1) = jjdr(izone) + 16 * nlayer(izone)
            if ( iphase(izone+1).eq.1 ) then
              kkdr(izone+1)
     &             = kkdr(izone) + 2 * nlayer(izone)
            else
              kkdr(izone+1)
     &             = kkdr(izone) + 2 * ( nlayer(izone)+1 )
            endif
          else
            jjdr(izone+1) = jjdr(izone) + 4 * nlayer(izone)
            if ( iphase(izone+1).eq.1 ) then
              kkdr(izone+1)
     &             = kkdr(izone) + ( nlayer(izone)+1 )
            else
              kkdr(izone+1)
     &             = kkdr(izone) + nlayer(izone)
            endif
          endif
  100   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,
     &rmax,sufzone )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real*8 pi
        parameter ( pi=3.1415926535897932d0 )
c
        integer l,nzone,sufzone
        real*8 omega,vrmin(*),vrmax(*),vmin(*),dzpar(*),rmax
        integer izone
        real*8 gtmp,tdzpar
c
        sufzone = 0
        do 110 izone=1,nzone
          gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) )
     &    - ( (dble(l)+0.5d0) * (dble(l)+0.5d0) )
     &    / ( vrmax(izone) * vrmax(izone) )
          if ( gtmp.gt.0.d0 ) then
            tdzpar = dsqrt( 1.d0/gtmp )
          else
            if ( vrmax(izone).gt.rmax*(1-2.d0*pi/(dble(l)+0.50)) )
     &      then
              tdzpar = 0.d0
            else
              sufzone = izone
              tdzpar = 0.d0
            endif
          endif
  110   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calu0( c0,bvec,u )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        complex*16 c0,bvec,u
c
        u = u + c0 * bvec
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calu( c0,lsq,bvec,u )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real*8 lsq
        complex*16 c0(2),bvec(3),u(3)
c
        u(1) = u(1) + c0(1) * bvec(1)
        u(2) = u(2) + c0(2) * bvec(2) / dcmplx(lsq)
        u(3) = u(3) + c0(2) * bvec(3) / dcmplx(lsq)
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calamp( g,l,lsuf,maxamp,ismall,ratl )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        integer ismall,l,lsuf
        real*8 maxamp,ratl
        complex*16 g(2)
        real*8 amp,ampratio
c
        ampratio = 0.d0
        amp = dsqrt( cdabs( g(1) )**2 + cdabs( g(2) )**2 )
        if ( amp.gt.maxamp ) maxamp = amp
        if ( (amp.ne.0.d0).and.(maxamp.ne.0.d0) )
     &  ampratio = amp / maxamp
        if ( ( ampratio.lt.ratl ).and.( l.ge.lsuf ) ) then
          ismall = ismall + 1
        else
          ismall = 0
        endif
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matinit( n1,n2,a )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer n1,n2,i,j
        real*8 a(n1,*)
c
        do 110 j=1,n2
          do 100 i=1,n1
            a(i,j) = 0.d0
  100     continue
  110   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmatinit( n1,n2,a )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer n1,n2,i,j
        complex*16 a(n1,*)
c
        do 110 j=1,n2
          do 100 i=1,n1
            a(i,j) = dcmplx( 0.d0 )
  100     continue
  110   continue

        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vecinit( nn,b )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Filling zero to the vector 'g'.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer nn,i
        real*8 b(*)
c
        do 100 i=1,nn
          b(i) = 0.d0
  100   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cvecinit( nn,b )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Filling zero to the vector 'g'.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer nn,i
        complex*16 b(*)
c
        do 100 i=1,nn
          b(i) = dcmplx( 0.d0 )
  100   continue
c
        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calcutd(nzone,nnl,tmpc,rat,nn,iphase,spo,spn,
     &                           ra,kkdr,kc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer nzone,nn,spn,kkdr(*),kc,iphase(*),nnl(*)
        complex*16 tmpc(*)
        real*8 rat,spo,ra(*)
        integer nc
c
        real*8 cU(nn),cV(nn),rc
        real*8 maxamp,amp(nn)
        integer iz,jz,jj,i,ml(nzone),tzone
c
        do 100 jj=1,nn
          cU(jj) = 0.d0
          cV(jj) = 0.d0
  100   continue
c
        iz = 2
        jz = 1
        do 110 jj=1,nn
          if(jj.eq.kkdr(iz)) then
            if(iphase(iz).ne.iphase(iz-1)) jz = jz - 1
            iz = iz + 1
          endif
          if(iphase(iz-1).eq.1) then
            if(mod((jj-kkdr(iz-1)),2).eq.1) then  ! U
              cU(jz) = cdabs(tmpc(jj))
              jz = jz + 1
            else  ! V
            endif
          else  ! U in fluid
            cU(jz) = cdabs(tmpc(jj))
            jz = jz + 1
          endif
  110   continue
c
        maxamp = -1.d0
        do 120 i=1,jz-1
          amp(i) = cU(i)
          if(maxamp.lt.amp(i)) maxamp = amp(i)
  120   continue
c
        maxamp = maxamp * rat  ! threshold value
c
        nc = 1
        do 130 i=1,jz-1
          if(amp(i).gt.maxamp) then
            nc = i
            goto 140
          endif
  130   continue
  140   continue
c
        i = 1
        do 150 jj=1,nzone
          i = i + nnl(jj)
          ml(jj) = i
  150   continue
c
        do 160 jj=nzone,1,-1
          if(ml(jj).gt.nc) tzone = jj
  160   continue
c
        rc = ra(nc)
c
        do 590 i=1,jz-1
          if( (ra(i).le.rc).and.(rc.lt.ra(i+1)) ) then
            nc = i
            if(tzone.eq.1) then ! case(tzone is innermost zone)
              if(iphase(tzone).eq.1) kc = 1 + 2 * nc
              if(iphase(tzone).eq.2) kc = 1 + nc
            else
              if(iphase(tzone).eq.1) then
                kc = kkdr(tzone) + 2 * (nc - ml(tzone-1))
              endif
              if(iphase(tzone).eq.2) then
                kc = kkdr(tzone) + nc - ml(tzone-1)
              endif
            endif
          endif
  590   continue

        return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit none
        integer nzone,lsuf
        real*8 omega,vrmax(*),vsv(4,*)
c
        real*8 tvs,coef
        integer i
c
        tvs = 0.d0
        do 100 i=1,4
          if(i.eq.1) then
            coef = 1.d0
          else
            coef = coef
          endif
          tvs = tvs + ( vsv(i,nzone) ) * coef
  100   continue
c
        lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
        return
      end
c
