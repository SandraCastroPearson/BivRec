c------------------------------------------------------------
c     nrowtmp     I  
c     ncoltmp     I  
c     tmp         I  
c------------------------------------------------------------
c     setup all variables

      subroutine changmdat(tugap1, tugap2, tdat, tmpin,   
     +ugapcols, tmpstart, tmpend, tugapstart, tugapend, 
     +idcount, nrowtdat, nrowugap, maxidcount, 
     +n, tmpcols, tdatcols)

      implicit none

c     given variables
      integer n, tmpcols
      integer tmpstart, tmpend, tugapstart, tugapend
      integer tdatcols, ugapcols
      integer idcount(n), nrowtdat, nrowugap, maxidcount
      double precision tmpin(maxidcount, tmpcols)
      double precision tdat(nrowtdat, tdatcols)
      double precision tugap1(nrowugap, ugapcols)
      double precision tugap2(nrowugap, ugapcols)

c     internal variables
      integer tdx, tdz, i, j, m, mstar, tc
      integer cti, a , b 
      double precision txij, tzij, tsum, comp
      double precision tmp(maxidcount, tmpcols)

c     Top loop
      do 10 i=1,n 

c        setup indexes for tmp and tugap

         if (i==1) then
            tmpstart = 1
            tmpend = idcount(1)
            tugapstart = 1
            tugapend = maxidcount
         else 
            tmpstart = tmpend + 1
            tmpend = tmpend + idcount(i)
            tugapstart = tugapend + 1
            tugapend = tugapend + maxidcount
         endif

c------------------------------------------------------------
c        setup tmp
c        column names for tmp: "id"   "epi"  "txij" "tzij" "tci"  "uxij" "uzij" "udx"  "udz" "covariates"

	 a = tmpstart
	 b = tmpend
         cti = idcount(i)

 	 tmp = tmpin
         tmp(1:cti, 1:5)=tdat(a:b,1:5)
c        uxij
         tmp(1,6)=minval((/ tmp(1,3), tmp(1, 5) /))
c        uzij
         tmp(1,7)=minval((/ tmp(1,4), tmp(1, 5) /))
c        udx
	 if (tmp(1,6).lt.tmp(1, 5)) then 
         	tmp(1,8) = dble(1)
	 else 
         	tmp(1,8) = dble(0)
	 end if
c        udz
	 if (tmp(1,7).lt.tmp(1, 5)) then 
         	tmp(1,9) = dble(1)
	 else 
         	tmp(1,9) = dble(0)	 
	 end if 
c        covariates
         tmp(1:cti, 10:tmpcols) = tdat(a:b, 6:tdatcols)

c------------------------------------------------------------
c        second loop: main for calculations
c        column names for tmp: "id"   "epi"  "txij" "tzij" "tci"  "uxij" "uzij" "udx"  "udz" 

         mstar=1
         if (cti.gt.1) then
            tdx=int(tmp(1,8))
            tdz=int(tmp(1,9))
            j=2

            do 20 while (j.le.cti.and.tdx.eq.1.and.tdz.eq.1)

c	        if (tdx.eq.1.and.tdz.eq.1) then

                   tsum=sum(tmp(1:cti,7))
                   comp = tmp(1,5)-tsum
                   txij=minval((/ tmp(j,3), comp /))
                   tzij=minval((/ tmp(j,4), comp /))

                   if ((txij+tsum).lt.tmp(1,5)) then
	              tdx = 1
                   else 
	              tdx = 0 
                   end if
                   if ((tzij+tsum).lt.tmp(1,5)) then
	              tdz = 1
                   else 
	              tdz = 0 
                   end if
 
                   if (tdx.eq.1.and.tdz.eq.1) then
                     tmp(j,6)=txij
                     tmp(j,7)=tzij
                     tmp(j,8)=dble(tdx)
                     tmp(j,9)=dble(tdz)
                     mstar = mstar + 1
                   end if

	        j =j+1
c	        end if   

20          continue

         end if

c------------------------------------------------------------
c        fill in tugap

         tc = tmpcols
         tugap1(tugapstart:tugapend,1) = tmp(:,6)
         tugap1(tugapstart:tugapend,2) = tmp(:,8)
         tugap1(tugapstart:tugapend,4:ugapcols)=tmp(:,10:tc)

         tugap2(tugapstart:tugapend,1)= tmp(:,7)
         tugap2(tugapstart:tugapend,2)= tmp(:,9)
         tugap2(tugapstart:tugapend,4:ugapcols)=tmp(:,10:tc)

         do 40 m = tugapstart, tugapend
            tugap2(m, 3) = mstar
            tugap1(m, 3) = mstar
40       continue

10    continue

      return
      end