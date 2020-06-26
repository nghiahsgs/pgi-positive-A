      subroutine fdm_correlation_time_ave(escale,eig0,zmod,itype,nshell,
     +                           it,myid,nz,new_comm,T0)
      implicit none
      include 'mpif.h'
      include 'openacc_lib.h'

      integer ierr, new_comm
      integer m0
      common/m0_start_disc/m0

      integer nshell_max
      parameter(nshell_max=100)
      double precision wz(nshell_max), pz(nshell_max)
      common/wz_and_pz/wz,pz

      integer mdim, ishell, ielec, ispin, itype, ne, ie,
     +     minae, maxae, minas, maxas,  nshell, mx, it, myid, nz,
     +     ivx, jvp, jvm, ivr, jvpr, jvmr, itime, ntime

      parameter(ne=112)
      parameter(ntime=1)
      double precision escale(nshell), eig0(nshell),zmod(nshell),
     +                 fdm_spectra1(ntime,2*ne), ene(2*ne),omega0, 
     +                 fdm_spectra2(ntime,2*ne), 
     +                 coset(2*ne,ntime), sinet(2*ne,ntime),
     +                 alpha, srpi, T0, xtime(ntime)

      parameter(mdim = 1280)
      double precision dm(mdim,mdim),dmp(mdim,mdim),dmm(mdim,mdim),
     +                 eig(mdim),eigp(mdim),eigm(mdim),
     +                 vp(mdim,mdim),vm(mdim,mdim),
     +                 gp(mdim,mdim),gm(mdim,mdim)

      integer ns, nsr
      integer t_x, t_m,t_p, n_m

      double precision collector1(ntime,2*ne), collector2(ntime,2*ne)

      character(len=8), parameter :: fmt1 = '(I3.3)'
      character(len=60) :: fmt2
      character*100 datfl,datfl2
      character*3 x1, x2

      double precision  temp, mu
      common/temp_only/temp,mu
        
      integer c_ivx,c_jvp,c_jvm

      integer c_ivx_less_10
      integer c_ivx_less_100
      integer c_ivx_less_1000

      integer c_jvp_less_10
      integer c_jvp_less_100
      integer c_jvp_less_1000
      
      integer c_jvm_less_10
      integer c_jvm_less_100
      integer c_jvm_less_1000
      
      integer c_ivx_gt_1000
      integer c_jvp_gt_1000
      integer c_jvm_gt_1000

      call mpi_allreduce(m0,mx,1,mpi_integer,mpi_max,new_comm,ierr)
      if (mx.ne.m0)stop !checking m0 at different iz

      write(x2,fmt1)it
      fmt2='(3E25.15,E25.15E3)'

      fdm_spectra1(:,:) = 0.d0
      fdm_spectra2(:,:) = 0.d0
      collector1(:,:) = 0.d0
      collector2(:,:) = 0.d0

      do itime = 1, ntime, +1
         xtime(itime) = 0.d0!1.d-4/T0 * 10.d0**((itime-1.d0)/1.d0)
      enddo

      do ie = 1, ne 
         ene(ie) = -1.d0*10.d0**(-(ie-1)/12.d0)!*1.0d4*TK
         ene(2*ne+1-ie) = -ene(ie)
         do itime = 1, ntime
            coset(ie,itime) = cos(2.d0*ene(ie)*xtime(itime))
            coset(2*ne+1-ie,itime) = coset(ie,itime)
            sinet(ie,itime) = sin(2.d0*ene(ie)*xtime(itime))
            sinet(2*ne+1-ie,itime) = -sinet(ie,itime)
         enddo
      enddo

      omega0 = min(0.5d0*T0,0.5d0*temp)
c      alpha = 0.03125d0
      alpha = 0.3d0
      srpi = 8.d0*atan(1.d0)

      !do ishell = m0, nshell, +1
      do ishell = 24,25, +1

      t_x = 0
      t_m = 0
      t_p = 0
      n_m = 0

      c_ivx_less_10=0
      c_ivx_less_100=0
      c_ivx_less_1000=0
      c_ivx_gt_1000=0

      c_jvp_less_10=0
      c_jvp_less_100=0
      c_jvp_less_1000=0
      c_jvp_gt_1000=0

      c_jvm_less_10=0
      c_jvm_less_100=0
      c_jvm_less_1000=0
      c_jvm_gt_1000=0

        c_ivx=0
        c_jvp=0
        c_jvm=0

         minae = 0
         maxae = 2*ishell

         do ielec = minae, maxae
            minas = max(mod(ielec,2),0)
            maxas = min(2 * ishell - ielec, ielec)

            do ispin = minas, maxas, +2
               ivx = ns(ishell,ielec,ispin)
               jvp = ns(ishell,ielec-1,ispin+1)
               jvm = ns(ishell,ielec-1,ispin-1)

               if(((ivx*jvp).gt.0).or.((ivx*jvm).gt.0))then
                      c_ivx=c_ivx+1

                      if ((ivx.gt.0).and.(ivx.lt.10))  then
                               c_ivx_less_10=c_ivx_less_10+ivx
                      endif
                      if ((ivx.ge.10).and.(ivx.lt.100))  then
                               c_ivx_less_100=c_ivx_less_100+ivx
                      endif
                      if ((ivx.ge.100).and.(ivx.lt.1000))  then
                               c_ivx_less_1000=c_ivx_less_1000+ivx
                      endif

                      if ((ivx.ge.100).and.(ivx.lt.1000))  then  
                               c_ivx_gt_1000=c_ivx_gt_1000+1
                      endif

                if (jvp.gt.0)  then
                      c_jvp=c_jvp+1

                      if (jvp.lt.10)  then
                               c_jvp_less_10=c_jvp_less_10+1
                      endif
                      if (jvp.lt.100)  then
                               c_jvp_less_100=c_jvp_less_100+1
                      endif
                      if (jvp.lt.1000)  then
                               c_jvp_less_1000=c_jvp_less_1000+1
                      endif
                      if (jvp.gt.1000)  then
                               c_jvp_gt_1000=c_jvp_gt_1000+1
                      endif
                endif
                if (jvm.gt.0)  then
                      c_jvm=c_jvm+1

                      if (jvm.lt.10)  then
                               c_jvm_less_10=c_jvm_less_10+1
                      endif
                      if (jvm.lt.100)  then
                               c_jvm_less_100=c_jvm_less_100+1
                      endif
                      if (jvm.lt.1000)  then
                               c_jvm_less_1000=c_jvm_less_1000+1
                      endif
                      if (jvm.gt.1000)  then
                               c_jvm_gt_1000=c_jvm_gt_1000+1
                      endif
                endif
               endif
               !  t_x=t_x+ivx
               !  t_m=t_m+jvp
               !  t_p=t_p+jvm

         !   if (ivx.gt.0) then
         !     if ((jvp.gt.0).or.(jvm.gt.0)) then
         !     ivr = nsr(ishell,ielec,ispin)
         !     jvpr = nsr(ishell,ielec-1,ispin+1)
         !     jvmr = nsr(ishell,ielec-1,ispin-1)
         !     call sub_fdm_correlation(dm,dmp,dmm,eig,eigp,eigm,vp,
         !       +vm,gp,gm,ivx,jvp,jvm,ivr,jvpr,jvmr,ishell,ielec,
         !        +ispin,escale(ishell),wz(ishell),eig0(ishell),
         !        +zmod(ishell),alpha,omega0,srpi,fdm_spectra1,myid,
         !        +fdm_spectra2,itype,ntime,ne,T0,coset,sinet)
         !     endif
         !   endif

            enddo
         enddo
        write(*,'(6I6)') ishell, c_ivx_less_10, c_ivx_less_100,
      + c_ivx_less_1000, c_ivx_gt_1000,c_ivx
        write(*,'(6I6)') ishell, c_jvp_less_10, c_jvp_less_100,
      +                    c_jvp_less_1000,  c_jvp_gt_1000,c_jvp
        write(*,'(6I6)') ishell, c_jvm_less_10, c_jvm_less_100,
        +                   c_jvm_less_1000,  c_jvm_gt_1000,c_jvm

       ! List to do:
       ! 1/build up mtrix
       ! 2/cacl comutator => in order to update fdm_spectra1

        
      enddo
      !   print *,"nghiahsgs"
      !   print *,t_x,t_m,t_p
      call MPI_Reduce(fdm_spectra1, collector1, 2*ntime*ne,
     +                  MPI_double_precision,MPI_SUM, 0, new_comm, ierr)
      call MPI_Reduce(fdm_spectra2, collector2, 2*ntime*ne,
     +                  MPI_double_precision,MPI_SUM, 0, new_comm, ierr)

      datfl2 = 'Spectra_temp'//trim(x2)//'.dat'
      if (mod(myid,nz).eq.0) then
          print *,'Greens function in ',datfl2
          open(1222, file = datfl2, form = 'formatted')

          do itime = 1, ntime

             do ie = 1, 2*ne 
                write(1222,fmt2)xtime(itime),ene(ie),
     +                   collector1(itime,ie)/(nz*1.d0)
             enddo

             write(1222,*)
          enddo

          close(1222)
      endif

c      stop

      return
      end

      subroutine sub_fdm_correlation(dm,dmp,dmm,eig,eigp,eigm,vp,
     +           vm,gp,gm,ivx,jvp,jvm,ivr,jvpr,jvmr,ishell,ielec,ispin,
     +        escale,wz,eig0,zmod,alpha,omega0,srpi,fdm_spectra1,myid,
     +       fdm_spectra2,itype0,ntime,ne,T0,coset,sinet)
      implicit none
      include 'openacc_lib.h'
      integer ivx,jvp,jvm,ivr,jvpr,jvmr,ishell,ielec,ispin, itype0,
     +        itype, itype2, ir, is, ik, ntime, itime,myid, ne
      double precision dm(ivx,ivx),dmp(jvp,jvp),dmm(jvm,jvm),
     +                 eig(ivx),eigp(jvp),eigm(jvm), T0,
     +                 fdm_spectra1(ntime,2*ne), escale,
     +                 fdm_spectra2(ntime,2*ne), 
     +                 coset(2*ne,ntime), sinet(2*ne,ntime),
     +                 vp(ivx,jvp), vm(ivx,jvm), cleb, omega0,
     +                 gp(jvm,ivx), gm(jvp,ivx), alpha, srpi,
     +                 wz, eig0, zmod,  fp1, fm1

      double precision  temp, mu
      common/temp_only/temp,mu

      integer mdim
      parameter(mdim = 1280)
      double precision dmc(ntime,mdim,mdim), dmmc(ntime,mdim,mdim), 
     +      dmpc(ntime,mdim,mdim)
      double precision dms(ntime,mdim,mdim), dmms(ntime,mdim,mdim),
     +      dmps(ntime,mdim,mdim)

      logical gotcha

      cleb = ispin*0.5d0+0.5d0
      fp1 = -0.5d0*sqrt((ispin+1.0d0)*(ispin+2.0d0))
      fm1 = +0.5d0*sqrt((ispin)*(ispin+1.0d0))

    
      itype = + 4
      call read_data(eig, ivx, 1, ishell, ielec, ispin, itype)
      if (jvp.gt.0)then
         call read_data(eigp, jvp, 1, ishell, ielec-1, ispin+1, itype)
         itype2 = +1
         call read_data(vp, ivx, jvp, ishell, ielec, ispin, itype2)
         itype2 = -2
         call read_data(gm, jvp, ivx, ishell, ielec-1, ispin+1, itype2)
      endif

      if (jvm.gt.0)then
         call read_data(eigm, jvm, 1, ishell, ielec-1, ispin-1, itype)
         itype2 = -1
         call read_data(vm, ivx, jvm, ishell, ielec, ispin, itype2)
         itype2 = +2
         call read_data(gp, jvm, ivx, ishell, ielec-1, ispin-1, itype2)
      endif

      call setz2d(dm,ivx,ivx)
      call setz2d(dmp,jvp,jvp)
      call setz2d(dmm,jvm,jvm)

      itype = +8
      call test_existence_rho0m(ishell,ielec,ispin, itype,gotcha)
      if (gotcha.eqv..true.) then
         call read_data(dm,ivx,ivx,ishell,ielec,ispin,itype)
      endif
      call test_existence_rho0m(ishell,ielec-1,ispin+1, itype,gotcha)
      if (gotcha.eqv..true.) then
      call read_data(dmp,jvp,jvp,ishell,ielec-1,ispin+1,itype)
      endif
      call test_existence_rho0m(ishell,ielec-1,ispin-1, itype,gotcha)
      if (gotcha.eqv..true.) then
      call read_data(dmm,jvm,jvm,ishell,ielec-1,ispin-1,itype)
      endif

      call commutator(dm,dmp,dmm,eig,eigp,eigm,vp,vm,gp,gm,fp1,fm1,
     +  ivx,jvp,jvm,escale,cleb,alpha,omega0,srpi,fdm_spectra1,
     +  fdm_spectra2,ntime,dmc,dmmc,dmpc,dms,dmms,dmps,
     +  ivr,jvpr,jvmr,ne,T0,coset,sinet,myid)
   
      return
      end


      subroutine commutator(dm,dmp,dmm,eig,eigp,eigm,vp,vm,gp,gm,fp1,
     +  fm1,ivx,jvp,jvm,escale,cleb,alpha,omega0,srpi,fdm_spectra1,
     +  fdm_spectra2,ntime,dmc,dmmc,dmpc,dms,dmms,dmps,
     +  ivr,jvpr,jvmr,ne,T0,coset,sinet,myid)
      implicit none
      include 'openacc_lib.h'

      integer ivx, jvp, jvm, ir, is, iq, ie, ntime, it, ivr, jvpr, jvmr,
     +        ne,myid
      double precision dm(ivx,ivx), dmp(jvp,jvp), dmm(jvm,jvm),
     +        eig(ivx), eigp(jvp), eigm(jvm), 
     +        vp(ivx,jvp), vm(ivx,jvm), gp(jvm,ivx), gm(jvp,ivx), 
     +        fdm_spectra1(ntime,2*ne), cleb,fp1,fm1,escale,ene(2*ne), 
     +        fdm_spectra2(ntime,2*ne),
     +                 coset(2*ne,ntime), sinet(2*ne,ntime),
     +        exc, weight, omega0, c1con, c0con, cef1, cef2, wtm,
     +        alpha, T0, eta, srpi, xtime(ntime), weight1, weight2,
     +     dmc(ntime,ivx,ivx), dmmc(ntime,ivx,jvm), dmpc(ntime,ivx,jvp),
     +     dms(ntime,ivx,ivx), dmms(ntime,ivx,jvm), dmps(ntime,ivx,jvp)
      character typ
      double precision wxt, etaxt, exp2met,zv2,ext, decay, cdecay
      double precision zw1, zw2, zw(2*ne)

      double precision  temp, mu
      common/temp_only/temp,mu

      typ = 'O'

      c0con = 0.d0
      c1con = 1.d0
      cef1 = cleb
c      cef1 = fm1
      cef2 = cleb
c      cef1 = fp1

      do ir = 1, ivx
      do is = 1, ivx
      dm(ir,is) = dm(ir,is)*cleb
      enddo
      enddo

!$acc data copy(fdm_spectra1(ntime,2*ne))
!$acc& copyin(vp(ivx,jvp), vm(ivx,jvm),dm(ivx,ivx))
!$acc& create(dmc(ntime,ivx,ivx),dmmc(ntime,ivx,jvm))
!$acc& create(dmms(ntime,ivx,jvm),dmpc(ntime,ivx,jvp))
!$acc& create(dmps(ntime,ivx,jvp),xtime(ntime),ene(2*ne))
!$acc& create(coset(2*ne,ntime),sinet(2*ne,ntime))
!$acc kernels

!$acc loop vector(32)
      do it = 1, ntime, +1
         xtime(it) = 1.d-4/T0 * 10.d0**((it-1.d0)/1.d0)
      enddo

!$acc loop 
      do ie = 1, ne 
         ene(ie) = -1.d0*10.d0**(-(ie-1)/12.d0)!*1.0d4*TK
         ene(2*ne+1-ie) = -ene(ie)
!$acc loop vector(32)
        do it = 1, ntime,+1
            coset(ie,it) = cos(2.d0*ene(ie)*xtime(it))
            coset(2*ne+1-ie,it) = coset(ie,it)
            sinet(ie,it) = sin(2.d0*ene(ie)*xtime(it))
            sinet(2*ne+1-ie,it) = -sinet(ie,it)
         enddo
      enddo

!$acc loop collapse(3) vector(32)
      do it = 1, ntime
      do ir = 1, ivx
      do is = 1, ivx
      ext = (eig(ir)-eig(is))*escale*xtime(it)
      dmc(it,ir,is) = cos(ext)
c      dms(it,ir,is) = sin(ext)
      enddo
      enddo
      enddo

!$acc loop collapse(3)  vector(32)
      do it = 1, ntime
      do ir = 1, ivx
      do is = 1, jvm
      ext = 2.d0*(eigm(is)-eig(ir))*escale*xtime(it)
      dmmc(it,ir,is)=cos(ext)
      dmms(it,ir,is)=sin(ext)
      enddo
      enddo
      enddo

!$acc loop collapse(3)  vector(32)
      do it = 1, ntime
      do ir = 1, ivx
      do is = 1, jvp
      ext = 2.d0*(eigp(is)-eig(ir))*escale*xtime(it)
      dmpc(it,ir,is)=cos(ext)
      dmps(it,ir,is)=sin(ext)
      enddo
      enddo

      enddo

!$acc loop collapse(2) 
!$acc& private(zw)
      do ir = 1, ivx, +1
       do is = 1, ivx, +1
         zw(:) = 0.d0
         do iq = 1, jvp, +1
          if((ir.gt.ivr).or.(is.gt.ivr).or.(iq.gt.jvpr))then
           zv2 = vp(is,iq)*vp(ir,iq) * dm(is,ir)
c           zv2 = gm(iq,is)*gm(iq,is) * cef1
           exc = ((eig(ir)+eig(is))*0.5d0-eigp(iq))*escale
           eta = alpha*abs(exc)
           if(eta.eq.0.d0)eta=T0*1.d-5
           do ie = 1, 2*ne, +1
            wtm = ene(ie) - exc
            weight = zv2/(wtm*wtm+eta*eta)
            zw(ie) = zw(ie) + weight * eta 
           enddo
          endif
         enddo

!$acc loop 
            do it = 1, ntime
!$acc loop 
           do ie = 1, 2*ne, +1
!$acc atomic update
            fdm_spectra1(it,ie)=fdm_spectra1(it,ie)+
     +                                       zw(ie)*dmc(it,ir,is)
c            fdm_spectra2(it,ie)=fdm_spectra2(it,ie)+
c     +                                       zw(ie)*dms(it,ir,is)
!$acc end atomic
            enddo
           enddo
       enddo
       enddo

c!$acc loop collapse(2) 
c!$acc& private(zw)
c      do ir = 1, ivx, +1
c       do is = 1, ivx, +1
c         zw(:) = 0.d0
c         do iq = 1, jvm, +1
c          if((ir.gt.ivr).or.(is.gt.ivr).or.(iq.gt.jvmr))then
c           zv2 = vm(is,iq)*vm(ir,iq) * dm(is,ir)
cc           zv2 = gp(iq,is)*gp(iq,ir) * cef1
c           exc = ((eig(ir)+eig(is))*0.5d0-eigm(iq))*escale
c           eta = alpha*abs(exc)
c           if(eta.eq.0.d0)eta=T0*1.d-5
c           do ie = 1, 2*ne, +1
c            wtm = ene(ie) - exc
c            weight = zv2/(wtm*wtm+eta*eta)
c            zw(ie) = zw(ie) + weight * eta 
c           enddo
c          endif
c         enddo
c
c!$acc loop 
c           do ie = 1, 2*ne, +1
c!$acc loop 
c            do it = 1, ntime
c!$acc atomic update
c            fdm_spectra1(it,ie)=fdm_spectra1(it,ie)+
c     +                                       zw(ie)*dmc(it,ir,is)
cc            fdm_spectra2(it,ie)=fdm_spectra2(it,ie)+
cc     +                                       zw(ie)*dms(it,ir,is)
c!$acc end atomic
c            enddo
c           enddo
c       enddo
c       enddo
c

c!$acc loop collapse(3) vector(32)
c       do iq = 1, jvp, +1
c      do ir = 1, ivx, +1
c         do is = 1, ivx, +1
cc          if((ir.gt.ivr).or.(is.gt.ivr).or.(iq.gt.jvpr))then
c           zv2 = vp(is,iq)*vp(ir,iq) * dm(is,ir)
cc           zv2 = gm(iq,is)*gm(iq,ir) * dm(is,ir)
c           exc = ((eig(ir)+eig(is))*0.5d0-eigp(iq))*escale
c           eta = alpha*abs(exc)
c!$acc loop collapse(2)
c           do it = 1, ntime
c            do ie = 1, 2*ne, +1
c            cdecay = 2.d0*xtime(it)*eta
c            decay=exp(-cdecay)
c             wtm = ene(ie) - exc
c             weight = zv2 * decay/(wtm*wtm+eta*eta)
c             zw1 = - weight * eta
c             zw2 =  weight * wtm
c!$acc atomic update
c            fdm_spectra1(it,ie)=fdm_spectra1(it,ie)+
c     +       (zw1*
c     +        (coset(ie,it)*dmpc(it,ir,iq)-sinet(ie,it)*dmps(it,ir,iq))+
c     +        zw2*
c     +        (sinet(ie,it)*dmpc(it,ir,iq)+coset(ie,it)*dmps(it,ir,iq)))
c!$acc end atomic
c            enddo
c           enddo
cc          endif
c         enddo
c        enddo
c       enddo

c
c!$acc loop   
c!$acc& private(zw1,zw2)
c      do ir = 1, ivx, +1
c        do iq = 1, jvm, +1
c         zw1(:,:) = 0.d0
c         zw2(:,:) = 0.d0
c         do is = 1, ivx, +1
c          if((ir.gt.ivr).or.(is.gt.ivr).or.(iq.gt.jvmr))then
c           zv2 = vm(is,iq)*vm(ir,iq) * dm(is,ir)
cc           zv2 = gp(iq,is)*gp(iq,ir) * dm(is,ir)
c           exc = ((eig(ir)+eig(is))*0.5d0-eigm(iq))*escale
c           eta = alpha*abs(exc)
c           if(eta.eq.0.d0)eta=T0*1.d-5
c           do it = 1, ntime
c           decay = 0.d0
c            cdecay = 2.d0*xtime(it)*eta
c            if(cdecay.le.36.d0)then
c             decay = 1.d0
c             if (cdecay.gt.4.d-16)decay=exp(-cdecay)
c            endif
c            do ie = 1, 2*ne, +1
c             wtm = ene(ie) - exc
c             weight = zv2/(wtm*wtm+eta*eta)
c             zw1(it,ie) = zw1(it,ie) - weight * eta * decay
c             zw2(it,ie) = zw2(it,ie) + weight * wtm * decay
c            enddo
c           enddo
c          endif
c         enddo
c
c!$acc loop 
c           do ie = 1, 2*ne, +1
c!$acc loop 
c           do it = 1, ntime
c!$acc atomic update
c            fdm_spectra1(it,ie)=fdm_spectra1(it,ie)+
c     +       (zw1(it,ie)*
c     +        (coset(ie,it)*dmmc(it,ir,iq)-sinet(ie,it)*dmms(it,ir,iq))+
c     +        zw2(it,ie)*
c     +        (sinet(ie,it)*dmmc(it,ir,iq)+coset(ie,it)*dmms(it,ir,iq)))
c!$acc end atomic
c           enddo
c           enddo
c        enddo
c       enddo
!$acc end kernels
!$acc exit data delete(dmc,dmmc,dmms,dmpc,dmps,xtime,ene,coset,sinet)
!$acc end data 

      return
      end

      subroutine broadening_t(ene,exc,weight,typ,alpha,spi,srpi,
     +                      eta)
      implicit none

      double precision ene, exc, weight, wtm, lgex, alpha, eta, spi,
     +                 omega0, srpi
      character*1 typ

      double precision  temp, mu
      common/temp_only/temp,mu


      if (typ.eq.'O') then
         wtm = (ene-exc)
         weight = 1.d0/(eta**2.d0+wtm**2.d0)!/(srpi**2.d0)
      endif

      if (typ.eq.'L') then
            wtm = exc/ene
            lgex = abs(wtm)
            lgex = log(lgex)/alpha
            lgex = lgex**2
            if (lgex.gt.37d0) then
               weight = 0.d0
            else
               wtm = 0.5d0*(1.0d0+sign(1.0d0,wtm))
               weight = wtm*spi*dexp(-lgex)/abs(exc)
            endif
      endif

      if (typ.eq.'G') then

        if (abs(ene).gt.5.d0*temp)then
         eta = 0.12d0*max(temp,abs(exc))*sqrt(2.d0)
        else
         eta = temp*sqrt(2.d0)
        endif

         wtm = (ene-exc)/eta
         lgex = wtm**2.d0
         if (lgex.gt.37d0)then
            weight = 0.d0
           else
            weight = dexp(-lgex)/srpi/eta
         endif
      endif

      return
      end
