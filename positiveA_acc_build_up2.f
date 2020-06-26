c
      subroutine fdm_correlation_time_ave(escale,eig0,zmod,itype,nshell,
     +                           it,myid,nz,new_comm,T0)
      implicit none
      include 'mpif.h'
      include 'openacc_lib.h'
      !==========Start khai bao

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
      integer:: t1,t2, count_rate, count_max
      
      !build up
      double precision, dimension(1000,1000) ::dm_up,dmp_up,dmm_up
      double precision, dimension(1000,1000) ::vp_up,vm_up
      double precision, dimension(1000,1000) ::gp_up,gm_up
      double precision, dimension(1000) ::eig_up,eigp_up,eigm_up
      integer:: dm_up_size_x,dm_up_size_y,dmp_up_size_x,dmp_up_size_y
      integer:: dmm_up_size_x,dmm_up_size_y,vp_up_size_x,vp_up_size_y
      integer::vm_up_size_x,vm_up_size_y,gp_up_size_x,gp_up_size_y
      integer::gm_up_size_x,gm_up_size_y,eig_up_size_x,eigp_up_size_x
      integer::eigm_up_size_x,i,j
      double precision:: fp1,fm1,cleb
      integer::ir,is,iq
      double precision zw(2*ne)
      
      ! integer mdim
      ! parameter(mdim = 1280)
      double precision dmc(ntime,mdim,mdim), dmmc(ntime,mdim,mdim), 
     +      dmpc(ntime,mdim,mdim)
      double precision dms(ntime,mdim,mdim), dmms(ntime,mdim,mdim),
     +      dmps(ntime,mdim,mdim)

      !==========End khai bao
      common/temp_only/temp,mu


      !start log time
      call system_clock(count_max=count_max, count_rate=count_rate)
      call system_clock(t1)

      

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

      !===========MAIN LOOP============
      do ishell = 24, 24, +1
         ! t_x = 0
         ! t_m = 0
         ! t_p = 0
         dm_up_size_x=0
         dm_up_size_y=0
         dmp_up_size_x=0
         dmp_up_size_y=0
         dmm_up_size_x=0
         dmm_up_size_y=0
         vp_up_size_x=0
         vp_up_size_y=0
         vm_up_size_x=0
         vm_up_size_y=0
         gp_up_size_x=0
         gp_up_size_y=0
         gm_up_size_x=0
         gm_up_size_y=0
         eig_up_size_x=0
         eigp_up_size_x=0
         eigm_up_size_x=0
         
         !update element in matrix dm_up
         call set_zeros_2D(dm_up)
         call set_zeros_2D(dmp_up)
         call set_zeros_2D(dmm_up)
         call set_zeros_2D(vp_up)
         call set_zeros_2D(vm_up)
         call set_zeros_2D(gp_up)
         call set_zeros_2D(gm_up)
         
         
         !update element in matrix eig_up
         call set_zeros_1D(eig_up)
         call set_zeros_1D(eigp_up)
         call set_zeros_1D(eigm_up)
         
         minae = 0
         maxae = 2*ishell

         do ielec = minae, maxae
            minas = max(mod(ielec,2),0)
            maxas = min(2 * ishell - ielec, ielec)

            do ispin = minas, maxas, +2
               ivx = ns(ishell,ielec,ispin)
               jvp = ns(ishell,ielec-1,ispin+1)
               jvm = ns(ishell,ielec-1,ispin-1)

               ! t_x = t_x + ivx
               ! t_p = t_p + jvp
               ! t_m = t_m + jvm

               if (ivx.gt.0) then
                  if ((jvp.gt.0).or.(jvm.gt.0)) then
                  ivr = nsr(ishell,ielec,ispin)
                  jvpr = nsr(ishell,ielec-1,ispin+1)
                  jvmr = nsr(ishell,ielec-1,ispin-1)

                  call sub_fdm_correlation(dm,dmp,dmm,eig,eigp,eigm,vp,
     +                 vm,gp,gm,ivx,jvp,jvm,ivr,jvpr,jvmr,ishell,ielec,
     +                 ispin,escale(ishell),wz(ishell),eig0(ishell),
     +                 zmod(ishell),alpha,omega0,srpi,fdm_spectra1,myid,
     +            fdm_spectra2,itype,ntime,ne,T0,coset,sinet,
     +fp1, fm1, cleb)

                  
                  do i = 1, ivx
                     do j = 1, ivx
                        dm(i,j) = dm(i,j)*cleb
                     enddo
                  enddo
                  !build up matrix
                  
                  !update element in matrix dm_up
                  do i=dm_up_size_x+1,dm_up_size_x+ivx
                     do j=dm_up_size_y+1,dm_up_size_y+ivx
                           dm_up(i,j)=dm(i-dm_up_size_x,j-dm_up_size_y)
                     end do
                  end do
                  dm_up_size_x=dm_up_size_x+ivx
                  dm_up_size_y=dm_up_size_y+ivx

                  !update element in matrix dmp_up
                  do i=dmp_up_size_x+1,dmp_up_size_x+jvp
                     do j=dmp_up_size_y+1,dmp_up_size_y+jvp
                           dmp_up(i,j)=dmp(i-dmp_up_size_x,
     +                                     j-dmp_up_size_y)
                     end do
                  end do
                  dmp_up_size_x=dmp_up_size_x+jvp
                  dmp_up_size_y=dmp_up_size_y+jvp
                  

                  !update element in matrix dmm_up
                  call update_large_matrix(dmm_up_size_x, jvm, 
     +    dmm_up_size_y, jvm,dmm,dmm_up)

                  !update element in matrix vp_up
                  call update_large_matrix(vp_up_size_x, ivx, 
     +    vp_up_size_y, jvp,vp,vp_up)
                  
                  !update element in matrix vm_up
                  call update_large_matrix(vm_up_size_x, ivx
     +    , vm_up_size_y, jvm, vm,vm_up)
                  
                  !update element in matrix gp_up
                  call update_large_matrix(gp_up_size_x, jvm
     +    , gp_up_size_y, ivx, gp,gp_up)
                  
                  
                  !update element in matrix gm_up
                  call update_large_matrix(gm_up_size_x, jvp
     +    , gm_up_size_y, ivx, gm,gm_up)
                  
                  !update element in matrix eig_up
                  do i=eig_up_size_x+1,eig_up_size_x+ivx
                     eig_up(i)=eig(i-eig_up_size_x)
                  end do
                  eig_up_size_x=eig_up_size_x+ivx
                  !update element in matrix eigp_up
                  do i=eigp_up_size_x+1,eigp_up_size_x+jvp
                     eigp_up(i)=eigp(i-eigp_up_size_x)
                  end do
                  eigp_up_size_x=eigp_up_size_x+jvp
                  !update element in matrix eigm_up
                  do i=eigm_up_size_x+1,eigm_up_size_x+jvm
                     eigm_up(i)=eigm(i-eigm_up_size_x)
                  end do
                  eigm_up_size_x=eigm_up_size_x+jvm
                  ! end build up matrix
                  
                  endif
               endif

            enddo
         enddo
         !voi moi m
         ! can xac dinh dau la B,C,E
         !Build up ma tran to hon
         ! write(*,*) ishell, t_x, t_p, t_m, n_m
         ivx=dm_up_size_x
         jvp=dmp_up_size_x
         jvm=dmm_up_size_x
         ! fp1,fm1: ko co vai tro gi ca
         !chi co 1 vong lap cua m tai 24 => escale
         !cleb khong co vai tro gi
         !alpha,omega0,srpi,ntime la hang so
         !ivr = nsr(ishell,ielec,ispin)
         
         call commutator(dm_up,dmp_up,dmm_up,eig_up,eigp_up,eigm_up,
     +                   vp_up,vm_up,gp_up,gm_up,fp1,fm1,
     +                   ivx,jvp,jvm,escale(24),cleb,alpha,omega0,srpi,
     +                   fdm_spectra1,fdm_spectra2,ntime,dmc,dmmc,dmpc,
     +                   dms,dmms,dmps,ivr,jvpr,jvmr,ne,T0,coset,sinet,
     +                   myid)
      enddo
      !end log time
      call system_clock(t2)
      ! print *,'total time', real(t2-t1)/real(count_rate),'seconds'
      write(*,*) 'total time'
      write(*,*) real(t2-t1)/real(count_rate),'seconds'

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
     +       fdm_spectra2,itype0,ntime,ne,T0,coset,sinet,
     + fp1,fm1,cleb) !nghia khai bao them
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

   !    integer mdim
   !    parameter(mdim = 1280)
   !    double precision dmc(ntime,mdim,mdim), dmmc(ntime,mdim,mdim), 
   !   +      dmpc(ntime,mdim,mdim)
   !    double precision dms(ntime,mdim,mdim), dmms(ntime,mdim,mdim),
   !   +      dmps(ntime,mdim,mdim)

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

   !    call commutator(dm,dmp,dmm,eig,eigp,eigm,vp,vm,gp,gm,fp1,fm1,
   !   +  ivx,jvp,jvm,escale,cleb,alpha,omega0,srpi,fdm_spectra1,
   !   +  fdm_spectra2,ntime,dmc,dmmc,dmpc,dms,dmms,dmps,
   !   +  ivr,jvpr,jvmr,ne,T0,coset,sinet,myid)
   
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

      ! do ir = 1, ivx
      ! do is = 1, ivx
      ! dm(ir,is) = dm(ir,is)*cleb
      ! enddo
      ! enddo

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
         !  if((ir.gt.ivr).or.(is.gt.ivr).or.(iq.gt.jvpr))then
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
         !  endif
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

      subroutine set_zeros_2D(A)
         double precision, dimension(1000,1000) ::A
         do i=1,1000
            do j=1,1000
                  A(i,j)=0
            end do
         end do
      end

      subroutine set_zeros_1D(A)
         double precision, dimension(1000) ::A
         do i=1,1000
            A(i)=0
         end do
      end

      subroutine update_large_matrix(size_x, delta_x, size_y, delta_y
     +    ,A,A_up)
         integer:: size_x,size_y
         double precision, dimension(1000,1000) ::A, A_up
         
         do i=size_x+1,size_x+delta_x
            do j=size_y+1,size_y+delta_x
               A_up(i,j)=A(i-size_x,j-size_y)
            end do
         end do
         size_x=size_x+delta_x
         size_y=size_y+delta_x

      end