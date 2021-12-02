      subroutine init_variables
      use inoutdata; use nrtype
      implicit none
      integer :: i,j
      include 'commands.inc'
      model=undefined_c
      input_weight=undefined_c
      dep_model=undefined_c
      probt=undefined
      lEXISTS=.TRUE.
      distr1=.FALSE.
      distr2=.FALSE.
      deprates=.FALSE.
      usr_z_dynpr=.FALSE.
      usr_z_c=.FALSE.
      usr_pcx_sol=.FALSE.
      only_deprates=.FALSE.
      pn_cut=.FALSE.
      ncomp=undefined_i
      mu=undefined
      dengas=undefined
      zlam=undefined
      zlams=undefined
      zlam_massive=undefined
      c0=undefined
      ks=undefined
      wtot_sample=undefined
      dens_ent=undefined
      dm_ent=undefined
      pnsavgguess=1.d0
      pnsmaxguess=1.d0
      pnsminguess=1.d0
      zsfavgguess=1.d0
      z0avgguess=0.0001d0
      z0minguess=0.0001d0
      z0maxguess=0.0001d0
      n_solutions=undefined_i
      do i=1,3
      rho_flow(i)=undefined
      ztot_flow(i)=undefined
      ush_flow(i)=undefined
      pns_flow(i)=undefined
      enddo
!      do i=1,2
!      phi50(i)=undefined
!      d50mm(i)=undefined
!      phi84(i)=undefined
!      d84mm(i)=undefined
!      phi16(i)=undefined
!      d16mm(i)=undefined
!      davgeqsph(i)=undefined
!      sorting(i)=undefined
!      enddo
      do i=1,6
      phi50(i)=undefined
      d50mm(i)=undefined
      phi84(i)=undefined
      d84mm(i)=undefined
      phi16(i)=undefined
      d16mm(i)=undefined
      davgeqsph(i)=undefined
      sorting(i)=undefined
      cdlaw(i)=undefined_c
      rholaw(i)=undefined_c
      rho_custom(i)=undefined_c
      fractal(i)=.FALSE.
      isometric(i)=.TRUE.
      merge_classes(i)=.FALSE.
      sensmerge(i)=undefined
      do j=0,30
      sphericity(i,j)=undefined
      longspher(i,j)=undefined
      crossspher(i,j)=undefined
      circularity(i,j)=undefined
      voleqsphd(i,j)=undefined
      circeqard(i,j)=undefined
      corey(i,j)=undefined
      voleqsphd(i,j)=undefined
      flatratio(i,j)=undefined
      shapefact(i,j)=undefined
      fractdim(i,j)=undefined
!      isometric(i,j)=.TRUE.
      shpar1read(i,j)=undefined
      shpar2read(i,j)=undefined
      shpar3read(i,j)=undefined
      rhos(i,j)=undefined
      enddo
      enddo
      do i=1,6
      dotestchi(i)=.TRUE.
      siglevchi(i)=undefined
      senschi(i)=0.05d0
      dphi(i)=undefined
      phimin(i)=undefined
      phimax(i)=undefined
      nclass(i)=undefined_i
      zdynpr(i)=undefined
      zc(i)=undefined
      pcx(i)=undefined
      enddo
      do i=1,6
      do j=1,30
      weight(i,j)=undefined
      enddo
      enddo
      slope_ground=undefined
      dep_median=undefined
      rhos_median=undefined
      end subroutine init_variables
