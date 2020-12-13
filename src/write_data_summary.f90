      subroutine write_data_summary
      USE inoutdata; USE nrtype
      implicit none
      if(only_deprates) return
      write(*,*)'Data summary'
                if(model.eq.'TWOLAYERS') then
                write(*,*)'*********TWO LAYERS METHOD*********'
                write(*,*)'ENTRAINED PARTICLE'
                write(*,190)dens_ent,dm_ent
                write(*,*)'REPRESENTATIVE PARTICLE IN THE OVERLYING LAYER'
                write(*,191)rhos(1,0),d50mm(1),sorting(1),nclass(1),cdlaw(1)
                        select case (cdlaw(1))
                        case ('HAIDLEV')
                        write(*,192)shpar1(1,0)
                        case ('SWAMOJ')
                        write(*,193)shpar1(1,0)
                        case ('GANSER')
!                        write(*,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1,0)
                        write(*,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1)
                        case ('CHIEN')
                        write(*,192)shpar1(1,0)
                        case ('TRANCONG')
                        write(*,195)shpar1(1,0),shpar2(1,0)
                        case ('DELLINO')
                        write(*,196)shpar1(1,0)
                        case ('HOLZSOMM')
                        write(*,197)shpar1(1,0),shpar2(1,0),shpar3(1,0)
                        case ('DIOGMELE')
                        write(*,196)shpar1(1,0)
                        case ('DIOG2016')
                             if(fractal(1)) then
                             write(*,199)shpar1(1,0)
                             else
                             write(*,192)shpar1(1,0)
                             endif
                        end select
                write(*,*)'OTHER DATA'
                write(*,198)probt,zlam,zlams,c0,ks
                else
                write(*,*)'*********TWO COMPONENTS METHOD*********'
                write(*,*)'COMPONENT 1'
                write(*,191)rhos(1,0),d50mm(1),sorting(1),nclass(1),cdlaw(1)
                        select case (cdlaw(1))
                        case ('HAIDLEV')
                        write(*,192)shpar1(1,0)
                        case ('SWAMOJ')
                        write(*,193)shpar1(1,0)
                        case ('GANSER')
!                        write(*,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1,0)
                        write(*,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1)
                        case ('CHIEN')
                        write(*,192)shpar1(1,0)
                        case ('TRANCONG')
                        write(*,195)shpar1(1,0),shpar2(1,0)
                        case ('DELLINO')
                        write(*,196)shpar1(1,0)
                        case ('HOLZSOMM')
                        write(*,197)shpar1(1,0),shpar2(1,0),shpar3(1,0)
                        case ('DIOGMELE')
                        write(*,196)shpar1(1,0)
                        case ('DIOG2016')
                             if(fractal(1)) then
                             write(*,199)shpar1(1,0)
                             else
                             write(*,192)shpar1(1,0)
                             endif
                        end select
                write(*,*)'COMPONENT 2'
                write(*,191)rhos(2,0),d50mm(2),sorting(2),nclass(2),cdlaw(2)
                        select case (cdlaw(2))
                        case ('HAIDLEV')
                        write(*,192)shpar1(2,0)
                        case ('SWAMOJ')
                        write(*,193)shpar1(2,0)!,dlong(2,0),dmed(2,0),dshort(2,0)
                        case ('GANSER')
!                        write(*,194)shpar1(2,0),shpar2(2,0),shpar3(2,0),shpar4(2,0)
                        write(*,194)shpar1(2,0),shpar2(2,0),shpar3(2,0),shpar4(2)
                        case ('CHIEN')
                        write(*,192)shpar1(2,0)
                        case ('TRANCONG')
                        write(*,195)shpar1(2,0),shpar2(2,0)
                        case ('DELLINO')
                        write(*,196)shpar1(2,0)
                        case ('HOLZSOMM')
                        write(*,197)shpar1(2,0),shpar2(2,0),shpar3(2,0)
                        case ('DIOGMELE')
                        write(*,196)shpar1(2,0)
                        case ('DIOG2016')
                             if(fractal(2)) then
                             write(*,199)shpar1(2,0)
                             else
                             write(*,192)shpar1(2,0)
                             endif
                        end select
                write(*,*)'OTHER DATA'
                write(*,198)probt,zlam,zlams,c0,ks
                endif

                write(50,*)'Data summary'
                    if(model.eq.'TWOLAYERS') then
                    write(50,*)'*********TWO LAYER METHOD*********'
                    write(50,*)'ENTRAINED PARTICLE'
                    write(50,190)dens_ent,dm_ent
                    write(50,*)'REPRESENTATIVE PARTICLE IN THE OVERLYING LAYER'
                    write(50,191)rhos(1,0),d50mm(1),sorting(1),nclass(1),cdlaw(1)
                        select case (cdlaw(1))
                        case ('HAIDLEV')
                        write(50,192)shpar1(1,0)
                        case ('SWAMOJ')
                        write(50,193)shpar1(1,0)!,dlong(1,0),dmed(1,0),dshort(1,0)
                        case ('GANSER')
!                        write(50,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1,0)
                        write(50,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1)
                        case ('CHIEN')
                        write(50,192)shpar1(1,0)
                        case ('TRANCONG')
                        write(50,195)shpar1(1,0),shpar2(1,0)
                        case ('DELLINO')
                        write(50,196)shpar1(1,0)
                        case ('HOLZSOMM')
                        write(50,197)shpar1(1,0),shpar2(1,0),shpar3(1,0)
                        case ('DIOGMELE')
                        write(50,196)shpar1(1,0)
                        case ('DIOG2016')
                             if(fractal(1)) then
                             write(50,199)shpar1(1,0)
                             else
                             write(50,192)shpar1(1,0)
                             endif
                        end select
                    write(50,*)'OTHER DATA'
                    write(50,198)probt,zlam,zlams,c0,ks
                    else
                    write(50,*)'*********TWO COMPONENT METHOD*********'
                    write(50,*)'COMPONENT 1'
                    write(50,191)rhos(1,0),d50mm(1),sorting(1),nclass(1),cdlaw(1)
                        select case (cdlaw(1))
                        case ('HAIDLEV')
                        write(50,192)shpar1(1,0)
                        case ('SWAMOJ')
                        write(50,193)shpar1(1,0)!,dlong(1,0),dmed(1,0),dshort(1,0)
                        case ('GANSER')
!                        write(50,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1,0)
                        write(50,194)shpar1(1,0),shpar2(1,0),shpar3(1,0),shpar4(1)
                        case ('CHIEN')
                        write(50,192)shpar1(1,0)
                        case ('TRANCONG')
                        write(50,195)shpar1(1,0),shpar2(1,0)
                        case ('DELLINO')
                        write(50,196)shpar1(1,0)
                        case ('HOLZSOMM')
                        write(50,197)shpar1(1,0),shpar2(1,0),shpar3(1,0)
                        case ('DIOGMELE')
                        write(50,196)shpar1(1,0)
                        case ('DIOG2016')
                             if(fractal(1)) then
                             write(50,199)shpar1(1,0)
                             else
                             write(50,192)shpar1(1,0)
                             endif
                        end select
                    write(50,*)'COMPONENT 2'
                    write(50,191)rhos(2,0),d50mm(2),sorting(2),nclass(2),cdlaw(2)
                        select case (cdlaw(2))
                        case ('HAIDLEV')
                        write(50,192)shpar1(2,0)
                        case ('SWAMOJ')
                        write(50,193)shpar1(2,0)!,dlong(2,0),dmed(2,0),dshort(2,0)
                        case ('GANSER')
!                        write(50,194)shpar1(2,0),shpar2(2,0),shpar3(2,0),shpar4(2,0)
                        write(50,194)shpar1(2,0),shpar2(2,0),shpar3(2,0),shpar4(2)
                        case ('CHIEN')
                        write(50,192)shpar1(2,0)
                        case ('TRANCONG')
                        write(50,195)shpar1(2,0),shpar2(2,0)
                        case ('DELLINO')
                        write(50,196)shpar1(2,0)
                        case ('HOLZSOMM')
                        write(50,197)shpar1(2,0),shpar2(2,0),shpar3(2,0)
                        case ('DIOGMELE')
                        write(50,196)shpar1(2,0)
                        case ('DIOG2016')
                             if(fractal(2)) then
                             write(50,199)shpar1(2,0)
                             else
                             write(50,192)shpar1(2,0)
                             endif
                        end select
                    write(50,*)'OTHER DATA'
                    write(50,198)probt,zlam,zlams,c0,ks,dep_median,rhos_median
                    if(slope_ground.ne.undefined) write(50,200) slope_ground
                    endif
  190 format('Density of the entrained particle (kg/m^3)',7x,f10.3,/,&
     &'Diameter of the entrained particle (m)',11x,f10.3,/)
  191 format('Particle density (kg/m^3)',25x,f10.3,/,&
     &'Particle equivalent diameter of the median size (mm)',f8.3,/,&
     &'Sorting particle grainsize distribution (phi)',7x,f8.3,/,&
     &'Classes number particle grainsize distribution (-)',3x,i4,/,&
     &'Drag law ',45x,a10)
  192 format('Sphericity (-)',38x,f8.3,/)
  193 format('Corey shape factor (-)',32x,f8.3,/)
  194 format('Sphericity (-)',38x,f8.3,/,&
     &'Volume equivalent sphere diameter (mm)',14x,f8.3,/,&
     &'Equal projected area circle diameter (mm)',11x,f8.3,/,&
     &'Isometric particle',34x,l8,/)
  195 format('Circularity (-)',38x,f8.3,/,&
      'Flatness ratio (-)',36x,f8.3,/)
  196 format('Shape factor (-)',36x,f8.3,/)
  197 format('Sphericity (-)',40x,f8.3,/,&
     &'Lengthwise sphericity (-)',29x,f8.3,/,&
     &'Crosswise sphericity (-)',30x,f8.3,/)
  198 format('Significance level t-test (-)',23x,f8.3,/,&
     &'Layer thickness (m)',35x,f8.5,/,&
     &'Sublayer thickness (m)',32x,f8.5,/,&
     &'Particle concentration in the layer (-)',13x,f8.3,/,&
     &'Substrate roughness (m)',31x,f8.5,/,&
     &'Deposit median (m)',36x,f8.5,/,&
     &'Deposit density (kg m-3)',30x,f8.3)
  199 format('Fractal dimension (-)',38x,f8.3,/)
  200 format('Slope ground (deg)',34x,f8.3)
      end subroutine write_data_summary
