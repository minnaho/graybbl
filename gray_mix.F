        module gray_mix
     !---------------------------------------------------------
     ! "Gray" bottom boundary layer for domains with small dz
     ! alternative to bottom KPP at high resolutions
     ! Akv = ustar*(z-zb)/(vonKar(1+((z-zb)/h)^4))
     ! h = n*dz, where n is number of grid points from the bottom
     ! h is the maximum between 1 m height above the bottom  
     ! (height is set in maxhbbl in gray_mix.opt)
     ! or number of grid points above bottom (np in gray_mix.opt) 
     !--------------------------------------------------------- 

#include "cppdefs.opt"
! CPP flag GRAY_MIX covering entire module, module empty if off
#if defined GRAY_MIX

        use param
        use hidden_mpi_vars
        use private_scratch

        use dimensions
        use ocean_vars
        use mixing
        use grid
        use coupling
        use scalars
        use random_output


        implicit none
#include "gray_mix.opt"
        private
        public gray_bmix
        contains

!-----------------------------------------------------------------

        subroutine gray_bmix(tind) ![
        implicit none
        integer :: i,j,k
        integer,intent(in) :: tind
        real,allocatable,dimension(:,:) :: ustar
        real,allocatable,dimension(:,:) :: hn
        !real hn

        allocate(ustar(GLOBAL_2D_ARRAY))
        ustar(:,:) = 0

        allocate(hn(GLOBAL_2D_ARRAY))
        hn(:,:) = 0

!#include "compute_extended_bounds.h"
        do j=0,ny
          do i=0,nx 
# ifdef MASKING
            if (rmask(i,j) > 0.5) then          
# endif
              hn(i,j) = max(maxhbbl,(z_w(i,j,np)+h(i,j))) ! parameters in .opt
              !hn = z_w(i,j,np)+h(i,j) ! parameters in .opt
              ustar(i,j)=sqrt(
     &                   sqrt((r_D(i,j)*(u(i-1,j,1,tind)+u(i,j,1,tind))/2)**2+
     &                        (r_D(i,j)*(v(i,j-1,1,tind)+v(i,j,1,tind))/2)**2))

              ustar_out(i,j) = ustar(i,j)
              hn_out(i,j) = hn(i,j)
              !rd_out(i,j) = r_D(i,j)

              do k=0,nz
                Akv(i,j,k) = ustar(i,j)*(z_w(i,j,k)+h(i,j))/
     &                       (vonKar*(1+((z_w(i,j,k)+h(i,j))/hn(i,j))**4))

                Akt(i,j,k,itemp) = Akv(i,j,k)

#  ifdef SALINITY
                Akt(i,j,k,isalt) = Akv(i,j,k)
#  endif 
              enddo
# ifdef MASKING
            else
              do k=0,nz
                Akv(i,j,k) = 0.
                Akt(i,j,k,itemp) = 0.
#  ifdef SALINITY
                Akt(i,j,k,isalt) = 0.
#  endif
              enddo
            endif
# endif
          enddo
        enddo

#  ifdef EXCHANGE
        call exchange_xxx(Akv(:,:,:),Akt(:,:,:,itemp))
#    ifdef SALINITY
        call exchange_xxx(Akt(:,:,:,isalt))
#    endif
#  endif

        end subroutine gray_bmix !]

# endif /* GRAY_MIX for whole module */

        end module gray_mix


