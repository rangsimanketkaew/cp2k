!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief sets the environment for optimization of exponents and contraction
!>        coefficients of the lri auxiliary
!>        lri : local resolution of the identity
!> \par History
!>      created Dorothea Golze [12.2014]
!> \authors Dorothea Golze
! **************************************************************************************************
MODULE lri_optimize_ri_basis_types

   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_optimize_ri_basis_types'
   PUBLIC :: lri_opt_type
   PUBLIC :: create_lri_opt, deallocate_lri_opt, get_original_gcc, &
             orthonormalize_gcc

! **************************************************************************************************

   TYPE lri_gcc_p_type
      ! gcc without normalization factor
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER                :: gcc_orig => NULL()
   END TYPE lri_gcc_p_type

   TYPE lri_subset_type
      ! amount of l quantum numbers per set
      INTEGER                                            :: nl = -1
      ! number of contraction per l quantum number for a given set
      INTEGER, DIMENSION(:), POINTER                     :: ncont_l => NULL()
   END TYPE lri_subset_type

   TYPE lri_opt_type
      LOGICAL                                            :: opt_exps = .FALSE.
      LOGICAL                                            :: opt_coeffs = .FALSE.
      LOGICAL                                            :: use_condition_number = .FALSE.
      LOGICAL                                            :: use_geometric_seq = .FALSE.
      LOGICAL                                            :: use_constraints = .FALSE.
      INTEGER                                            :: nexp = -1
      INTEGER                                            :: ncoeff = -1
      REAL(KIND=dp)                                      :: cond_weight = 0.0_dp
      REAL(KIND=dp)                                      :: scale_exp = 0.0_dp
      REAL(KIND=dp)                                      :: fermi_exp = 0.0_dp
      REAL(KIND=dp)                                      :: rho_diff = 0.0_dp
      ! array holding the variables that are optimized
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x => NULL()
      ! initial exponents
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zet_init => NULL()
      ! holds the original contraction coeff of the lri basis
      TYPE(lri_gcc_p_type), DIMENSION(:), POINTER        :: ri_gcc_orig => NULL()
      TYPE(lri_subset_type), DIMENSION(:), POINTER      :: subset => NULL()
   END TYPE lri_opt_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief creates lri_opt
!> \param lri_opt optimization environment
! **************************************************************************************************
   SUBROUTINE create_lri_opt(lri_opt)

      TYPE(lri_opt_type), POINTER                        :: lri_opt

      ALLOCATE (lri_opt)

      NULLIFY (lri_opt%ri_gcc_orig)
      NULLIFY (lri_opt%subset)
      NULLIFY (lri_opt%x)
      NULLIFY (lri_opt%zet_init)

      lri_opt%opt_exps = .FALSE.
      lri_opt%opt_coeffs = .FALSE.
      lri_opt%use_condition_number = .FALSE.
      lri_opt%use_geometric_seq = .FALSE.
      lri_opt%use_constraints = .FALSE.

      lri_opt%nexp = 0
      lri_opt%ncoeff = 0

   END SUBROUTINE create_lri_opt

! **************************************************************************************************
!> \brief deallocates lri_opt
!> \param lri_opt optimization environment
! **************************************************************************************************
   SUBROUTINE deallocate_lri_opt(lri_opt)

      TYPE(lri_opt_type), POINTER                        :: lri_opt

      INTEGER                                            :: i

      IF (ASSOCIATED(lri_opt)) THEN
         IF (ASSOCIATED(lri_opt%subset)) THEN
            DO i = 1, SIZE(lri_opt%subset)
               DEALLOCATE (lri_opt%subset(i)%ncont_l)
            END DO
            DEALLOCATE (lri_opt%subset)
         END IF
         IF (ASSOCIATED(lri_opt%x)) THEN
            DEALLOCATE (lri_opt%x)
         END IF
         IF (ASSOCIATED(lri_opt%zet_init)) THEN
            DEALLOCATE (lri_opt%zet_init)
         END IF
         IF (ASSOCIATED(lri_opt%ri_gcc_orig)) THEN
            DO i = 1, SIZE(lri_opt%ri_gcc_orig)
               DEALLOCATE (lri_opt%ri_gcc_orig(i)%gcc_orig)
            END DO
            DEALLOCATE (lri_opt%ri_gcc_orig)
         END IF
         DEALLOCATE (lri_opt)
      END IF
   END SUBROUTINE deallocate_lri_opt

! **************************************************************************************************
!> \brief primitive Cartesian Gaussian functions are normalized. The normalization
!>        factor is included in the Gaussian contraction coefficients.
!>        Division by this factor to get the original gcc.
!> \param gcc_orig original contraction coefficient
!> \param gto_basis_set gaussian type basis set
!> \param lri_opt optimization environment
! **************************************************************************************************
   SUBROUTINE get_original_gcc(gcc_orig, gto_basis_set, lri_opt)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc_orig
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis_set
      TYPE(lri_opt_type), POINTER                        :: lri_opt

      INTEGER                                            :: il, ipgf, iset, ishell, l, maxpgf, &
                                                            maxshell, nl, nset
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, ncont_l
      REAL(KIND=dp)                                      :: expzet, gcca, prefac, zeta

      maxpgf = SIZE(gto_basis_set%gcc, 1)
      maxshell = SIZE(gto_basis_set%gcc, 2)
      nset = SIZE(gto_basis_set%gcc, 3)

      ALLOCATE (gcc_orig(maxpgf, maxshell, nset))
      gcc_orig = 0.0_dp

      DO iset = 1, gto_basis_set%nset
         DO ishell = 1, gto_basis_set%nshell(iset)
            l = gto_basis_set%l(ishell, iset)
            expzet = 0.25_dp*REAL(2*l + 3, dp)
            prefac = 2.0_dp**l*(2.0_dp/pi)**0.75_dp
            DO ipgf = 1, gto_basis_set%npgf(iset)
               gcca = gto_basis_set%gcc(ipgf, ishell, iset)
               zeta = gto_basis_set%zet(ipgf, iset)
               gcc_orig(ipgf, ishell, iset) = gcca/(prefac*zeta**expzet)
            END DO
         END DO
      END DO

      IF (lri_opt%opt_coeffs) THEN
         ! **** get number of contractions per quantum number
         CALL get_gto_basis_set(gto_basis_set=gto_basis_set, &
                                lmax=lmax, lmin=lmin)
         ALLOCATE (lri_opt%subset(nset))
         DO iset = 1, gto_basis_set%nset
            nl = lmax(iset) - lmin(iset) + 1
            lri_opt%subset(iset)%nl = nl
            il = 1
            ALLOCATE (lri_opt%subset(iset)%ncont_l(nl))
            ncont_l => lri_opt%subset(iset)%ncont_l
            ncont_l = 1
            DO ishell = 2, gto_basis_set%nshell(iset)
               l = gto_basis_set%l(ishell, iset)
               IF (l == gto_basis_set%l(ishell - 1, iset)) THEN
                  ncont_l(il) = ncont_l(il) + 1
               ELSE
                  il = il + 1
                  ncont_l(il) = 1
               END IF
            END DO
         END DO
      END IF

   END SUBROUTINE get_original_gcc

! **************************************************************************************************
!> \brief orthonormalize contraction coefficients using Gram-Schmidt
!> \param gcc contraction coefficient
!> \param gto_basis_set gaussian type basis set
!> \param lri_opt optimization environment
! **************************************************************************************************
   SUBROUTINE orthonormalize_gcc(gcc, gto_basis_set, lri_opt)

      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis_set
      TYPE(lri_opt_type), POINTER                        :: lri_opt

      INTEGER                                            :: il, iset, ishell, ishell1, ishell2, &
                                                            istart, nset
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      REAL(KIND=dp)                                      :: gs_scale

      CALL get_gto_basis_set(gto_basis_set=gto_basis_set, nset=nset, nshell=nshell)

      DO iset = 1, nset
         istart = 1
         DO il = 1, lri_opt%subset(iset)%nl
            DO ishell1 = istart, istart + lri_opt%subset(iset)%ncont_l(il) - 2
               DO ishell2 = ishell1 + 1, istart + lri_opt%subset(iset)%ncont_l(il) - 1
                  gs_scale = DOT_PRODUCT(gcc(:, ishell2, iset), gcc(:, ishell1, iset))/ &
                             DOT_PRODUCT(gcc(:, ishell1, iset), gcc(:, ishell1, iset))
                  gcc(:, ishell2, iset) = gcc(:, ishell2, iset) - &
                                          gs_scale*gcc(:, ishell1, iset)
               END DO
            END DO
            istart = istart + lri_opt%subset(iset)%ncont_l(il)
         END DO

         DO ishell = 1, gto_basis_set%nshell(iset)
            gcc(:, ishell, iset) = gcc(:, ishell, iset)/ &
                                   SQRT(DOT_PRODUCT(gcc(:, ishell, iset), gcc(:, ishell, iset)))
         END DO
      END DO

   END SUBROUTINE orthonormalize_gcc

END MODULE lri_optimize_ri_basis_types
