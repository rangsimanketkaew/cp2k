!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief module that contains the algorithms to perform an itrative
!>         diagonalization by the block-Davidson approach
!>         P. Blaha, et al J. Comp. Physics, 229, (2010), 453-460
!>         \Iterative diagonalization in augmented plane wave based
!>              methods in electronic structure calculations\
!> \par History
!>      05.2011 created [MI]
!> \author MI
! **************************************************************************************************
MODULE qs_block_davidson_types

   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_block_davidson_types'

   PUBLIC :: block_davidson_allocate, block_davidson_deallocate, block_davidson_release, &
             block_davidson_env_create, davidson_type

   TYPE davidson_type
      INTEGER :: max_iter = -1, prec_type = -1, solver_type = -1, niter_new_prec = -1, first_prec = -1
      LOGICAL :: use_sparse_mos = .FALSE.
      REAL(KIND=dp) :: conv_percent = -1.0_dp, energy_gap = -1.0_dp, eps_iter = -1.0_dp
      TYPE(cp_fm_type), POINTER :: H_block_mat => NULL(), H_block_vec => NULL(), &
                                   matrix_z => NULL(), matrix_pz => NULL(), S_block_mat => NULL(), W_block_mat => NULL()
   END TYPE davidson_type

CONTAINS

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param bdav_env ...
!> \param nspins ...
!> \param scf_section ...
! **************************************************************************************************
   SUBROUTINE block_davidson_env_create(bdav_env, nspins, scf_section)

      TYPE(davidson_type), DIMENSION(:), POINTER         :: bdav_env
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(section_vals_type), POINTER                   :: scf_section

      INTEGER                                            :: ispin

      CPASSERT(.NOT. ASSOCIATED(bdav_env))
      ALLOCATE (bdav_env(nspins))
      DO ispin = 1, nspins
         NULLIFY (bdav_env(ispin)%H_block_mat)
         NULLIFY (bdav_env(ispin)%H_block_vec)
         NULLIFY (bdav_env(ispin)%S_block_mat)
         NULLIFY (bdav_env(ispin)%W_block_mat)
         NULLIFY (bdav_env(ispin)%matrix_z)
         NULLIFY (bdav_env(ispin)%matrix_pz)

         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DAVIDSON%PRECONDITIONER", &
                                   i_val=bdav_env(ispin)%prec_type)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DAVIDSON%PRECOND_SOLVER", &
                                   i_val=bdav_env(ispin)%solver_type)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DAVIDSON%ENERGY_GAP", &
                                   r_val=bdav_env(ispin)%energy_gap)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DAVIDSON%NEW_PREC_EACH", &
                                   i_val=bdav_env(ispin)%niter_new_prec)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%MAX_ITER", &
                                   i_val=bdav_env(ispin)%max_iter)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%EPS_ITER", &
                                   r_val=bdav_env(ispin)%eps_iter)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DAVIDSON%FIRST_PREC", &
                                   i_val=bdav_env(ispin)%first_prec)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DAVIDSON%CONV_MOS_PERCENT", &
                                   r_val=bdav_env(ispin)%conv_percent)
         CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DAVIDSON%SPARSE_MOS", &
                                   l_val=bdav_env(ispin)%use_sparse_mos)

      END DO

   END SUBROUTINE block_davidson_env_create

! **************************************************************************************************
!> \brief ...
!> \param bdav_env ...
!> \param mo_coeff ...
!> \param nao ...
!> \param nmo ...
! **************************************************************************************************
   SUBROUTINE block_davidson_allocate(bdav_env, mo_coeff, nao, nmo)

      TYPE(davidson_type)                                :: bdav_env
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff
      INTEGER, INTENT(IN)                                :: nao, nmo

      CHARACTER(len=*), PARAMETER :: routineN = 'block_davidson_allocate'

      INTEGER                                            :: handle, nmox2
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp

      CALL timeset(routineN, handle)
      NULLIFY (fm_struct_tmp)

      nmox2 = 2*nmo

      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmo, &
                               para_env=mo_coeff%matrix_struct%para_env, &
                               context=mo_coeff%matrix_struct%context)
      NULLIFY (bdav_env%matrix_z, bdav_env%matrix_pz)
      ALLOCATE (bdav_env%matrix_z, bdav_env%matrix_pz)
      CALL cp_fm_create(bdav_env%matrix_z, fm_struct_tmp, name="Z_mat")
      CALL cp_fm_create(bdav_env%matrix_pz, fm_struct_tmp, name="Z_mat")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL timestop(handle)

   END SUBROUTINE block_davidson_allocate

! **************************************************************************************************
!> \brief ...
!> \param bdav_env ...
! **************************************************************************************************
   SUBROUTINE block_davidson_deallocate(bdav_env)

      TYPE(davidson_type), DIMENSION(:), POINTER         :: bdav_env

      INTEGER                                            :: ispin, nspins

      IF (ASSOCIATED(bdav_env)) THEN

         nspins = SIZE(bdav_env)
         DO ispin = 1, nspins

            CALL cp_fm_release(bdav_env(ispin)%matrix_z)
            CALL cp_fm_release(bdav_env(ispin)%matrix_pz)
            DEALLOCATE (bdav_env(ispin)%matrix_z, bdav_env(ispin)%matrix_pz)
            NULLIFY (bdav_env(ispin)%matrix_z, bdav_env(ispin)%matrix_pz)

         END DO

      END IF

   END SUBROUTINE block_davidson_deallocate

! **************************************************************************************************
!> \brief ...
!> \param bdav_env ...
! **************************************************************************************************
   SUBROUTINE block_davidson_release(bdav_env)

      TYPE(davidson_type), DIMENSION(:), POINTER         :: bdav_env

      INTEGER                                            :: ispin, nspins

      IF (ASSOCIATED(bdav_env)) THEN

         nspins = SIZE(bdav_env)
         DO ispin = 1, nspins

            IF (ASSOCIATED(bdav_env(ispin)%matrix_z)) THEN
               CALL cp_fm_release(bdav_env(ispin)%matrix_z)
               CALL cp_fm_release(bdav_env(ispin)%matrix_pz)
               DEALLOCATE (bdav_env(ispin)%matrix_z, bdav_env(ispin)%matrix_pz)
               NULLIFY (bdav_env(ispin)%matrix_z, bdav_env(ispin)%matrix_pz)
            END IF

         END DO
         DEALLOCATE (bdav_env)

      END IF

   END SUBROUTINE block_davidson_release

END MODULE qs_block_davidson_types
