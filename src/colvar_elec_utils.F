!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief evaluations of colvar for electronic properties
!> \par History
!>      09.2024 created [Rangsiman Ketkaew]
!> \author Rangsiman Ketkaew - University of Zurich (2024)
! **************************************************************************************************
MODULE colvar_elec_utils
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mult_fm_fm_fm, &
             mult_fm_sm_fm

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_elec_utils'

CONTAINS

! **************************************************************************************************
!> \brief Matrix multiplication of fm, sm, and fm for the first and second term in dqdr contribution
!> \param mat_A Full matrix
!> \param mat_B DBCSR sparse matrix
!> \param mat_C Full matrix
!> \param mat_tmp Temporary full matrix
!> \param product Product matrix
!> \param nao Number of AOs
!> \param nmo Number of MOs
!> \date  09.2024
!> \author Rangsiman Ketkaew
! **************************************************************************************************

   SUBROUTINE mult_fm_sm_fm(mat_A, mat_B, mat_C, mat_tmp, product, nao, nmo)

      TYPE(cp_fm_type)                                   :: mat_A, mat_C, mat_tmp, product
      TYPE(dbcsr_type)                                   :: mat_B
      INTEGER                                            :: nao, nmo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mult_fm_sm_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! A*B*C = A*(B*C)
      !       = A*(tmp)

      CALL cp_dbcsr_sm_fm_multiply(mat_B, mat_C, &
                                   mat_tmp, nmo)
      CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                         1.0_dp, mat_A, mat_tmp, &
                         0.0_dp, product)

      CALL timestop(handle)

   END SUBROUTINE mult_fm_sm_fm

! **************************************************************************************************
!> \brief Matrix multiplication of fm, fm, and fm for the third term in dqdr contribution
!> \param mat_A Full matrix
!> \param mat_B Full matrix
!> \param mat_C Full matrix
!> \param mat_tmp Temporary full matrix
!> \param product Product matrix
!> \param nao Number of AOs
!> \param nmo Number of MOs
!> \date  09.2024
!> \author Rangsiman Ketkaew
! **************************************************************************************************

   SUBROUTINE mult_fm_fm_fm(mat_A, mat_B, mat_C, mat_tmp, product, nao, nmo)

      TYPE(cp_fm_type)                                   :: mat_A, mat_B, mat_C, mat_tmp, product
      INTEGER                                            :: nao, nmo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mult_fm_fm_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! A*B*C = A*(B*C)
      !       = A*(tmp)

      CALL parallel_gemm("N", "N", nao, nmo, nao, &
                         1.0_dp, mat_B, mat_C, &
                         0.0_dp, mat_tmp)
      CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                         1.0_dp, mat_A, mat_tmp, &
                         0.0_dp, product)

      ! We can (probably) also convert fm to sm and use mult_fm_sm_fm:
      !  CALL dbcsr_set(dsdr_dbcsr, 0.0_dp)
      !  CALL copy_fm_to_dbcsr(dsdr(i,ii)%matrix, dsdr_dbcsr(i,ii)%matrix)
      !  CALL cp_dbcsr_sm_fm_multiply(dsdr_dbcsr(i,ii)%matrix, mo_coeff, tmp, nmo)

      CALL timestop(handle)

   END SUBROUTINE mult_fm_fm_fm

END MODULE colvar_elec_utils
