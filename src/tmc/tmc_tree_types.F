!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief module handles definition of the tree nodes for the global and
!>      the subtrees binary tree
!>                   parent element
!>                      /      \
!>      accepted (acc) /        \  not accepted (nacc)
!>                    /          \
!>                  child       child
!>                   / \         / \
!>
!>      tree creation assuming acceptance (acc) AND rejectance (nacc)
!>        of configuration
!>      if configuration is accepted: new configuration (child on acc) on basis
!>        of last configuration (one level up)
!>      if configuration is rejected: child on nacc on basis of last accepted
!>        element (last element which is on acc brach of its parent element)
!>      The global tree handles all configurations of different subtrees.
!>      The structure element "conf" is an array related to the temperature
!>        (sorted) and points to the subtree elements.
!> \par History
!>      11.2012 created [Mandes Schoenherr]
!> \author Mandes
! **************************************************************************************************

MODULE tmc_tree_types
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_tree_types'

   PUBLIC :: tree_type, global_tree_type
   PUBLIC :: elem_list_type, elem_array_type, gt_elem_list_type
   PUBLIC :: add_to_list, clean_list
   PUBLIC :: read_subtree_elem_unformated, write_subtree_elem_unformated

   !-- tree element status
   INTEGER, PARAMETER, PUBLIC :: status_created = 100
   INTEGER, PARAMETER, PUBLIC :: status_calculate_energy = 101
   INTEGER, PARAMETER, PUBLIC :: status_calc_approx_ener = 102

   INTEGER, PARAMETER, PUBLIC :: status_calculate_NMC_steps = 111
   INTEGER, PARAMETER, PUBLIC :: status_calculate_MD = 112
   INTEGER, PARAMETER, PUBLIC :: status_calculated = 113

   INTEGER, PARAMETER, PUBLIC :: status_accepted_result = 123
   INTEGER, PARAMETER, PUBLIC :: status_accepted = 122
   INTEGER, PARAMETER, PUBLIC :: status_rejected = 121
   INTEGER, PARAMETER, PUBLIC :: status_rejected_result = 120

   INTEGER, PARAMETER, PUBLIC :: status_cancel_nmc = 133
   INTEGER, PARAMETER, PUBLIC :: status_cancel_ener = 132
   INTEGER, PARAMETER, PUBLIC :: status_canceled_nmc = 131
   INTEGER, PARAMETER, PUBLIC :: status_canceled_ener = 130

   INTEGER, PARAMETER, PUBLIC :: status_deleted = 140
   INTEGER, PARAMETER, PUBLIC :: status_deleted_result = 141

   !-- dimension status (for e.g. dividing atoms in sub box)
   INTEGER, PARAMETER, PUBLIC :: status_ok = 42
   INTEGER, PARAMETER, PUBLIC :: status_frozen = -1
   INTEGER, PARAMETER, PUBLIC :: status_proton_disorder = 1

   !-- subtree element
   TYPE tree_type
      TYPE(tree_type), POINTER                :: parent => NULL() ! points to element one level up
      !-- acc..accepted goes to next level (next step),
      !   nacc..not accepted takes an alternative configutation
      TYPE(tree_type), POINTER                :: acc => NULL(), nacc => NULL()
      !-- type of MC move (swap is handled only in global tree)
      INTEGER                                  :: move_type = -1
      !-- status (e.g. calculated, MD calculation, accepted...)
      INTEGER                                  :: stat = status_created
      REAL(KIND=dp), DIMENSION(:), POINTER     :: subbox_center => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER     :: pos => NULL() ! position array
      INTEGER, DIMENSION(:), POINTER           :: mol => NULL() ! specifies the molecules the atoms participate
      REAL(KIND=dp), DIMENSION(:), POINTER     :: vel => NULL() ! velocity array
      REAL(KIND=dp), DIMENSION(:), POINTER     :: frc => NULL() ! force array
      REAL(KIND=dp), DIMENSION(:), POINTER     :: dipole => NULL() ! dipole moments array
      INTEGER, DIMENSION(:), POINTER           :: elem_stat => NULL() ! status for every dimension
      INTEGER                                  :: nr = -1 ! tree node number
      REAL(KIND=dp), DIMENSION(3, 2, 3)        :: rng_seed = 0 ! random seed for childs
      !-- remembers which subtree number element is from
      INTEGER                                  :: sub_tree_nr = -1
      !-- remembers the temperature the configurational change (NMC) is done with
      INTEGER                                  :: temp_created = 0
      !-- pointer to counter of next subtree element number
      INTEGER, POINTER                         :: next_elem_nr => NULL()
      !-- for calculating the NPT ensamble, variable box sizes are necessary.
      REAL(KIND=dp), DIMENSION(:), POINTER     :: box_scale => NULL()
      REAL(KIND=dp)                            :: potential = 0.0_dp ! potential energy
      !-- potential energy calculated using (MD potential) cp2k input file
      REAL(KIND=dp)                            :: e_pot_approx = 0.0_dp
      !-- kinetic energy (espacially for HMC, where the velocities are respected)
      REAL(KIND=dp)                            :: ekin = 0.0_dp
      !-- kinetic energy before md steps (after gaussian velocity change)
      REAL(KIND=dp)                            :: ekin_before_md = 0.0_dp
      !-- estimated energies are stored in loop order in this array
      REAL(KIND=dp), DIMENSION(4)              :: scf_energies = 0.0_dp
      !-- counter to get last position in the array loop
      INTEGER                                  :: scf_energies_count = 0
      !-- list of global tree elements referint to that node (reference back to global tree)
      !   if no reference exist anymore, global tree element can be deleted
      TYPE(gt_elem_list_type), POINTER         :: gt_nodes_references => NULL()
   END TYPE tree_type

   ! type for global tree element list in tree elements
   TYPE gt_elem_list_type
      TYPE(global_tree_type), POINTER         :: gt_elem => NULL()
      TYPE(gt_elem_list_type), POINTER        :: next => NULL()
   END TYPE gt_elem_list_type

   TYPE elem_list_type
      TYPE(tree_type), POINTER      :: elem => NULL()
      TYPE(elem_list_type), POINTER :: next => NULL()
      INTEGER                        :: temp_ind = 0
      INTEGER                        :: nr = -1
   END TYPE elem_list_type

   !-- array with subtree elements
   TYPE elem_array_type
      TYPE(tree_type), POINTER :: elem => NULL()
      LOGICAL                   :: busy = .FALSE.
      LOGICAL                   :: canceled = .FALSE.
      REAL(KIND=dp)             :: start_time = 0.0_dp
   END TYPE elem_array_type

   !-- global tree element
   TYPE global_tree_type
      TYPE(global_tree_type), POINTER :: parent => NULL() ! points to element one level up
      !-- acc..accepted goes to next level (next step),
      !   nacc..not accepted takes an alternative configutation
      TYPE(global_tree_type), POINTER :: acc => NULL(), nacc => NULL()
      !-- status (e.g. calculated, MD calculation, accepted...)
      INTEGER                                      :: stat = -99
      !-- remember if configuration in node are swaped
      LOGICAL                                      :: swaped = .FALSE.
      !-- stores the index of the configuration (temperature)
      !   which is changed
      INTEGER                                      :: mv_conf = -54321
      !-- stores the index of the configuration (temp.) which should change next
      INTEGER                                      :: mv_next_conf = -2345
      !-- list of pointes to subtree elements (Temp sorting)
      TYPE(elem_array_type), DIMENSION(:), ALLOCATABLE :: conf
      !-- remembers if last configuration is assumed to be accepted or rejected (next branc in tree);
      !   In case of swaping, it shows if the configuration of a certain temperature is assumed
      !   to be acc/rej (which branch is followed at the last modification of the conf of this temp.
      !TODO store conf_n_acc in a bitshifted array to decrease the size (1Logical = 1Byte)
      LOGICAL, DIMENSION(:), ALLOCATABLE           :: conf_n_acc
      INTEGER :: nr = 0 ! tree node number
      REAL(KIND=dp), DIMENSION(3, 2, 3)            :: rng_seed = 0.0_dp ! random seed for childs
      !-- random number for acceptance check
      REAL(KIND=dp)                                :: rnd_nr = 0.0_dp
      !-- approximate probability of acceptance will be adapted while calculating the exact energy
      REAL(KIND=dp)                                :: prob_acc = 0.0_dp ! estimated acceptance probability
      REAL(KIND=dp)                                :: Temp = 0.0_dp ! temperature for simulated annealing
   END TYPE global_tree_type

CONTAINS

! **************************************************************************************************
!> \brief add a certain element to the specified element list at the beginning
!> \param elem the sub tree element, to be added
!> \param list  ...
!> \param temp_ind ...
!> \param nr ...
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE add_to_list(elem, list, temp_ind, nr)
      TYPE(tree_type), POINTER                           :: elem
      TYPE(elem_list_type), POINTER                      :: list
      INTEGER, OPTIONAL                                  :: temp_ind, nr

      TYPE(elem_list_type), POINTER                      :: last, list_elem_tmp

      NULLIFY (list_elem_tmp, last)

      CPASSERT(ASSOCIATED(elem))

      ALLOCATE (list_elem_tmp)
      list_elem_tmp%elem => elem
      list_elem_tmp%next => NULL()
      IF (PRESENT(temp_ind)) THEN
         list_elem_tmp%temp_ind = temp_ind
      ELSE
         list_elem_tmp%temp_ind = -1
      END IF

      IF (PRESENT(nr)) THEN
         list_elem_tmp%nr = nr
      ELSE
         list_elem_tmp%nr = -1
      END IF

      IF (ASSOCIATED(list) .EQV. .FALSE.) THEN
         list => list_elem_tmp
      ELSE
         last => list
         DO WHILE (ASSOCIATED(last%next))
            last => last%next
         END DO
         last%next => list_elem_tmp
      END IF

   END SUBROUTINE add_to_list

! **************************************************************************************************
!> \brief clean a certain element element list
!> \param list  ...
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE clean_list(list)
      TYPE(elem_list_type), POINTER                      :: list

      TYPE(elem_list_type), POINTER                      :: list_elem_tmp

      NULLIFY (list_elem_tmp)

      DO WHILE (ASSOCIATED(list))
         list_elem_tmp => list%next
         DEALLOCATE (list)
         list => list_elem_tmp
      END DO
   END SUBROUTINE clean_list

! **************************************************************************************************
!> \brief prints out the TMC sub tree structure element unformated in file
!> \param elem ...
!> \param io_unit ...
!> \param
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE write_subtree_elem_unformated(elem, io_unit)
      TYPE(tree_type), POINTER                           :: elem
      INTEGER                                            :: io_unit

      CPASSERT(ASSOCIATED(elem))
      CPASSERT(io_unit .GT. 0)
      WRITE (io_unit) elem%nr, &
         elem%sub_tree_nr, &
         elem%stat, &
         elem%rng_seed, &
         elem%move_type, &
         elem%temp_created, &
         elem%potential, &
         elem%e_pot_approx, &
         elem%ekin, &
         elem%ekin_before_md
      CALL write_subtree_elem_darray(elem%pos, io_unit)
      CALL write_subtree_elem_darray(elem%vel, io_unit)
      CALL write_subtree_elem_darray(elem%frc, io_unit)
      CALL write_subtree_elem_darray(elem%box_scale, io_unit)
      CALL write_subtree_elem_darray(elem%dipole, io_unit)
   END SUBROUTINE write_subtree_elem_unformated

! **************************************************************************************************
!> \brief reads the TMC sub tree structure element unformated in file
!> \param elem ...
!> \param io_unit ...
!> \param
!> \author Mandes 11.2012
! **************************************************************************************************
   SUBROUTINE read_subtree_elem_unformated(elem, io_unit)
      TYPE(tree_type), POINTER                           :: elem
      INTEGER                                            :: io_unit

      CPASSERT(ASSOCIATED(elem))
      CPASSERT(io_unit .GT. 0)

      READ (io_unit) elem%nr, &
         elem%sub_tree_nr, &
         elem%stat, &
         elem%rng_seed, &
         elem%move_type, &
         elem%temp_created, &
         elem%potential, &
         elem%e_pot_approx, &
         elem%ekin, &
         elem%ekin_before_md
      CALL read_subtree_elem_darray(elem%pos, io_unit)
      CALL read_subtree_elem_darray(elem%vel, io_unit)
      CALL read_subtree_elem_darray(elem%frc, io_unit)
      CALL read_subtree_elem_darray(elem%box_scale, io_unit)
      CALL read_subtree_elem_darray(elem%dipole, io_unit)
   END SUBROUTINE read_subtree_elem_unformated

! **************************************************************************************************
!> \brief ...
!> \param array ...
!> \param io_unit ...
! **************************************************************************************************
   SUBROUTINE write_subtree_elem_darray(array, io_unit)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: array
      INTEGER                                            :: io_unit

      WRITE (io_unit) ASSOCIATED(array)
      IF (ASSOCIATED(array)) THEN
         WRITE (io_unit) SIZE(array)
         WRITE (io_unit) array
      END IF
   END SUBROUTINE write_subtree_elem_darray

! **************************************************************************************************
!> \brief ...
!> \param array ...
!> \param io_unit ...
! **************************************************************************************************
   SUBROUTINE read_subtree_elem_darray(array, io_unit)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: array
      INTEGER                                            :: io_unit

      INTEGER                                            :: i_tmp
      LOGICAL                                            :: l_tmp

      READ (io_unit) l_tmp
      IF (l_tmp) THEN
         READ (io_unit) i_tmp
         IF (ASSOCIATED(array)) THEN
            CPASSERT(SIZE(array) .EQ. i_tmp)
         ELSE
            ALLOCATE (array(i_tmp))
         END IF
         READ (io_unit) array
      END IF
   END SUBROUTINE read_subtree_elem_darray

END MODULE tmc_tree_types
