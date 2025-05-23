!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains the mapping ATOM_KIND -> ELEMENT for the most
!>      common cases in CHARMM and AMBER
!>      This should avoid in most cases the need to provide the element
!>      column in the PDB if in the atom_name column of the PDB is provided
!>      instead the atom kind
!> \par History
!>      10.2006 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
MODULE mm_mapping_library

   USE kinds,                           ONLY: default_string_length
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mm_mapping_library'

! **************************************************************************************************
   TYPE ff_map_type
      CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: kind => NULL()
      CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: element => NULL()
   END TYPE ff_map_type

   TYPE(ff_map_type), POINTER, PUBLIC :: amber_map, &
                                         charmm_map, &
                                         gromos_map
   PUBLIC :: create_ff_map, &
             destroy_ff_map
CONTAINS

! **************************************************************************************************
!> \brief Initialize arrays for mapping KINDS <-> ELEMENTS
!>      for major elements in AMBER and CHARMM
!> \param ff_type ...
!> \par History
!>      10.2006 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE create_ff_map(ff_type)
      CHARACTER(LEN=*), INTENT(IN)                       :: ff_type

      INTEGER, PARAMETER                                 :: amb_imax = 66, chm_imax = 108, &
                                                            grm_imax = 22

      SELECT CASE (ff_type)
      CASE ("AMBER")
         ALLOCATE (amber_map)
         ! allocate substructures
         ALLOCATE (amber_map%kind(amb_imax))
         ALLOCATE (amber_map%element(amb_imax))

         amber_map%kind = (/ &
                          "BR", "C ", "CA", "CB", "CC", "CD", "CI", "CK", "CP", "CM", &
                          "CS", "CN", "CQ", "CR", "CT", "CV", "CW", "C*", "CX", "CY", &
                          "CZ", "CO", "TG", "2C", "3C", "C0", "F ", "H ", "HC", "HA", &
                          "HO", "HS", "HW", "HP", "HZ", "I ", "IM", "IP", "IB", "MG", &
                          "N ", "NA", "NB", "NC", "N*", "NT", "NY", "ND", "NL", "TN", &
                          "OA", "O ", "OW", "OH", "OD", "OS", "OP", "P ", "S ", "SH", &
                          "CU", "FE", "Li", "K ", "Rb", "Cs"/)
         amber_map%element = (/ &
                             "Br", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", &
                             "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", &
                             "C ", "C ", "C ", "C ", "C ", "Ca", "F ", "H ", "H ", "H ", &
                             "H ", "H ", "H ", "H ", "H ", "I ", "Cl", "Na", "Na", "Mg", &
                             "N ", "N ", "N ", "N ", "N ", "N ", "N ", "N ", "N ", "N ", &
                             "O ", "O ", "O ", "O ", "O ", "O ", "O ", "P ", "S ", "S ", &
                             "Cu", "Fe", "Li", "K ", "Rb", "Cs"/)

      CASE ("CHARMM")
         ALLOCATE (charmm_map)
         ! allocate substructures
         ALLOCATE (charmm_map%kind(chm_imax))
         ALLOCATE (charmm_map%element(chm_imax))

         charmm_map%kind = (/ &
                           "CA  ", "CAL ", "C   ", "CC  ", "CD  ", "CE  ", "CEL ", "CES ", "CLA ", "CL  ", &
                           "CM  ", "CN1A", "CN1T", "CN3A", "CN3B", "CN3C", "CN3D", "CN3T", "CN5G", "CN7B", &
                           "CN7C", "CN7D", "CN8B", "CNA ", "CN  ", "CNE ", "CPA ", "CPB ", "CP  ", "CPH ", &
                           "CPM ", "CPT ", "CS  ", "CT  ", "CTL ", "CY  ", "DUM ", "FE  ", "FNA ", "FN  ", &
                           "HA  ", "HAL ", "HB  ", "HC  ", "HCL ", "HE  ", "HE  ", "HEL ", "H   ", "HL  ", &
                           "HN3B", "HN3C", "HNE ", "HN  ", "HNP ", "HOL ", "HP  ", "HR  ", "HS  ", "HT  ", &
                           "MG  ", "NC  ", "NE  ", "NH3L", "NH  ", "N   ", "NN1C", "NN2B", "NN2C", "NN2G", &
                           "NN2U", "NN3A", "NN3G", "NN3I", "NN  ", "NPH ", "NP  ", "NR  ", "NR  ", "NTL ", &
                           "NY  ", "O2L ", "OBL ", "OB  ", "OCL ", "OC  ", "OHL ", "OH  ", "OM  ", "ON1C", &
                           "ON2b", "ON6B", "ON  ", "O   ", "OSL ", "OS  ", "OT  ", "P2  ", "PL  ", "POT ", &
                           "P   ", "SL  ", "SM  ", "SOD ", "S   ", "SS  ", "ZN  ", "CN6B"/)

         charmm_map%element = (/ &
                              "C ", "Ca", "C ", "C ", "C ", "C ", "C ", "Cs", "Cl", "C ", &
                              "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", &
                              "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", "C ", &
                              "C ", "C ", "C ", "C ", "C ", "C ", "H ", "Fe", "F ", "F ", &
                              "H ", "H ", "H ", "H ", "H ", "H ", "He", "H ", "H ", "H ", &
                              "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", "H ", &
                              "Mg", "N ", "Ne", "N ", "N ", "N ", "N ", "N ", "N ", "N ", &
                              "N ", "N ", "N ", "N ", "N ", "N ", "N ", "N ", "N ", "N ", &
                              "N ", "O ", "O ", "O ", "O ", "O ", "O ", "O ", "O ", "O ", &
                              "O ", "O ", "O ", "O ", "O ", "O ", "O ", "P ", "P ", "K ", &
                              "P ", "S ", "S ", "Na", "S ", "S ", "Zn", "C "/)

      CASE ("GROMOS")
         ALLOCATE (gromos_map)
         ! allocate substructures
         ALLOCATE (gromos_map%kind(grm_imax))
         ALLOCATE (gromos_map%element(grm_imax))

         gromos_map%kind = (/ &
                           "H1  ", "CA  ", "HA  ", "SD  ", "OG  ", "CG  ", "HD  ", "SG  ", &
                           "CZ  ", "OH  ", "C2  ", "H8  ", "O4* ", "O5* ", "PA  ", "O1A ", &
                           "MG  ", "Na  ", "Cl  ", "OW  ", "HW  ", "HW  "/)

         gromos_map%element = (/ &
                              "H ", "C ", "H ", "S ", "O ", "C ", "H ", "S ", &
                              "C ", "O ", "C ", "H ", "O ", "O ", "P ", "O ", &
                              "Mg", "Na", "Cl", "O ", "H ", "H "/)

      END SELECT

   END SUBROUTINE create_ff_map

! **************************************************************************************************
!> \brief Deallocates the arrays used for mapping
!> \param ff_type ...
!> \par History
!>      10.2006 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE destroy_ff_map(ff_type)
      CHARACTER(LEN=*), INTENT(IN)                       :: ff_type

      SELECT CASE (ff_type)
      CASE ("AMBER")
         ! deallocate substructures
         DEALLOCATE (amber_map%kind)
         DEALLOCATE (amber_map%element)
         ! deallocate main
         DEALLOCATE (amber_map)
      CASE ("CHARMM")
         ! deallocate substructures
         DEALLOCATE (charmm_map%kind)
         DEALLOCATE (charmm_map%element)
         ! deallocate main
         DEALLOCATE (charmm_map)
      CASE ("GROMOS")
         ! deallocate substructures
         DEALLOCATE (gromos_map%kind)
         DEALLOCATE (gromos_map%element)
         ! deallocate main
         DEALLOCATE (gromos_map)
      END SELECT
   END SUBROUTINE destroy_ff_map

END MODULE mm_mapping_library
