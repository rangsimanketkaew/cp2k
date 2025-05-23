!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief builds the input structure for the FORCE_EVAL section of cp2k
!> \par History
!>      06.2004 created [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_force_eval
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              debug_print_level,&
                                              high_print_level,&
                                              low_print_level,&
                                              medium_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        do_eip, do_embed, do_fist, do_ipi, do_mixed, do_nnp, do_qmmm, do_qs, do_sirius, &
        do_stress_analytical, do_stress_diagonal_anal, do_stress_diagonal_numer, do_stress_none, &
        do_stress_numerical, numerical
   USE input_cp2k_dft,                  ONLY: create_bsse_section,&
                                              create_dft_section
   USE input_cp2k_eip,                  ONLY: create_eip_section
   USE input_cp2k_embed,                ONLY: create_embed_section
   USE input_cp2k_mixed,                ONLY: create_mix_section
   USE input_cp2k_mm,                   ONLY: create_mm_section
   USE input_cp2k_nnp,                  ONLY: create_nnp_section
   USE input_cp2k_properties_dft,       ONLY: create_properties_section
   USE input_cp2k_pwdft,                ONLY: create_pwdft_section
   USE input_cp2k_qmmm,                 ONLY: create_qmmm_section
   USE input_cp2k_subsys,               ONLY: create_subsys_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              lchar_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_force_eval'

   PUBLIC :: create_force_eval_section

CONTAINS

! **************************************************************************************************
!> \brief creates the force_eval section
!> \param section the section to be created
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_force_eval_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="force_eval", &
                          description="parameters needed to calculate energy and forces and"// &
                          " describe the system you want to analyze.", &
                          n_keywords=1, n_subsections=10, repeats=.TRUE.)

      NULLIFY (subsection)
      NULLIFY (keyword)
      CALL keyword_create(keyword, __LOCATION__, name="METHOD", &
                          description="Which method should be used to compute forces", &
                          usage="METHOD <STRING>", &
                          enum_c_vals=s2a("QS", &
                                          "SIRIUS", &
                                          "FIST", &
                                          "QMMM", &
                                          "EIP", &
                                          "QUICKSTEP", &
                                          "NNP", &
                                          "MIXED", &
                                          "EMBED", &
                                          "IPI"), &
                          enum_desc=s2a("Alias for QUICKSTEP", &
                                        "PW DFT using the SIRIUS library", &
                                        "Molecular Mechanics", &
                                        "Hybrid quantum classical", &
                                        "Empirical Interatomic Potential", &
                                        "Electronic structure methods (DFT, ...)", &
                                        "Neural Network Potentials", &
                                        "Use a combination of two of the above", &
                                        "Perform an embedded calculation", &
                                        "Recieve forces from i–PI client"), &
                          enum_i_vals=(/do_qs, do_sirius, do_fist, do_qmmm, do_eip, do_qs, do_nnp, do_mixed, do_embed, do_ipi/), &
                          default_i_val=do_qs)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRESS_TENSOR", &
                          description="Controls the calculation of the stress tensor. The combinations defined below"// &
                          " are not implemented for all methods.", &
                          usage="stress_tensor (NONE|ANALYTICAL|NUMERICAL|DIAGONAL_ANA|DIAGONAL_NUM)", &
                          default_i_val=do_stress_none, &
                          enum_c_vals=s2a("NONE", "ANALYTICAL", "NUMERICAL", "DIAGONAL_ANALYTICAL", "DIAGONAL_NUMERICAL"), &
                          enum_i_vals=(/do_stress_none, do_stress_analytical, do_stress_numerical, &
                                        do_stress_diagonal_anal, do_stress_diagonal_numer/), &
                          enum_desc=s2a("Do not compute stress tensor", &
                                        "Compute the stress tensor analytically (if available).", &
                                        "Compute the stress tensor numerically.", &
                                        "Compute the diagonal part only of the stress tensor analytically (if available).", &
                                        "Compute the diagonal part only of the stress tensor numerically"))

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_ext_pot_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_rescale_force_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mix_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_embed_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_dft_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_pwdft_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mm_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_nnp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_qmmm_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_eip_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_bsse_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_subsys_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_properties_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_f_env_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_force_eval_section

! **************************************************************************************************
!> \brief Creates the section for applying an external potential
!> \param section ...
!> \date 03.2008
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_ext_pot_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXTERNAL_POTENTIAL", &
                          description="Section controlling the presence of an external potential dependent "// &
                          "on the atomic positions (X,Y,Z)", &
                          n_keywords=7, n_subsections=0, repeats=.TRUE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOMS_LIST", &
                          description="Specifies the atoms on which the external potential will act", &
                          usage="ATOMS_LIST {INT} {INT} ..", repeats=.TRUE., &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FUNCTION", &
                          description="Specifies the functional form in mathematical notation. Variables must be the atomic "// &
                          "coordinates (X,Y,Z).", usage="FUNCTION  X^2+Y^2+Z^2+LOG(ABS(X+Y))", &
                          type_of_var=lchar_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAMETERS", &
                          description="Defines the parameters of the functional form", &
                          usage="PARAMETERS a b D", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VALUES", &
                          description="Defines the values of  parameter of the functional form", &
                          usage="VALUES ", type_of_var=real_t, &
                          n_var=-1, repeats=.TRUE., unit_str="internal_cp2k")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Optionally, allows to define valid CP2K unit strings for each parameter value. "// &
                          "It is assumed that the corresponding parameter value is specified in this unit.", &
                          usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DX", &
                          description="Parameter used for computing the derivative with the Ridders' method.", &
                          usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERROR_LIMIT", &
                          description="Checks that the error in computing the derivative is not larger than "// &
                          "the value set; in case error is larger a warning message is printed.", &
                          usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ext_pot_section

! **************************************************************************************************
!> \brief Creates the section controlling the rescaling of forces
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_rescale_force_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RESCALE_FORCES", &
                          description="Section controlling the rescaling of forces. Useful when"// &
                          " starting from quite bad geometries with unphysically large forces.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_FORCE", &
                          description="Specify the Maximum Values of the force. If the force"// &
                          " of one atom exceed this value it's rescaled to the MAX_FORCE"// &
                          " value.", &
                          default_r_val=cp_unit_to_cp2k(value=50.0_dp, &
                                                        unit_str="kcalmol*angstrom^-1"), &
                          unit_str="hartree*bohr^-1")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_rescale_force_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE create_f_env_print_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      NULLIFY (keyword)
      NULLIFY (print_key)

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, &
                          name="PRINT", &
                          description="Properties that you want to output and that are common to all methods", &
                          n_keywords=0, n_subsections=10, repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="PROGRAM_RUN_INFO", &
                                       description="Controls the printing of basic information generated by FORCE_EVAL", &
                                       print_level=low_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, &
                          name="ENERGY_UNIT", &
                          description="Specifies the physical unit used for the printing of the total energy. "// &
                          "Note that the meaningfulness of the unit is not checked.", &
                          usage="ENERGY_UNIT eV", &
                          default_c_val="hartree", &
                          repeats=.FALSE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="FORCES", &
                                       description="Controls the printing of the forces after each force evaluation", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, &
                          name="NDIGITS", &
                          description="Specifies the number of digits used "// &
                          "for the printing of the forces", &
                          usage="NDIGITS 6", &
                          default_i_val=8, &
                          repeats=.FALSE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="FORCE_UNIT", &
                          variants=(/"UNIT"/), & ! add old keyword name for backward compatibility
                          description="Specifies the physical unit used for the printing of the forces. "// &
                          "Note that the meaningfulness of the unit is not checked.", &
                          usage="FORCE_UNIT eV/angstrom", &
                          default_c_val="hartree/bohr", &
                          repeats=.FALSE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="GRID_INFORMATION", &
                                      description="Controls the printing of information regarding the PW and RS grid structures.", &
                                       print_level=medium_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="TOTAL_NUMBERS", &
                                       description="Controls the printing of the total number of atoms, kinds, ...", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="DISTRIBUTION", &
                                       description="Controls the printing of the distribution of molecules, atoms, ...", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="DISTRIBUTION2D", &
                                       description="Controls the printing of the distribution of matrix blocks, ...", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="DISTRIBUTION1D", &
                                       description="Each node prints out its distribution info ...", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="STRESS_TENSOR", &
                                       description="Controls the printing of the stress tensor", &
                                       print_level=high_print_level, filename="__STD_OUT__")
      CALL keyword_create(keyword, __LOCATION__, &
                          name="COMPONENTS", &
                          description="Print all GPW/GAPW components contributing to the stress tensor", &
                          usage="COMPONENTS", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, &
                          name="STRESS_UNIT", &
                          description="Specifies the physical unit used for the printing of the stress tensor. "// &
                          "Note that the meaningfulness of the unit is not checked.", &
                          usage="STRESS_UNIT kbar", &
                          default_c_val="bar", &
                          repeats=.FALSE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="GRRM", &
                                       description="Controls the printing of the GRRM interface file", &
                                       print_level=debug_print_level + 1, filename="CP2K_GRRM")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, &
                                       name="SCINE", &
                                       description="Controls the printing of the SCINE interface file", &
                                       print_level=debug_print_level + 1, filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_f_env_print_section

END MODULE input_cp2k_force_eval
