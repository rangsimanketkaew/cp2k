!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A simple hash table of integer keys, using hash function:
!>          H(k) = (k*p) mod n + 1
!>        where:
!>          k = key
!>          p = a prime number >= n
!>          n = size of the hash table
!>         And collision resolvation is done by open addressing with linear
!>         probing.
!>
!>         The table consists of an array of (key,val) pairs, and
!>         there are no intermediate buckets. For every new entry (k,v):
!>         We first look up slot H(k), and if it already contains an entry,
!>         then move to the next empty slot using a predefined linear probing
!>         sequence (e.g. iterate from slots H(k) to n, and then 1 to H(k)-1).
!>         When we look up, we use the same probing sequence.
! **************************************************************************************************
MODULE qs_fb_hash_table_types

   USE kinds,                           ONLY: int_8
   USE qs_hash_table_functions,         ONLY: hash_table_matching_prime
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! public types
   PUBLIC :: fb_hash_table_obj

! public methods
!API
   PUBLIC :: fb_hash_table_add, &
             fb_hash_table_create, &
             fb_hash_table_get, &
             fb_hash_table_has_data, &
             fb_hash_table_nullify, &
             fb_hash_table_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_hash_table_types'

! key value indicating an empty slot
   INTEGER(KIND=int_8), PARAMETER, PRIVATE :: EMPTY_KEY = -1_int_8
! Parameters related to automatic resizing of the hash_table:
! Resize by EXPAND_FACTOR if total no. slots / no. of filled slots < ENLARGE_RATIO
   INTEGER, PARAMETER, PRIVATE :: ENLARGE_RATIO = 1
   INTEGER, PARAMETER, PRIVATE :: REDUCE_RATIO = 3
   INTEGER, PARAMETER, PRIVATE :: EXPAND_FACTOR = 2
   INTEGER, PARAMETER, PRIVATE :: SHRINK_FACTOR = 2

! **************************************************************************************************
!> \brief hash table entry data type
!> \param key       : key of the entry
!> \param val       : value of the entry
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   TYPE fb_hash_table_element
      INTEGER(KIND=int_8) :: key = -1_int_8
      INTEGER :: val = -1
   END TYPE fb_hash_table_element

! **************************************************************************************************
!> \brief data defining a hash table using open addressing for collision
!>        resolvation. Uses simple entry structure to be memory efficient
!>        as well as small overhead
!> \param table     : hash table data area
!> \param nelements : number of non-empty slots in table
!> \param nmax      : max number of slots in table
!> \param prime     : prime number used in the hash function
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   TYPE fb_hash_table_data
      TYPE(fb_hash_table_element), DIMENSION(:), POINTER :: table => NULL()
      INTEGER :: nelements = -1
      INTEGER :: nmax = -1
      INTEGER :: prime = -1
   END TYPE fb_hash_table_data

! **************************************************************************************************
!> \brief the object container which allows for the creation of an array
!>        of pointers to fb_hash_table objects
!> \param obj : pointer to the fb_hash_table object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   TYPE fb_hash_table_obj
      TYPE(fb_hash_table_data), POINTER, PRIVATE :: obj => NULL()
   END TYPE fb_hash_table_obj

CONTAINS

! **************************************************************************************************
!> \brief Add element to a hash table, auto resize if necessary
!> \param hash_table : the fb_hash_table object
!> \param key        : key of the element
!> \param val        : value of the element
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   RECURSIVE SUBROUTINE fb_hash_table_add(hash_table, key, val)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key
      INTEGER, INTENT(IN)                                :: val

      INTEGER                                            :: islot
      LOGICAL                                            :: check_ok

      check_ok = fb_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      ! check hash table size, if too small rehash in a larger table
      IF (hash_table%obj%nelements*ENLARGE_RATIO .GE. &
          hash_table%obj%nmax) THEN
         CALL fb_hash_table_rehash(hash_table=hash_table, &
                                   nmax=hash_table%obj%nmax*EXPAND_FACTOR)
      END IF
      ! find the right slot for the given key
      islot = fb_hash_table_linear_probe(hash_table, key)
      CPASSERT(islot > 0)
      ! we are adding a new entry only if islot points to an empty slot,
      ! otherwise just change the val of the existing entry
      IF (hash_table%obj%table(islot)%key == EMPTY_KEY) THEN
         hash_table%obj%nelements = hash_table%obj%nelements + 1
         hash_table%obj%table(islot)%key = key
      END IF
      hash_table%obj%table(islot)%val = val
   END SUBROUTINE fb_hash_table_add

! **************************************************************************************************
!> \brief Creates and initialises an empty fb_hash_table object
!> \param hash_table : the fb_hash_table object, its content must be NULL
!>                     and cannot be UNDEFINED
!> \param nmax       : total size of the table, optional. If absent default
!>                     size is 1.
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE fb_hash_table_create(hash_table, nmax)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(IN), OPTIONAL                      :: nmax

      INTEGER                                            :: my_nmax
      LOGICAL                                            :: check_ok

      check_ok = .NOT. fb_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      ALLOCATE (hash_table%obj)
      NULLIFY (hash_table%obj%table)
      hash_table%obj%nmax = 0
      hash_table%obj%nelements = 0
      hash_table%obj%prime = 2
      my_nmax = 1
      IF (PRESENT(nmax)) my_nmax = nmax
      CALL fb_hash_table_init(hash_table=hash_table, &
                              nmax=my_nmax)

   END SUBROUTINE fb_hash_table_create

! **************************************************************************************************
!> \brief Retrieve value from a key from a hash table
!> \param hash_table : the fb_hash_table object
!> \param key        : input key
!> \param val        : output value, equals to 0 if key not found
!> \param found      : .TRUE. if key is found, .FALSE. otherwise
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE fb_hash_table_get(hash_table, key, val, found)
      TYPE(fb_hash_table_obj), INTENT(IN)                :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key
      INTEGER, INTENT(OUT)                               :: val
      LOGICAL, INTENT(OUT)                               :: found

      INTEGER                                            :: islot
      LOGICAL                                            :: check_ok

      check_ok = fb_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      found = .FALSE.
      val = 0
      islot = fb_hash_table_linear_probe(hash_table, key)
      IF (islot > 0) THEN
         IF (hash_table%obj%table(islot)%key == key) THEN
            val = hash_table%obj%table(islot)%val
            found = .TRUE.
         END IF
      END IF
   END SUBROUTINE fb_hash_table_get

! **************************************************************************************************
!> \brief check if the object has data associated to it
!> \param hash_table : the fb_hash_table object in question
!> \return : true if hash_table%obj is associated, false otherwise
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   PURE FUNCTION fb_hash_table_has_data(hash_table) RESULT(res)
      TYPE(fb_hash_table_obj), INTENT(IN)                :: hash_table
      LOGICAL                                            :: res

      res = ASSOCIATED(hash_table%obj)
   END FUNCTION fb_hash_table_has_data

! **************************************************************************************************
!> \brief Initialises a fb_hash_table object
!> \param hash_table : the fb_hash_table object, its content must be NULL
!>                     and cannot be UNDEFINED
!> \param nmax       : new size of the table, optional. If absent use the
!>                     old size
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE fb_hash_table_init(hash_table, nmax)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(IN), OPTIONAL                      :: nmax

      INTEGER                                            :: ii, my_nmax, power
      LOGICAL                                            :: check_ok

      check_ok = fb_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      my_nmax = hash_table%obj%nmax
      IF (PRESENT(nmax)) my_nmax = nmax
      ! table length should always be power of 2. Find the least
      ! power that is greater or equal to my_nmax
      power = 0
      DO WHILE (2**power .LT. my_nmax)
         power = power + 1
      END DO
      my_nmax = 2**power
      IF (ASSOCIATED(hash_table%obj%table)) THEN
         IF (SIZE(hash_table%obj%table) .NE. my_nmax) THEN
            DEALLOCATE (hash_table%obj%table)
            ALLOCATE (hash_table%obj%table(my_nmax))
         END IF
      ELSE
         ALLOCATE (hash_table%obj%table(my_nmax))
      END IF
      hash_table%obj%nmax = my_nmax
      hash_table%obj%prime = hash_table_matching_prime(my_nmax)
      ! initiate element to be "empty"
      DO ii = 1, hash_table%obj%nmax
         hash_table%obj%table(ii)%key = EMPTY_KEY
         hash_table%obj%table(ii)%val = 0
      END DO
      hash_table%obj%nelements = 0
   END SUBROUTINE fb_hash_table_init

! **************************************************************************************************
!> \brief Nullifies a fb_hash_table object
!> \param hash_table : the fb_hash_table object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   PURE SUBROUTINE fb_hash_table_nullify(hash_table)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table

      NULLIFY (hash_table%obj)
   END SUBROUTINE fb_hash_table_nullify

! **************************************************************************************************
!> \brief Rehash table. If nmax is present, then also change the table size
!>        to MAX(nmax, number_of_non_empty_elements).
!> \param hash_table      : the fb_hash_table object
!> \param nmax [OPTIONAL] : maximum size of the rehashed table
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   RECURSIVE SUBROUTINE fb_hash_table_rehash(hash_table, nmax)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(IN), OPTIONAL                      :: nmax

      INTEGER                                            :: ii, my_nmax
      TYPE(fb_hash_table_element), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_table

      IF (.NOT. fb_hash_table_has_data(hash_table)) THEN
         CALL fb_hash_table_create(hash_table, nmax)
         RETURN
      END IF
      IF (PRESENT(nmax)) THEN
         my_nmax = MAX(nmax, hash_table%obj%nelements)
      ELSE
         my_nmax = hash_table%obj%nmax
      END IF
      ALLOCATE (tmp_table(hash_table%obj%nmax))
      tmp_table(:) = hash_table%obj%table(:)
      CALL fb_hash_table_release(hash_table)
      CALL fb_hash_table_create(hash_table=hash_table, &
                                nmax=my_nmax)
      DO ii = 1, SIZE(tmp_table)
         IF (tmp_table(ii)%key .NE. EMPTY_KEY) THEN
            CALL fb_hash_table_add(hash_table=hash_table, &
                                   key=tmp_table(ii)%key, &
                                   val=tmp_table(ii)%val)
         END IF
      END DO
      DEALLOCATE (tmp_table)
   END SUBROUTINE fb_hash_table_rehash

! **************************************************************************************************
!> \brief releases given object
!> \param hash_table : the fb_hash_table object in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE fb_hash_table_release(hash_table)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table

      IF (ASSOCIATED(hash_table%obj)) THEN
         IF (ASSOCIATED(hash_table%obj%table)) THEN
            DEALLOCATE (hash_table%obj%table)
         END IF
         DEALLOCATE (hash_table%obj)
      ELSE
         NULLIFY (hash_table%obj)
      END IF
   END SUBROUTINE fb_hash_table_release

! **************************************************************************************************
!> \brief Remove element from a table, automatic resize if necessary
!> \param hash_table : the fb_hash_table object
!> \param key        : key of the element to be removed
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE fb_hash_table_remove(hash_table, key)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key

      INTEGER                                            :: islot
      LOGICAL                                            :: check_ok

      check_ok = fb_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      islot = fb_hash_table_linear_probe(hash_table, key)
      ! we are only removing an entry if the key is found
      IF (islot > 0) THEN
         IF (hash_table%obj%table(islot)%key == key) THEN
            hash_table%obj%table(islot)%key = EMPTY_KEY
            hash_table%obj%nelements = hash_table%obj%nelements - 1
            ! must rehash after setting a filled slot to empty, otherwise the
            ! table will not work. Automatic resize if required
            IF (hash_table%obj%nelements*REDUCE_RATIO .LT. &
                hash_table%obj%nmax) THEN
               CALL fb_hash_table_rehash(hash_table=hash_table, &
                                         nmax=hash_table%obj%nmax/SHRINK_FACTOR)
            ELSE
               CALL fb_hash_table_rehash(hash_table=hash_table)
            END IF
         END IF
      END IF
   END SUBROUTINE fb_hash_table_remove

! **************************************************************************************************
!> \brief outputs the current information about the table
!> \param hash_table : the fb_hash_table object in question
!> \param nelements  : number of non-empty slots in the table
!> \param nmax       : maximum number of slots in the table
!> \param prime      : the prime used in the hash function
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE fb_hash_table_status(hash_table, nelements, nmax, prime)
      TYPE(fb_hash_table_obj), INTENT(INOUT)             :: hash_table
      INTEGER, INTENT(OUT), OPTIONAL                     :: nelements, nmax, prime

      LOGICAL                                            :: check_ok

      check_ok = fb_hash_table_has_data(hash_table)
      CPASSERT(check_ok)
      IF (PRESENT(nelements)) nelements = hash_table%obj%nelements
      IF (PRESENT(nmax)) nmax = hash_table%obj%nmax
      IF (PRESENT(prime)) prime = hash_table%obj%prime
   END SUBROUTINE fb_hash_table_status

! **************************************************************************************************
!> \brief Linear probing algorithm for the hash table
!> \param hash_table : the fb_hash_table object
!> \param key        : key to locate
!> \return : slot location in the table correspond to key, 0 if key
!>                     not found
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   PURE FUNCTION fb_hash_table_linear_probe(hash_table, key) &
      RESULT(islot)
      TYPE(fb_hash_table_obj), INTENT(IN)                :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key
      INTEGER                                            :: islot

      INTEGER                                            :: guess

! first guess is mapped by the hash_function

      guess = fb_hash_table_hash_function(hash_table, key)
      ! then search for key and stop at first empty slot from guess to
      ! nmax.  using the same linear probe for adding and retrieving
      ! makes all non-empty keys being put before the first empty slot.
      DO islot = guess, hash_table%obj%nmax
         IF ((hash_table%obj%table(islot)%key == key) .OR. &
             (hash_table%obj%table(islot)%key == EMPTY_KEY)) RETURN
      END DO
      ! if unsuccessful, search from 1 to guess
      DO islot = 1, guess - 1
         IF ((hash_table%obj%table(islot)%key == key) .OR. &
             (hash_table%obj%table(islot)%key == EMPTY_KEY)) RETURN
      END DO
      ! if not found and table is full set islot to 0
      islot = 0
   END FUNCTION fb_hash_table_linear_probe

! **************************************************************************************************
!> \brief Hash function
!> \param hash_table : the fb_hash_table object
!> \param key        : key to locate
!> \return : slot location in the table correspond to key, 0 if key
!>                     not found
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   PURE FUNCTION fb_hash_table_hash_function(hash_table, key) RESULT(hash)
      TYPE(fb_hash_table_obj), INTENT(IN)                :: hash_table
      INTEGER(KIND=int_8), INTENT(IN)                    :: key
      INTEGER                                            :: hash

      INTEGER(KIND=int_8)                                :: hash_8, nmax_8, prime_8

      nmax_8 = INT(hash_table%obj%nmax, int_8)
      prime_8 = INT(hash_table%obj%prime, int_8)
      ! IAND with nmax-1 is equivalent to MOD nmax if nmax is alway a power of 2.
      hash_8 = IAND(key*prime_8, nmax_8 - 1) + 1_int_8
      hash = INT(hash_8)
   END FUNCTION fb_hash_table_hash_function

END MODULE qs_fb_hash_table_types
