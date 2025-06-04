!  ┏┓┏┓┏┓  Licensed under the MIT License
!  ┃ ┗┓┣   Copyright (c) 2025 Rodrigo Castro 
!  ┗┛┗┛┻   https://github.com/rodpcastro/colspecf

module csf_kinds
!! ColSpecF kinds for numeric types.

  use, intrinsic :: iso_fortran_env, only: int8, int16, int32, real32, real64

  implicit none
  private
  public :: i1, i2, i4, sp, dp, wp

  integer(int8), parameter :: i1 = int8   !! 1-byte integer
  integer(int8), parameter :: i2 = int16  !! 2-byte integer
  integer(int8), parameter :: i4 = int32  !! 4-byte integer

  integer(int8), parameter :: sp = real32  !! single precision real
  integer(int8), parameter :: dp = real64  !! double precision real
  integer(int8), parameter :: wp = real64  !! working precision real

end module csf_kinds
