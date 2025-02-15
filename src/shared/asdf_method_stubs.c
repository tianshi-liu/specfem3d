/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

// for ASDF reader setup

void FC_FUNC_(asdf_setup,ASDF_SETUP)(void) {
  fprintf(stderr,"ERROR: ASDF_FORMAT enabled without ASDF Support. To enable support, reconfigure with --with-asdf flag.\n");
  exit(1);
}
void FC_FUNC_(asdf_cleanup,ASDF_CLEANUP)(void) {}

// for ASDF writer

void FC_FUNC_(write_output_asdf,WRITE_OUTPUT_ASDF)(void) {
  fprintf(stderr,"ERROR: ASDF_FORMAT enabled without ASDF Support. To enable support, reconfigure with --with-asdf flag.\n");
  exit(1);
}

void FC_FUNC_(init_asdf_data,INIT_ASDF_DATA)(void) {}
void FC_FUNC_(store_asdf_data,STORE_ASDF_DATA)(void) {}
void FC_FUNC_(close_asdf_data,CLOSE_ASDF_DATA)(void) {}
void FC_FUNC_(write_asdf,WRITE_ASDF)(void) {}

// for ASDF reader

void FC_FUNC_(read_adjoint_sources_asdf,READ_ADJOINT_SOURCES_ASDF)(void) {}
void FC_FUNC_(check_adjoint_sources_asdf,CHECK_ADJOINT_SOURCES_ASDF)(void) {}
