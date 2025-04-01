c
c     Copyright 2025 Simon Glover      
c
c     This file is part of SGPDR
c      
c     SGPDR is free software: you can redistribute it and/or modify it under the terms
c     of the GNU General Public License as published by the Free Software Foundation, either
c     version 3 of the License, or (at your option) any later version.
c
c     This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
c     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c     See the GNU General Public License for more details.
c      
c     You should have received a copy of the GNU General Public License along with this program.
c     If not, see <https://www.gnu.org/licenses/>.
c
#include "config.h"
#ifdef SGCHEM
c
c     Dummy routine -- in the mode in which it is used in SGPDR, the DVODE solver requires a routine
c     called 'jac' to be present in order to compile and link properly, but does not use it to compute
c     the Jacobian of the set of rate equations, which is instead computed internally using finite
c     differencing    
c
      SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
      IMPLICIT NONE
#include "cool.h"
      integer NEQ, IPAR
      REAL T, Y(NEQ), RPAR
      integer ML, MU, NROWPD
      REAL PD(NROWPD,NEQ)
      return
      end
#endif
