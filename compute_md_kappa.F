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
      subroutine compute_md_kappa(temp, rho, kappa)
#ifdef SGCHEM
      implicit none
#include "cool.h"
#include "kappa-planck.h"
      REAL temp, kappa
      REAL rho, logT, dT, logrho, rho_offset
      REAL drho, k1, k2, k3, k4, ka, kb
      integer iT, irho, idx, idx1, idx2, idx3, idx4

      logT = log10(temp)
      iT = int((logT - 1.8) / 0.1) + 1
      if (iT .lt. 1) then
c At low T, use the opacities for log(T) = 1.8, in absence of anything
c better. NB 10**1.8 ~ 63K, so we probably won't have gas falling off the
c bottom of the table in the pop III case
        iT = 1
        dT = 0d0
      elseif (iT .ge. 29) then
c At high T, use upper limit in table
        iT = 29
        dT = 0d0
      else
        dT = (logT - (0.1*(iT-1) + 1.8)) / 0.1
      endif

      logrho = log10(rho)
      irho = int(logrho + 16.00) + 1
      rho_offset = 0d0

      if (irho .lt. 1) then
c At v. low density, assume constant opacity for log(T) < 3.1, 
c opacity ~ density at higher T. This isn't quite right, but should
c be accurate to within a factor of a few
        irho    = 1
        drho    = 0d0
        if (iT .gt. 14) then
          rho_offset = logrho + 16
        endif
      elseif (irho .ge. 15) then
c Keep constant if log(T) > 3.8, otherwise scale with density
        irho = 15
        drho = 0d0
        if (iT .le. 21) then
          rho_offset = logrho + 2
        endif
      else
         drho = (logrho - ((irho-1) - 16.))
      endif
c 
c Now do lookup
c
      if (dT .eq. 0d0) then
        if (drho .eq. 0d0) then
          idx = 15 * (iT - 1) + irho            
          kappa = kappa_p(idx) + rho_offset
        else
          idx1 = 15 * (iT - 1) + irho
          idx2 = idx1 + 1
          k1 = kappa_p(idx1)
          k2 = kappa_p(idx2)
          kappa = k1 + drho * (k2 - k1)
        endif
      else
        if (drho .eq. 0d0) then
          idx1 = 15 * (iT - 1) + irho            
          idx2 = 15 * iT + irho
          k1 = kappa_p(idx1)
          k2 = kappa_p(idx2)
          kappa = k1 + dT * (k2 - k1) + rho_offset
        else
          idx1 = 15 * (iT - 1) + irho
          idx2 = idx1 + 1
          k1 = kappa_p(idx1)
          k2 = kappa_p(idx2)
          ka = k1 + drho * (k2 - k1)
c
          idx3 = idx1 + 15
          idx4 = idx3 + 1
          k3 = kappa_p(idx3)
          k4 = kappa_p(idx4)
          kb = k3 + drho * (k4 - k3)
c
          kappa = ka + dT * (kb - ka)
        endif
      endif

c Finally, convert from log -> actual value
      kappa = 1d1**kappa

#endif /* SGCHEM */
      return
      end

