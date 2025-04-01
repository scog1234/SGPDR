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
      subroutine compute_heating(flux, temp, rho, heating_rate)
#ifdef SGCHEM
      implicit none
#include "cool.h"
      REAL flux, temp, temp_local, rho, heating_rate
      REAL kappa, specific_heating_rate
c
c Cap the temperature used for the heating rate calculation.
c This is a kludge necessitated by the fact that we currently
c assume that T_rad = T_gas when computing the heating rate.
c This is reasonable when T_gas is a few thousand K, but gives
c us far too hard a radiation field for T_gas > 1e4 K, and hence
c an extremely large heating rate.
c
c Once I rework the opacity calculation to remove the T_gas = T_rad
c assumption, the need for this hack should go away.
c
      temp_local = temp
      if (temp_local .gt. 6d3) then
        temp_local = 6d3
      endif
c
c Compute planck-mean kappa, given T, n   [cm^2 g^-1]
c
      call compute_md_kappa(temp_local, rho, kappa)
c
c Compute specific heating rate
c
      specific_heating_rate = flux * kappa
c
c Convert to volumetric heating rate
c
      heating_rate = rho * specific_heating_rate
c
#endif /* SGCHEM */
      return
      end

