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
      program SGPDR
#ifdef SGCHEM
#define NPIX  1
#include "sgchem_def.h"
      implicit NONE
      real*8 abundD, abundc, abundo, abundM, Zatom, dust_temp, G0,
     $       LWBGz, D2G, crion, initz, h2ex, h2kin, AVeff, time,
     $       dl, yn, divv, energy, current_redshift,
     $       temperature_init
      real*8 non_eq_abundances(SGCHEM_NUM_SPECIES)
      real*8 non_eq_abundances_init(SGCHEM_NUM_SPECIES)
      real*8 photorates_direct(NDIRECT)
      real*8 fshield_H2, fshield_CO, fshield_C, AV_mean, ntot
      real*8 flux_mean_temp
      real*8 L, NH2, NCO, NCI, temp, dNH2
      real*8 partition, fu, flux, dflux
      integer LWBG, iphoto, isrf, iatom, I, J, iflag_h2_opacity

c Main adjustable parameters
      open(10, file='SGPDR_params.txt', status='old')
c Elemental abundances of D, C, O and "M" relative to H, by number
c "M" is a representative low ionisation potential metal (e.g. silicon)
      read(10,*) abundD    
      read(10,*) abundc
      read(10,*) abundo
      read(10,*) abundM
c Gas-phase metallicity (in units of solar)
      read(10,*) Zatom
c Dust:gas ratio (relative to local ISM)    
      read(10,*) D2G
c Strength of interstellar radiation field in Habing units
      read(10,*) G0
c Cosmic ray ionisation rate of **atomic** hydrogen [value for H2 is ~ 2x larger]
      read(10,*) crion
c Atomic cooling option -- set to 3 unless you know what you're doing
      read(10,*) iatom      
      close(10, status = "keep")

c A bunch of other parameters that we need to set up for the run. I don't advise
c varying these unless you know what they do
      current_redshift = 0.0
      dust_temp = 2d0
      dl     = 0.0
      LWBG   = 0
      LWBGz  = 0d0
      initz  = 0d0
      h2ex   = 0.84d0
      h2kin  = 0.12d0
      AVeff  = 0d0
      iphoto = 6
      isrf   = 1
      iflag_h2_opacity = 0
      
      call init_chemistry_parameters(abundD,
#ifndef SGCHEM_VARIABLE_Z
     $ abundc, abundo, abundM, Zatom,
#endif
     $ dust_temp, G0, LWBG, LWBGz,
#ifndef SGCHEM_VARIABLE_Z
     $ D2G,
#endif
#ifndef SGCHEM_VARIABLE_CRION
     $ crion,
#endif
     $ initz, AVeff, h2ex, h2kin, iphoto, isrf, iatom, iflag_h2_opacity)

      call coolinmo
      call cheminmo
      call init_tolerances

#ifdef SGCHEM_VARIABLE_Z
      call set_local_dust_abundance(D2G)
      call set_local_element_abundances(abundc, abundo, abundM,
     $                                  Zatom)
#endif

c Ionization/dissociation rates for point source contribution
      do I = 1, 8
        photorates_direct(I) = 0d0
      enddo
      
c     Initial conditions for chemical update
      open(11, file='SGPDR_IC.txt', status='old')
c Size of chemistry timestep
      read(11,*) time
c Velocity divergence -- sensible default value is 1 km/s/pc ~ 3.2e-14
c Note: only influences [CII] and CO cooling in optically thick limit
      read(11,*) divv
c Fractional abundances of the chemical species for which we're solving ODEs
c For our modified version of the Gong et al. (2017) chemistry network, these
c are: 1 -- H2, 2 -- H+, 3 -- C+, 4 -- CHx, 5 -- OHx, 6 -- CO, 7 -- HCO+, 8 -- He+, 9 -- M+.     
      read(11,*) non_eq_abundances_init
c Initial temperature
      read(11,*) temperature_init
c     Slab depth (in pc)
      read(11,*) L
      close(11, status = "keep")

c Convert slab depth to cgs, compute width of each slice      
      L  = L * 3.086e18
      dl = L / 3d3

c Open files for output
      open(14, file='SGPDR_output.txt', status='replace')
#ifdef OUTPUT_MEAN_TEMPS
      open(15, file='SGPDR_temps.txt', status='replace')
#endif
      
      do J = 1, 31
         yn = 1d1**(2.0 + 0.1 * (J-1))
         print*, "Solving slab for density n = ", yn
         do I = 1, SGCHEM_NUM_SPECIES
            non_eq_abundances(I) = non_eq_abundances_init(I)
         enddo
         
         energy = 2.07e-13 * temperature_init * (yn / 1d3) 
     $        * (1.1 + non_eq_abundances(IHP+1)
     $        + non_eq_abundances(IHEP+1)
     $        - non_eq_abundances(IH2+1))

         NH2 = 0d0
         NCO = 0d0
         NCI = 0d0
         flux = 0.0
         flux_mean_temp = 0.0
      
         do I = 1, 3000
            AV_mean = 5.348d-22 * yn * dl * I * D2G
            ntot = yn * (1.1 + non_eq_abundances(IHP+1)
     $                       + non_eq_abundances(IHEP+1)
     $                       - non_eq_abundances(IH2+1))
            temp = energy / (1.5 * 1.38d-16 * ntot)
         
            call calc_H2_shielding_factor(temp, 0.0, NH2, fshield_H2) 
            call calc_CO_shielding_factor(NH2, NCO, 0.0, 0.0,
     $                                    fshield_CO)
            call calc_C_shielding_factor(NH2, NCI, 0.0, 0.0, fshield_C)

          
            call evolve_abundances(time, dl, yn, divv, energy,
     $                             current_redshift, non_eq_abundances,
     $                             photorates_direct, dust_temp,
     $                             fshield_H2, fshield_CO, fshield_C,
     $                             AV_mean)

            dNH2 = non_eq_abundances(IH2+1) * yn * dl
            NH2 = NH2 + dNH2
            NCO = NCO + non_eq_abundances(ICO+1) * yn * dl
            NCI = NCI + (abundC - non_eq_abundances(ICP+1)
     $                - non_eq_abundances(ICO+1)) * yn * dl         
            ntot = yn * (1.1 + non_eq_abundances(IHP+1)
     $                    - non_eq_abundances(IH2+1))
            temp = energy / (1.5 * 1.38d-16 * ntot)

            write(14, *) yn, I, non_eq_abundances, energy, temp

#ifdef OUTPUT_MEAN_TEMPS            
c Compute S(1) emission from this slice of slab
            partition = 0.024 * temp / (1 - exp(-6d3/temp))
            fu = 21 * exp(-1015/temp) / partition
c Note: A_31 * E_31 = 5.53e-23 erg/s
            dflux = fu * dNH2 * 5.53e-23 / (4 * 3.1415927) 
            flux = flux + dflux
c Flux-weighted mean temperature -- need to normalise at end by total flux
            flux_mean_temp = flux_mean_temp + dflux * temp
#endif
         enddo
#ifdef OUTPUT_MEAN_TEMPS         
       write(15, *) yn, flux_mean_temp/flux
#endif
      enddo

      close(14, status='keep')
#ifdef OUTPUT_MEAN_TEMPS
      close(15, status='keep')
#endif

      
#endif /* SGCHEM */
      end
