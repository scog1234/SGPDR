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
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////     E V O L V E _ A B U N D A N C E S     \\\\\\\\\\
c
c=======================================================================
c
      subroutine evolve_abundances(time, dl, yn, divv, energy, 
     $                             current_redshift, non_eq_abundances,
     $                             photorates_direct, dust_temp,
     $                             fshield_H2, fshield_CO, fshield_C,
     $                             AV_mean)
#ifdef SGCHEM
      implicit NONE
#include "cool.h"
#include "non_eq.h"
c
      REAL    time, dl, yn, divv, energy, current_redshift
      REAL    non_eq_abundances(num_non_eq_species)
      REAL    photorates_direct(NDIRECT)
      REAL    dust_temp, fshield_H2, fshield_CO, fshield_C
      
      REAL    abe, abh2, abhep, abhepp, abco, abHI, abcI
      REAL    yntot, ekn, temp, gamma
      REAL    t_start, t_end, dt_est, dt, t_in, t_out
      REAL    tdust_old
      REAL    fshield_dust, AV_mean,
     $        chi_mean, fshield_dust_H2, fshield_dust_CO
      REAL    RH2, Rdust, uvdust
      integer new_step,ntry
c
      integer nstep
      parameter (nstep = 5000000)
c
      REAL safety(nspec)
      integer used_for_stepsize_control(nspec)
#if CHEMISTRYNETWORK == 1
      DATA safety /1d-1, 0d0, 0d0, 0d0, 0d0, 0d0,
     $             DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 0, 0, 0, 0, 0, 1/
#endif
#if CHEMISTRYNETWORK == 4
      DATA safety /1d-1, 1d-1, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1/
#endif
#if CHEMISTRYNETWORK == 5
      DATA safety /1d-1, 1d-1, 1d-1, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1, 1/
#endif
#if CHEMISTRYNETWORK == 7
      DATA safety /1d-1, 1d-1, 1d-1, 1d-1, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1, 1, 1/
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      DATA safety /1d-1, 0d0, 1d-1, 0d0, 0d0, 1d-1, 0d0, 0d0, 
     $             0d0, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 0, 1, 0, 0, 1, 0, 0, 0, 1/
#endif
c
      integer I, J
c
c DVODE control variables -- NB atol declared in non_eq.h
c
      REAL    rtol
      integer itask, iopt, itol, mf, istate, iflag_dvi
c
c Other DVODE variables
c
      integer lrw, liw
      parameter (lrw   = 2000)
      parameter (liw   = 2000)
      REAL    rpar(nrpar), rwork(lrw)
      integer ipar(nipar), iwork(liw)
      REAL    y(nspec), ydot(nspec), y_old(nspec)
c
      REAL chi
      external jac, rate_eq, chi
c
c Sanity check for our input time
c
      if (time .le. 0d0) then
#ifdef DEBUG_EVOLVE
        if (time .lt. -1d0) then
          print*, "Warning: negative timestep in evolve_abundances",
     $             time
        endif
#endif
        return
      endif      
c
c Set redshift
c
      redshift = current_redshift
c
c Initialize abundance vector, treating the internal energy as 
c an additional abundance. 
c
c N.B. We assume that abundances are stored consecutively, with
c no gaps, and that itmp == nspec (i.e. that the energy is the last
c value stored)
c
      do I = 1, nspec-1
        y(I) = non_eq_abundances(I)
      enddo
      y(itmp)  = energy

      call check_advection(y)
c
c Setup initial parameters for DVODE integrator. These don't vary over 
c the course of the timestep.
c
c Optional input supplied -- see iwork, rwork init. below:
c
      iopt = 1
c
c Evolve ODEs from t -> tout, allowing overshooting & interpolation
c [See DVODE documentation for more details].
c
      itask = 1
c
c Use same relative tolerance for each species, different absolute tolerances
c (since magnitudes may vary by many orders of magnitude). 
c
      itol = 2
c
c Set relative tolerance -- for simplicity we keep this fixed at 0.1% for now  
c
      rtol = RTOL
c
c ODE integration method: 22 => BDF with an internally generated Jacobian
c
      mf = 22
c
c Start DVODE in start-up mode (which forces it to initialize all of its
c internal arrays). Subsequent calls within the loop below use istate = 2
c
      istate = 1
c
c Optional input -- increase maximum number of steps that DVODE can take 
c internally without returning an to 20000, as 500 doesn't appear to be 
c enough for portions of our parameter space
c
      do I = 5, 10      
        rwork(I) = 0d0
        iwork(I) = 0
      enddo
c      
      iwork(6) = 200000
c
c Compute shielding factors for photochemistry
c
      abH2 = y(ih2)
#if CHEMISTRYNETWORK == 1
      abhep  = y(ihep)
      abhepp = y(ihepp)
      abco   = 0d0
      abcI   = 0d0
#endif
#if CHEMISTRYNETWORK == 4
      abco   = 0d0
      abhep  = 0d0
      abhepp = 0d0
      abcI   = 0d0
#endif
#if CHEMISTRYNETWORK == 5
      abco   = y(ico)
      abhep  = 0d0
      abhepp = 0d0
      abcI   = 0d0
#endif
#if CHEMISTRYNETWORK == 7
      abco   = y(ico)
      abhep  = y(ihep)
      abhepp = 0d0
      abcI   = 0d0
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      abco   = y(ico)
      abhep  = y(ihep)
      abhepp = 0d0
      abcI   = abundc - y(ico) - y(icp) - y(ichx) - y(ihcop)
#endif
c
c Compute temperature
c
c Unless VARIABLE_GAMMA is defined, we assume that gamma = 5/3. This is a
c good approximation for gas with a small molecular fraction, or for
c fully molecular gas with a low temperature (as found in most Galactic GMCs).
c
      abe = y(ihp) + abhep + abhepp
      yntot = (1d0 + ABHE - abh2 + abe) * yn
#ifdef VARIABLE_GAMMA
      ekn = energy / (kboltz * yntot)
      call lookup_temp_gamma(abh2, ekn, temp, gamma)
#else
      gamma = GAMMA
      temp = (gamma - 1d0) * energy / (yntot * kboltz)
#endif

      if (temp.gt.1.e10) then
         print*, 'SGCHEM WARNING! Resetting the temperature',yn,energy
         temp = 1.e8
         energy = temp * (yntot * kboltz) / (gamma - 1d0) 
      endif
c     
#if CHEMISTRYNETWORK == 16
        fshield_dust    = exp(-1.87d0 * AV_mean)
        fshield_dust_H2 = exp(-4.18d0 * AV_mean)
        fshield_dust_CO = exp(-3.88d0 * AV_mean)
#else
        fshield_dust    = exp(-2.5d0 * AV_mean)
        fshield_dust_H2 = exp(-3.5d0 * AV_mean)
        fshield_dust_CO = 0d0
#endif
c
        chi_mean = chi(AV_mean)
c
c Parameters to pass to rate_eq
c
      rpar(1) = yn            ! H nuclei number density
      rpar(2) = 0.0           ! Cell size -- not used in this version
      rpar(3) = divv          ! Velocity divergence
      rpar(4) = fshield_H2    ! H2 self-shielding factor
      rpar(5) = fshield_CO    ! CO self-shielding factor
      rpar(6) = fshield_dust  ! Scaling used for CO, PE effect etc.
      rpar(7) = chi_mean      ! ISRF heating rate of dust, relative to value in unshielded gas
      rpar(8) = fshield_dust_H2     ! Dust shielding factor (H2; only for Gong et al network)
      rpar(9)  = fshield_dust_CO    ! Dust shielding factor (CO; only for Gong et al network)
      rpar(10) = AV_mean      ! Visual extinction
      rpar(11) = 0.0          ! Accretion luminosity from sink; not used in standalone version
      rpar(20) = fshield_C    ! CI self-shielding factor

c Photoionization / photodissociation from point sources
c Computed in RT module (e.g. SWEEP) before input here
c For standalone calculations, often want to set these to zero    
c
      rpar(12) = photorates_direct(1)  ! phih:   HI photoionization rate
      rpar(13) = photorates_direct(2)  ! hvphih: Heating rate (H photoionization)
      rpar(14) = photorates_direct(3)  ! h2phih: H2 photoionization
      rpar(15) = photorates_direct(4)  ! h2vphih: Heating rate (H2 photoionization)
      rpar(16) = photorates_direct(5)  ! uvphih: H2 photodissociation rate
      rpar(17) = photorates_direct(6)  ! uvdust: Dust heating rate
      rpar(18) = photorates_direct(7)  ! hephih:  HeI photoionization rate
      rpar(19) = photorates_direct(8)  ! hevphih: HeI photoionization heating rate
c
      ipar(1) = 0
c
      t_start  = 0d0
      t_end    = time
      new_step = 1
      ntry     = 1

      do I = 1, nstep
c
c If we're starting a new step, then figure out a reasonable timestep to
c pass to DVODE -- we need to do this here to avoid DVODE producing -ve
c abundances or -ve temperatures, which can occur if it starts with too
c large a timestep.
c
c If we're not starting a new step, then dt has already been set. 
c
        if (new_step .eq. 1) then
c
c First, save a copy of the current abundances
c
          do J = 1, nspec
            y_old(J) = y(J)
          enddo
          tdust_old = tdust

c Debug a particular particle 
#ifdef DEBUG_PARTICLE_ID
          if (id_current .eq. DEBUG_PARTICLE_ID) then
             print*, "-- Debug Particle ---"
             print*, "ID: ", id_current
             print*, "Timestep: ", dt
             print*, "rpar: ", rpar
             print*, "Current abundances & energy: ", y
             print*, "ydot: ", ydot
             print*, "-- Debug Particle ---"
          endif
#endif
          
c
c Now compute ydot. On the first step, we must use rate_eq directly,
c but on subsequent steps, we can just call dvindy [if we're using
c the dvode integrator]
c
          if (I .eq. 1) then
            call rate_eq(nspec, t_start, y, ydot, rpar, ipar)
          else
            call dvindy(t_start, 1, rwork(21), nspec, ydot, iflag_dvi)
          endif
c
c Init. dt_est to an unrealistic value
c
          dt_est = 1d18

          do J = 1, nspec
            if (used_for_stepsize_control(J) .eq. 1) then
              if (ydot(J)   .ne. 0d0 .and.
     $            safety(J) .ne. 0d0 .and.
     $            y(J)      .gt. atol(J)) then
                dt_est = min(dt_est, abs(safety(J) * y(J) / ydot(J)))
              endif
            endif
          enddo
c
          ipar(2) = 0 
c
c Don't evolve past the end of the global timestep
c
          dt = min(dt_est, t_end - t_start)
        endif
c
        t_in  = t_start
        t_out = t_start + dt
        if (dt .lt. 0d0) then
          print*, 'Negative timestep!'
          ABORT(1)
        endif
        if (dt .eq. 0d0) then 
          if (t_end .eq. t_start) then
            return
          else
            print*, 'Zero timestep!'
            ABORT(2)
          endif
        endif
c

c
c Debug code -- print out rpar, y and ydot each timestep and each time
c through the loop.
c
#if defined(SGCHEM_DEBUG_PARTICLE)
        if(id_current .eq. SGCHEM_DEBUG_PARTICLE) then
           print *, 'SGCHEM DEBUG id, i, y', id_current, I, y
           print *, 'SGCHEM DEBUG id, i, ydot', id_current, I, ydot 
        end if
#endif


        call dvode(rate_eq, nspec, y, t_in, t_out, itol, rtol, atol, 
     $             itask, istate, iopt, rwork, lrw, iwork, liw, jac, 
     $             mf, rpar, ipar)
c
c Ensure that dvode has produced sane output values -- if not, we
c need to shrink the timestep. Note that if ipar(1) is non-zero, 
c or if istate != 2 (indicating an integrator error), we're going
c to shrink the timestep anyway, so we can skip the check.
c
        if (ipar(1) .eq. 0 .and. istate .eq. 2) then
          call check_abundances(y, ipar)
        endif
c
        if ((ipar(1) .ne. 0) .or. (istate .ne. 2)) then 
#ifdef DEBUG_EVOLVE
         if (istate .ne. 2) then
           print*, "Error: istate = ", istate
         endif
#endif
          do J = 1, nspec
            y(J) = y_old(J)
          enddo
          tdust    = tdust_old
          dt       = dt / 2d0
          if (dt .lt. 1d-10 .and. ntry .gt. 80) then
            print*, "Excessively small timestep! ", dt, dt_est, 
     $              t_start, t_end
            print*, "ID: ", id_current
            print*, "yn",yn,"dl", dl,"divv", divv
            print*, "y",y
	    print*, "ydot",ydot
            print*, "rpar",rpar
c Print values from start of timestep
            print*, "Initial energy: ", energy
            print*, "Initial abundances: ", non_eq_abundances
            print*, "Full timestep: ", time
#ifdef SGCHEM_VARIABLE_Z
            print*, "Carbon abundance", abundc
            print*, "Oxygen abundance", abundo
            print*, "M abundance", abundM
            print*, "Total metallicity", Z_atom
            print*, "Dust-to-gas ratio", dust_to_gas_ratio
#endif
#ifdef SGCHEM_VARIABLE_CRION
            print*, "CR ion rate", cosmic_ray_ion_rate
#endif

#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
            ABORT(100)
#endif
          endif
          new_step = 0
          ipar(1)  = 0
          istate   = 1
          ntry     = ntry + 1
        else
          t_start = t_start + dt
          new_step = 1
          ntry     = 1
        endif
c
        if (t_start .ge. t_end) then
          goto 100
        endif
      enddo
c
      write(*,'(A)') 'EVOLVE_ABUNDANCES: nstep exceeded,' 
      print*, "ID: ", id_current
      print*, "Timestep: ", dt
      print*, "rpar: ", rpar
      print*, "Current abundances & energy: ", y
      print*, "ydot: ", ydot
c Current temperature (assuming gamma = 5/3)
      yntot = (1d0 + ABHE - y(ih2) + y(ihp)) * yn
      temp  = (2d0 / 3d0) * y(itmp) / (yntot * kboltz)
      print*, "Current temperature: ", temp
c Print values from start of timestep
      print*, "Initial energy: ", energy
      print*, "Initial abundances: ", non_eq_abundances
      yntot = (1d0 + ABHE - non_eq_abundances(ih2)
     $      +  non_eq_abundances(ihp)) * yn
      temp  = (2d0 / 3d0) * energy / (yntot * kboltz)
      print*, "Initial temperature: ", temp
      print*, "Full timestep: ", time
#ifdef SGCHEM_VARIABLE_Z
      print*, "Carbon abundance", abundc
      print*, "Oxygen abundance", abundo
      print*, "M abundance", abundM
      print*, "Total metallicity", Z_atom
      print*, "Dust-to-gas ratio", dust_to_gas_ratio
#endif
#ifdef SGCHEM_VARIABLE_CRION
      print*, "CR ion rate", cosmic_ray_ion_rate
#endif
      ABORT(0)
      return

100   continue

      do I = 1, nspec-1
        non_eq_abundances(I) = y(I)
      enddo
      energy = y(itmp)

c Calculate final dust temp.
#if CHEMISTRYNETWORK != 1
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
      abHI   = max(0d0, 1d0 - 2d0 * y(ih2) - y(ihp))
      abhep  = 0d0
#endif
#if CHEMISTRYNETWORK == 7
      abHI   = max(0d0, 1d0 - 2d0 * y(ih2) - y(ihp))
      abhep  = y(ihep)
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
      abHI  = max(0d0, 1d0 - 2d0 * y(ih2) - y(ihp) - y(ichx) - y(ihcop))
      abhep = y(ihep)
#endif
      abe   = y(ihp) + abhep
      yntot = (1d0 + ABHE - y(ih2) + abe) * yn
#ifdef VARIABLE_GAMMA
      ekn = energy / (yntot * kboltz)
      call calc_temp(abh2, ekn, temp)
#else
      temp = (gamma - 1d0) * energy / (yntot * kboltz)
#endif
      uvdust = photorates_direct(6)
      call calc_dust_temp(yn, uvdust, chi_mean, temp, abHI,
     $                    RH2, Rdust, 0)
      dust_temp = tdust
#endif

#endif /* SGCHEM */
      return      
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////     E V O L V E _ A B U N D A N C E S     \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////      C H E C K _ A B U N D A N C E S    \\\\\\\\\\
c
c=======================================================================
c
      subroutine check_abundances(y, ipar)
#ifdef SGCHEM
      implicit NONE
#include "cool.h"
#include "non_eq.h"
      REAL y(nspec), max_abundance(nspec)
      REAL current_carbon, current_oxygen
      REAL diff
      integer ipar(nipar)
      integer J
c
#if CHEMISTRYNETWORK == 1
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
        max_abundance(idp)   = abundD
        max_abundance(ihd)   = abundD
        max_abundance(ihep)  = abhe
        max_abundance(ihepp) = abhe
#endif
c
#if CHEMISTRYNETWORK == 4
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
#endif
c
#if CHEMISTRYNETWORK == 5
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
        max_abundance(ico)   = min(abundc, abundo)
#endif
c
#if CHEMISTRYNETWORK == 7
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
        max_abundance(ico)   = min(abundc, abundo)
        max_abundance(ihep)  = abhe
#endif
c
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
        max_abundance(ih2) = 0.5d0
        max_abundance(ihp) = 1.0d0
        max_abundance(icp)  = abundc
        max_abundance(ichx) = abundc
        max_abundance(iohx) = abundo
        max_abundance(ico)   = min(abundo, abundc)
        max_abundance(ihcop) = min(abundo, abundc)
        max_abundance(ihep) = abhe
        max_abundance(imp)  = abundM  ! Silicon when using network 16 
#endif
c
c Check for NaNs
c
      do J = 1, nspec
        if (y(J) .ne. y(J)) then
#ifdef DEBUG_EVOLVE
          print*, "NaN for ID ", id_current, y
#endif
          ipar(1) = 1
          return
        endif
      enddo
c
c Check for negative abundances. Small negative values (which have
c absolute values less than the corresponding absolute tolerance)
c are accepted, but are clipped to zero. For large negative values,
c we repeat the step with a smaller dt.
c
      do J = 1, nspec
        if (y(J) .lt. 0d0) then
          if (y(J) .lt. -atol(J)) then
#ifdef DEBUG_EVOLVE
            print*, "Negative abundance for ID ", id_current, y
#endif
            ipar(1) = 1
            return
          else
            y(J) = 0d0
          endif
        endif
      enddo
c
c Don't allow ion abundances to exceed corresponding total metal
c abundances. Again, small errors are silently corrected, large
c errors lead to us repeating with a smaller timestep
c
      do J = 1, nspec
        if (J .ne. itmp .and. y(J) .gt. max_abundance(J)) then
          diff = y(J) - max_abundance(J)
          if (diff .lt. RTOL * y(J)) then
            y(J) = max_abundance(J)
          else
#ifdef DEBUG_EVOLVE
            print*, "Excess abundance for ID ", id_current, y
#endif
            ipar(1) = 1
            return
          endif
        endif
      enddo
c
c Add check that He+ and He++ abundances together don't exceed total helium
c abundance
c
#if CHEMISTRYNETWORK == 1
      diff = y(ihep) + y(ihepp) - abhe
      if (diff .gt. 0d0) then
        if (diff .gt. RTOL * abhe) then
#ifdef DEBUG_EVOLVE
          print*, "Excessive ionized helium for ID ", id_current, y
#endif
          ipar(1) = 1
          return
        else
c Rescale to ensure HeI abundance non-negative
          y(ihep)  = y(ihep)  / (1d0 + diff / abhe)
          y(ihepp) = y(ihepp) / (1d0 + diff / abhe)
        endif
      endif
#endif
#endif /* SGCHEM */
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////      C H E C K _ A B U N D A N C E S    \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////       I N I T _ T O L E R A N C E S       \\\\\\\\\\
c
c=======================================================================
c
      subroutine init_tolerances
#ifdef SGCHEM
      implicit NONE
#include "non_eq.h"
#include "cool.h"
c
c Set-up absolute tolerances
c
#if CHEMISTRYNETWORK == 1
        atol(ihp)   = ATOL_HP
        atol(ih2)   = ATOL_H2
        atol(idp)   = ATOL_DP
        atol(ihd)   = ATOL_HD
        atol(ihep)  = ATOL_HEP
        atol(ihepp) = ATOL_HEPP
#endif
c
#if CHEMISTRYNETWORK == 4
        atol(ihp) = ATOL_HP
        atol(ih2) = ATOL_H2
#endif
c
#if CHEMISTRYNETWORK == 5
        atol(ihp) = ATOL_HP
        atol(ih2) = ATOL_H2
        atol(ico) = ATOL_CO
#endif
c
#if CHEMISTRYNETWORK == 7
        atol(ihp) = ATOL_HP
        atol(ih2) = ATOL_H2
        atol(ico) = ATOL_CO
        atol(ihep) = ATOL_HEP
#endif
c
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
        atol(ih2)   = ATOL_H2
        atol(ihp)   = ATOL_HP
        atol(icp)   = ATOL_CP
        atol(ichx)  = ATOL_CHX
        atol(iohx)  = ATOL_OHX
        atol(ico)   = ATOL_CO
        atol(ihcop) = ATOL_HCOP
        atol(ihep)  = ATOL_HEP
        atol(imp)   = ATOL_MP
#endif
c
      atol(itmp)  = ATOL_TMP
c
#endif /* SGCHEM */
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////       I N I T _ T O L E R A N C E S       \\\\\\\\\\
c
c=======================================================================
c
      subroutine check_advection(y)
#ifdef SGCHEM
      implicit NONE
#include "cool.h"
#include "non_eq.h"
      REAL y(nspec), max_abundance(nspec)
      REAL current_carbon, current_oxygen
      REAL diff, tiny, tot_helium
      integer J

      parameter(tiny = 1d-10)
c
#ifdef SGCHEM_VARIABLE_Z
c Check that we've got sane values for total C, O abundances.
c
c We don't both with this check when these are constant, as in this
c case they can only be negative if this is explicitly requested in
c the parameter file.
c
      if (abundc .lt. 0d0) then
        if (abundc .lt. -tiny) then
          print*, "Error: negative elemental abundance of C", abundc
          ABORT(101)
        else
c Allow for a tiny amount of undershooting
          abundc = tiny
        endif
      endif
c
      if (abundo .lt. 0d0) then
        if (abundo .lt. -tiny) then
          print*, "Error: negative elemental abundance of O", abundo
          ABORT(102)
        else
c Allow for a tiny amount of undershooting
          abundo = tiny
        endif
      endif

      if (abundM .lt. 0d0) then
        if (abundM .lt. -tiny) then
          print*, "Error: negative elemental abundance of M", abundo
          ABORT(103)
        else
c Allow for a tiny amount of undershooting
          abundM = tiny
        endif
      endif

      if (Z_atom .lt. 0d0) then
        if (Z_atom .lt. -tiny) then
          print*, "Error: negative elemental abundance of Z", Z_atom
          ABORT(104)
        else
c Allow for a tiny amount of undershooting
          Z_atom = 0d0
        endif
      endif

      if (dust_to_gas_ratio.lt. 0d0) then
        if (dust_to_gas_ratio.lt. -tiny) then
          print*, "Error: negative dust to gas", dust_to_gas_ratio
          ABORT(105)
        else
c Allow for a tiny amount of undershooting
          dust_to_gas_ratio = 0d0
        endif
      endif

#endif /* SGCHEM_VARIABLE_Z */
c
#if CHEMISTRYNETWORK == 1
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
        max_abundance(idp)   = abundD
        max_abundance(ihd)   = abundD
        max_abundance(ihep)  = abhe
        max_abundance(ihepp) = abhe
#endif
c
#if CHEMISTRYNETWORK == 4
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
#endif
c
#if CHEMISTRYNETWORK == 5
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
        max_abundance(ico)   = min(abundc, abundo)
#endif
c
#if CHEMISTRYNETWORK == 7
        max_abundance(ih2)   = 0.5d0
        max_abundance(ihp)   = 1.0d0
        max_abundance(ico)   = min(abundc, abundo)
        max_abundance(ihep)  = abhe
#endif
c
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
        max_abundance(ih2) = 0.5d0
        max_abundance(ihp) = 1.0d0
        max_abundance(icp)  = abundc
        max_abundance(ichx) = abundc
        max_abundance(iohx) = abundo
        max_abundance(ico)   = min(abundo, abundc)
        max_abundance(ihcop) = min(abundo, abundc)
        max_abundance(ihep) = abhe
        max_abundance(imp)  = abundM
#endif
c
c
c Check for negative abundances. Small negative values (which have
c absolute values less than the corresponding absolute tolerance)
c are accepted, but are clipped to zero. For large negative values,
c we repeat the step with a smaller dt.
c
      do J = 1, nspec
        if (y(J) .lt. 0d0) then
          if (y(J) .lt. -atol(J)) then
            print*, "Negative abundance for ID greater than tolerance",
     &            id_current, y, J, atol(J)
            y(J) = 0d0
          else
            y(J) = 0d0
          endif
        endif
      enddo
c
c Don't allow ion abundances to exceed corresponding total metal
c abundances. Again, small errors are silently corrected, large
c errors lead to us repeating with a smaller timestep
c
      do J = 1, nspec
        if (J .ne. itmp .and. y(J) .gt. max_abundance(J)) then
          diff = y(J) - max_abundance(J)
          if (diff .lt. RTOL * y(J)) then
            y(J) = max_abundance(J)
          else
#ifdef DEBUG_EVOLVE
             print*, "Excess abundance for ID greater than tolerance", 
     &            id_current, y, J, max_abundance(J)
#endif
             y(J) = max_abundance(J)
          endif
       endif
      enddo

#if CHEMISTRYNETWORK == 1
      tot_helium = y(ihep) + y(ihepp)
#ifdef DEBUG_EVOLVE
      if (tot_helium .gt. abhe + tiny) then
c Only print the warning if the excess is more than a tiny amount
        print*, "Excess total helium abundance", id_current,
     $          y(ihep), y(ihepp), tot_helium
      endif
#endif
c
      diff = tot_helium - abhe
      if (diff .gt. 0d0) then
c Rescale to ensure HeI abundance non-negative
        y(ihep)  = y(ihep)  / (1d0 + diff / abhe)
        y(ihepp) = y(ihepp) / (1d0 + diff / abhe)
      endif
#endif

#endif /* SGCHEM */
      return
      end

      subroutine compute_Jeans_length(energy, yn, gamma, dl_Jeans)
#ifdef SGCHEM
#include "cool.h"
      REAL energy, yn, gamma, dl_Jeans
      REAL rho, cs2
      REAL G_n
      parameter (G_n = 6.672e-8)
      REAL pi
      parameter (pi = 3.14159265358979323846)

      rho = (1d0 + 4d0 * abhe) * mh * yn
      cs2 = (gamma * (gamma - 1d0) * energy / rho)
      dl_Jeans  = dsqrt(pi * cs2 / (G_n * rho))
#endif
      return
      end      

