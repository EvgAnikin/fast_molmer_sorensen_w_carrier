      subroutine ppoly_carrier_phase_eqf_at_bp(phase_at_bp, t_bp, oomegas, n_seg, mu, psi, pmax)
        use ppoly_alpha_chi_func
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t_bp(n_seg+1), mu, psi
        integer pmax, j, n_seg
        double precision phase_at_bp(n_seg+1)

!F2PY  intent(out)  :: phase_at_bp
!F2PY  intent(hide) :: pmax, n_seg

        phase_at_bp = ppoly_carrier_phase_eqf_at_bp_func(oomegas, t_bp, n_seg, mu, psi, pmax)
      end 

      subroutine ppoly_carrier_phase_eqf(phase_func, t2, t1, t_bp, oomegas, mu, psi, n_seg, carrier_phases_at_bp, pmax)
        use ppoly_alpha_chi_func
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t2, t1, t_bp(n_seg + 1), mu, psi
        integer n_seg, pmax, j, k
        double precision carrier_phases_at_bp(n_seg + 1)
        double precision phase_func

!F2PY  intent(out)  :: phase_func
!F2PY  intent(hide) :: pmax, n_seg

        phase_func = ppoly_carrier_phase_eqf_func(oomegas, t2, t1, t_bp, mu, psi, n_seg, carrier_phases_at_bp, pmax)
      end


      subroutine ppoly_alpha_eqf_at_bp(alphas_at_bp, t_bp, oomegas, mu, psi, omegas, n_seg, n_modes, pmax)
!
        use ppoly_alpha_chi_func
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alphas_at_bp(n_seg + 1, n_modes)

!F2PY  intent(out)  :: alphas_at_bp
!F2PY  intent(hide) :: pmax, n_seg, n_modes

        alphas_at_bp = ppoly_alpha_eqf_at_bp_func(oomegas, t_bp, mu, psi, omegas, n_seg, n_modes, pmax)
      end

      subroutine ppoly_alpha_eqf(alphas, t2, t1, t_bp, oomegas, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, pmax)
        use ppoly_alpha_chi_func
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t2, t1, t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alphas_at_bp(n_seg + 1, n_modes)
        double complex alphas(n_modes)
  

!F2PY  intent(out)  :: alphas
!F2PY  intent(hide) :: pmax, n_seg, n_modes

        alphas = ppoly_alpha_eqf_func(oomegas, t2, t1, t_bp, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, pmax)

      end 


      subroutine ppoly_chi_eqf_at_bp(chis_at_bp, t_bp, oomegas, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, pmax)
        use ppoly_alpha_chi_func
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alphas_at_bp(n_seg + 1, n_modes)
        double precision chis_at_bp(n_seg + 1, n_modes)

!F2PY  intent(out)  :: chis_at_bp
!F2PY  intent(hide) :: pmax, n_seg, n_modes

        chis_at_bp = ppoly_chi_eqf_at_bp_func(oomegas, t_bp, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, pmax)
      end

      subroutine ppoly_chi_eqf(chis, t2, t1, t_bp, oomegas, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, chis_at_bp, pmax)
        use ppoly_alpha_chi_func
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t2, t1
        double precision t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alphas_at_bp(n_seg + 1, n_modes)
        double precision chis(n_modes), chis_at_bp(n_seg + 1, n_modes)

!F2PY  intent(out)  :: chis
!F2PY  intent(hide) :: pmax, n_seg, n_modes

        chis = ppoly_chi_eqf_func(oomegas, t2, t1, t_bp, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, chis_at_bp, pmax)
      end
