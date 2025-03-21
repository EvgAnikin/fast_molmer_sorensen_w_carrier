


      function ppoly_chi_multitone_eqf_at_bp_func(oomegas_arr, t_bp, mu_arr, psi_arr, omegas, n_seg, n_modes, n_tones, alphas_at_bp,
        pmax) result (chis_at_bp)
        implicit none
        double precision oomegas_arr(n_tones, pmax+1, n_seg) ! array of pulse envelopes for each tone
        double precision t_bp(n_seg + 1)
        double precision mu_arr(n_tones)
        double precision psi(n_tones)
        integer pmax, n_seg, n_modes, n_tones
        double precision omegas(n_modes)  ! normal mode frequencies
         

      end
