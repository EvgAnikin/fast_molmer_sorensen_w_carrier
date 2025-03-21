      module ppoly_alpha_chi_func
          implicit none
          public :: poly_carrier_phase_eqf_func, poly_alpha_eqf_func, ppoly_carrier_phase_eqf_at_bp_func, &
            ppoly_carrier_phase_eqf_func0, ppoly_carrier_phase_eqf_func, ppoly_alpha_eqf_at_bp_func, &
            ppoly_alpha_eqf_func, ppoly_chi_eqf_at_bp_func, ppoly_chi_eqf_func0, poly_chi_eqf_func, tensordot, segment_number
      contains

      function poly_carrier_phase_eqf_func(oomegas, t2, t1, mu, psi, pmax) result(phase_func)
        use alpha_chi_expr
        implicit none
        double precision oomegas(pmax+1)
        double precision t2, t1, mu, psi
        integer pmax
        double complex alpha_tilde(pmax+1)
        double precision phase_func

        alpha_tilde = (0,1)*alpha_tilde_func(t2, t1, mu, 0D0, psi, pmax)
        phase_func = real(dot_product(oomegas, alpha_tilde))

      end function

      function poly_alpha_eqf_func(oomegas, t2, t1, mu, psi, n_modes, omegas, pmax) result(alphas)
         use alpha_chi_expr
         implicit none
         double precision oomegas(pmax+1)
         double precision t2, t1, mu, psi
         integer pmax, n_modes
         double precision omegas(n_modes)
         double complex alphas(n_modes)
         double complex alpha_tilde(pmax+1)
         integer m

         do m = 1, n_modes
            alpha_tilde = alpha_tilde_func(t2, t1, mu, omegas(m), psi, pmax)
            alphas(m) = dot_product(oomegas, alpha_tilde)
         end do
      end

      function poly_chi_eqf_func(oomegas, t2, t1, mu, psi, n_modes, omegas, pmax) result(chis)
         use alpha_chi_expr
         implicit none
         double precision oomegas(pmax+1)
         double precision t2, t1, mu, psi
         integer pmax, n_modes
         double precision omegas(n_modes)
         double precision chis (n_modes)
         double complex x_part(pmax+1, pmax+1)
         integer j, k, m

         chis = 0

         do m = 1, n_modes
            x_part = x_tilde_func(t2, t1, mu, omegas(m), psi, psi, pmax)
            chis(m) = dot_product(oomegas, matmul(aimag(x_part + transpose(x_part)), oomegas))
         end do
      end

      function segment_number(t, t_bp, n_seg) ! Returns the number of segment numbered from 1
        implicit none 
        double precision t
        integer n_seg
        double precision t_bp(n_seg + 1)
        integer segment_number, i

        if (t < t_bp(1)) then
            segment_number = 1 
            return 
        end if 

        do i = 1, n_seg
          if (t < t_bp(i+1)) then
            segment_number = i
            return 
          end if 
        end do

        if (t >= t_bp(n_seg + 1)) then 
          segment_number = n_seg
          return
        end if 
      end


      function ppoly_oomega_func(t, t_bp, oomegas, n_seg, pmax) result(oomega)
        implicit none 
        double precision oomega, oomegas(pmax+1, n_seg)
        double precision t, t_bp(n_seg + 1)
        integer n_seg, pmax, p, k

        k = segment_number(t, t_bp, n_seg)

        oomega = 0
        do p = 0, pmax 
          oomega = oomega + oomegas(pmax + 1 - p, k)*(t - t_bp(k))**p
        end do

        return 
      end function


      function ppoly_carrier_phase_eqf_at_bp_func(oomegas, t_bp, n_seg, mu, psi, pmax) result(phase_at_bp)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t_bp(n_seg+1), mu, psi
        integer pmax, j, n_seg
        double precision phase_at_bp(n_seg+1)

        phase_at_bp(1) = 0

        do j = 1, n_seg
          phase_at_bp(j+1) = phase_at_bp(j) + poly_carrier_phase_eqf_func(oomegas(:, j), t_bp(j+1), t_bp(j), mu, psi, pmax)
        end do
      end function


      function ppoly_carrier_phase_eqf_func0(oomegas, t, t_bp, mu, psi, n_seg, carrier_phases_at_bp, pmax) result(phase_func)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t, t_bp(n_seg + 1), mu, psi
        integer n_seg, pmax, j
        double precision carrier_phases_at_bp(n_seg + 1)
        double precision d_phase
        double precision phase_func

        j = segment_number(t, t_bp, n_seg)

        d_phase = poly_carrier_phase_eqf_func(oomegas(:, j), t, t_bp(j), mu, psi, pmax)

        phase_func = d_phase + carrier_phases_at_bp(j)
      end

      function ppoly_carrier_phase_eqf_func(oomegas, t2, t1, t_bp, mu, psi, n_seg, carrier_phases_at_bp, pmax) result(phase_func)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t2, t1, t_bp(n_seg + 1), mu, psi
        integer n_seg, pmax, j, k
        double precision carrier_phases_at_bp(n_seg + 1)
        double precision d_phase_2, d_phase_1
        double precision phase_func

        j = segment_number(t1, t_bp, n_seg)
        k = segment_number(t2, t_bp, n_seg)

        d_phase_1 = poly_carrier_phase_eqf_func(oomegas(:, j), t1, t_bp(j), mu, psi, pmax)
        d_phase_2 = poly_carrier_phase_eqf_func(oomegas(:, k), t2, t_bp(k), mu, psi, pmax)

        phase_func = d_phase_2 + carrier_phases_at_bp(k) - d_phase_1 - carrier_phases_at_bp(j)
      end


      function ppoly_alpha_eqf_at_bp_func(oomegas, t_bp, mu, psi, omegas, n_seg, n_modes, pmax) result(alphas_at_bp)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alphas_at_bp(n_seg + 1, n_modes)
        integer j

        alphas_at_bp(1, :) = 0
        do j = 1, n_seg
          alphas_at_bp(j+1, :) = alphas_at_bp(j, :) + poly_alpha_eqf_func(oomegas(:, j), t_bp(j+1), t_bp(j), mu, psi, n_modes, omegas, pmax)
        end do
      end 

      function ppoly_alpha_eqf_func(oomegas, t2, t1, t_bp, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, pmax) result(alphas)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t1, t2, t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alphas_at_bp(n_seg + 1, n_modes)
        double complex alphas(n_modes)
        double complex dalpha1(n_modes), dalpha2(n_modes)
        integer j,k
 
        j = segment_number(t1, t_bp, n_seg)
        k = segment_number(t2, t_bp, n_seg)
        dalpha1 = poly_alpha_eqf_func(oomegas(:, j), t1, t_bp(j), mu, psi, n_modes, omegas, pmax)
        dalpha2 = poly_alpha_eqf_func(oomegas(:, k), t2, t_bp(k), mu, psi, n_modes, omegas, pmax)
        alphas = alphas_at_bp(k, :) - alphas_at_bp(j, :) + dalpha2 - dalpha1
      end 

      function tensordot(a, b)
        implicit none
        double precision a(:), b(:)
        double precision :: tensordot(size(a),size(b))
        integer j,k

        do j = 1, size(a)
          do k = 1, size(b)
            tensordot(j,k) = a(j)*b(k)       
          end do
        end do
      end function

      function ppoly_chi_eqf_at_bp_func(oomegas, t_bp, mu, psi, omegas, n_seg, n_modes, alphas_at_bp, pmax) result(chis_at_bp)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alphas_at_bp(n_seg + 1, n_modes)
        double complex alpha_prod(n_modes)
        double precision chis_at_bp(n_seg + 1, n_modes)
        double precision dchi(n_modes) 
        integer j

        chis_at_bp(1, :) = 0
        do j = 1, n_seg
          dchi = poly_chi_eqf_func(oomegas(:, j), t_bp(j+1), t_bp(j), mu, psi, n_modes, omegas, pmax)
          alpha_prod = alphas_at_bp(j, :)*conjg(alphas_at_bp(j+1, :) - alphas_at_bp(j, :))
          chis_at_bp(j+1, :) = chis_at_bp(j, :) + dchi + 2*aimag(alpha_prod)
        end do
      end

      function ppoly_chi_eqf_func0(oomegas, t, t_bp, mu, psi, omegas, n_seg, n_modes, alpha_at_bp, chi_at_bp, pmax) result(chi)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t
        double precision t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alpha_at_bp(n_seg + 1, n_modes)
        double precision chi_at_bp(n_seg + 1, n_modes)
        double complex dalpha(n_modes)
        double complex alpha_prod(n_modes)
        double precision chi(n_modes)
        double precision dchi(n_modes)
        integer j

        j = segment_number(t, t_bp, n_seg)
        dalpha = poly_alpha_eqf_func(oomegas(:, j), t, t_bp(j), mu, psi, n_modes, omegas, pmax)
        dchi   =   poly_chi_eqf_func(oomegas(:, j), t, t_bp(j), mu, psi, n_modes, omegas, pmax)

        alpha_prod = (dalpha + alpha_at_bp(j, :))*conjg(dalpha)
        chi = chi_at_bp(j,:) + dchi + 2*aimag(alpha_prod)

      end

      function ppoly_chi_eqf_func(oomegas, t2, t1, t_bp, mu, psi, omegas, n_seg, n_modes, alpha_at_bp, chi_at_bp, pmax) result(chi)
        implicit none
        double precision oomegas(pmax+1, n_seg)
        double precision t2, t1
        double precision t_bp(n_seg + 1)
        double precision mu
        integer pmax, n_seg, n_modes
        double precision psi
        double precision omegas(n_modes)
        double complex alpha_at_bp(n_seg + 1, n_modes)
        double precision chi_at_bp(n_seg + 1, n_modes)
        double complex alpha10(n_modes), alpha21(n_modes)
        double complex  alpha_prod(n_modes)
        double precision chi1(n_modes), chi2(n_modes), chi(n_modes)

        alpha10 = ppoly_alpha_eqf_func(oomegas, t1, t_bp(1), t_bp, mu, psi, omegas, n_seg, n_modes, alpha_at_bp, pmax)
        alpha21 = ppoly_alpha_eqf_func(oomegas, t2, t1,      t_bp, mu, psi, omegas, n_seg, n_modes, alpha_at_bp, pmax)

        alpha_prod = alpha10*conjg(alpha21)
        chi1 = ppoly_chi_eqf_func0(oomegas, t1, t_bp, mu, psi, omegas, n_seg, n_modes, alpha_at_bp, chi_at_bp, pmax)
        chi2 = ppoly_chi_eqf_func0(oomegas, t2, t_bp, mu, psi, omegas, n_seg, n_modes, alpha_at_bp, chi_at_bp, pmax)

        chi = chi2 - chi1 - 2*aimag(alpha_prod)
      end

      end module
