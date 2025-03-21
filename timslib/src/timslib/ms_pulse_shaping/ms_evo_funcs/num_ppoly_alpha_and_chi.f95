       subroutine num_alpha_and_chi_eqf_w_carrier(alpha_and_chi, n_vals, oomegas, t_bp, mu, psi, omegas, carrier_phases_at_bp, n_seg, n_modes, pmax) 
         use ppoly_alpha_chi_func
         implicit none
         double precision oomegas(pmax+1, n_seg)
         double precision t_bp(n_seg + 1)
         double precision mu
         integer pmax, n_seg, n_modes, n_vals, j
         double precision psi
         double precision omegas(n_modes)
         double precision t, t_i, dt
         double precision oomega_1, oomega_2, oomega_4, phi_1, phi_2, phi_4
         double precision carrier_phases(n_vals), carrier_phases_at_bp(n_seg + 1)
         double precision t_range(n_vals)
         double complex f1(n_modes), f2(n_modes), f4(n_modes), alpha(n_vals, n_modes), alpha_prev(n_modes)
         double precision chi(n_vals, n_modes)
         double complex alpha_and_chi(2, n_vals, n_modes)
         double complex i 

!F2PY intent(out)  :: alpha_and_chi
!F2PY intent(hide) :: n_seg, n_modes, pmax

         i = (0, 1)
         dt = (t_bp(n_seg + 1) - t_bp(1))/(n_vals - 1)

         alpha = 0
         chi = 0

         oomega_4 = 0
         phi_4 = 0
         f4 = oomega_4*cos(mu*t_bp(1) + psi)*cos(2*phi_4)*exp(i*omegas*(t_bp(1)))
         do j = 2, n_vals
           t = t_bp(1) + (j - 2)*dt
           oomega_1 = oomega_4
           oomega_2 = ppoly_oomega_func(t+dt/2, t_bp, oomegas, n_seg, pmax)
           oomega_4 = ppoly_oomega_func(t+dt,   t_bp, oomegas, n_seg, pmax)

           phi_1 = phi_4
           phi_2 = ppoly_carrier_phase_eqf_func0(oomegas, t + dt/2, t_bp, mu, psi, n_seg, carrier_phases_at_bp, pmax) 
           phi_4 = ppoly_carrier_phase_eqf_func0(oomegas, t + dt, t_bp, mu, psi, n_seg, carrier_phases_at_bp, pmax) 

           f1 = f4!oomega_1*cos(mu*t        + psi)*cos(2*phi_1)*exp(i*omegas*t)
           f2 = oomega_2*cos(mu*(t+dt/2) + psi)*cos(2*phi_2)*exp(i*omegas*(t+dt/2))
           f4 = oomega_4*cos(mu*(t+dt  ) + psi)*cos(2*phi_4)*exp(i*omegas*(t+dt  ))

           alpha_prev = alpha(j-1,:)
           alpha(j,:) = alpha_prev - I/6*(f1 + 4*f2 + f4)*dt

           chi(j,:) = chi(j-1,:) + dt/3*real(alpha_prev*conjg(f1)& 
                               + 2*(alpha_prev - I*dt/2*f1)*conjg(f2) &
                               + 2*(alpha_prev - I*dt/2*f2)*conjg(f2) &
                               +   (alpha_prev - I*dt  *f2)*conjg(f4))
         end do

         alpha_and_chi(1, :, :) = alpha
         alpha_and_chi(2, :, :) = chi

      end


       subroutine num_alpha_and_chi_eqf(alpha_and_chi, n_vals, oomegas, t_bp, mu, psi, omegas, n_seg, n_modes, pmax) 
         use ppoly_alpha_chi_func
         implicit none
         double precision oomegas(pmax+1, n_seg)
         double precision t_bp(n_seg + 1)
         double precision mu
         integer pmax, n_seg, n_modes, n_vals, j
         double precision psi
         double precision omegas(n_modes)
         double precision t, t_i, dt
         double precision oomega_1, oomega_2, oomega_4, phi_1, phi_2, phi_4
         double precision t_range(n_vals)
         double complex f1(n_modes), f2(n_modes), f4(n_modes), alpha(n_vals, n_modes), alpha_prev(n_modes)
         double precision chi(n_vals, n_modes)
         double complex alpha_and_chi(2, n_vals, n_modes)
         double complex i 

!F2PY intent(out)  :: alpha_and_chi
!F2PY intent(hide) :: n_seg, n_modes, pmax

         i = (0, 1)
         dt = (t_bp(n_seg + 1) - t_bp(1))/(n_vals - 1)

         alpha = 0
         chi = 0

         oomega_4 = 0
         phi_4 = 0
         do j = 2, n_vals
           t = t_bp(1) + (j - 2)*dt
           oomega_1 = oomega_4
           oomega_2 = ppoly_oomega_func(t+dt/2, t_bp, oomegas, n_seg, pmax)
           oomega_4 = ppoly_oomega_func(t+dt,   t_bp, oomegas, n_seg, pmax)

           f1 = oomega_1*cos(mu*t        + psi)*exp(i*omegas*t)
           f2 = oomega_2*cos(mu*(t+dt/2) + psi)*exp(i*omegas*(t+dt/2))
           f4 = oomega_4*cos(mu*(t+dt  ) + psi)*exp(i*omegas*(t+dt  ))

           alpha_prev = alpha(j-1,:)
           alpha(j,:) = alpha_prev - I/6*(f1 + 4*f2 + f4)*dt

           chi(j,:) = chi(j-1,:) + dt/3*real(alpha_prev*conjg(f1)& 
                               + 2*(alpha_prev - I*dt/2*f1)*conjg(f2) &
                               + 2*(alpha_prev - I*dt/2*f2)*conjg(f2) &
                               +   (alpha_prev - I*dt  *f2)*conjg(f4))
         end do

         alpha_and_chi(1, :, :) = alpha
         alpha_and_chi(2, :, :) = chi

      end
