
      subroutine spin_flip_prob(val, s, t_vals, n_ions, n_modes,  &
                        oomega_vals, eta, omegas, mu, psi, &
                        carrier_phases, alphas_tilde, chis_tilde, n_vals) 

          implicit none
          double precision val
          integer s(n_ions), s1(n_ions), n_ions, n_modes, n_vals, delta_m1_m2
          double precision mu, psi, eta(n_ions, n_modes), omegas(n_modes)
          double complex  alphas_tilde(n_vals, n_modes)
          double precision chis_tilde(n_vals, n_modes)
          double precision t_vals(n_vals)
          double precision oomega_vals(n_vals)
          double precision carrier_phases(n_vals)
          double complex alphas(n_vals, n_ions, n_modes)
          double precision chis(n_vals, n_ions, n_ions)
          double complex integrand_vals(n_vals, n_vals)
          double precision oomega1, oomega2, phase1, phase2, chi_t1(n_ions, n_ions), chi_t2(n_ions, n_ions)
          double complex alpha_t1(n_ions, n_modes), alpha_t2(n_ions, n_modes), dalpha(n_modes)
          double complex alpha_t1_s(n_modes), alpha_t2_s(n_modes), alpha_1(n_modes), alpha_2(n_modes), alpha_3(n_modes)
          integer i, j, k, m, m1, m2
          double precision dt, t1, t2

!F2PY intent(out)  :: val
!F2PY intent(hide) :: n_ions, n_modes, n_vals

         dt = t_vals(2) - t_vals(1)

         do i = 1, n_vals
           do j = 1, n_ions
             do m = 1, n_modes
               alphas(i, j, m) = eta(j, m)*alphas_tilde(i, m)
             end do
           end do
         end do

         do i = 1, n_vals
           do j = 1, n_ions
             do k = 1, n_ions
               chis(i, j, k) = sum(eta(j, :)*chis_tilde(i, :)*eta(k, :))
             end do
           end do
         end do

         integrand_vals = 0
         do i = 1, n_vals
           do j = 1, n_vals
             t1 = t_vals(i)
             t2 = t_vals(j)

             oomega1 = oomega_vals(i)
             oomega2 = oomega_vals(j)

             phase1 = carrier_phases(i)
             phase2 = carrier_phases(j)

             alpha_t1 = alphas(i, :, :)
             alpha_t2 = alphas(j, :, :)

             alpha_t1_s = matmul(transpose(alpha_t1), s)
             alpha_t2_s = matmul(transpose(alpha_t2), s)

             chi_t1 = chis(i, :, :)
             chi_t2 = chis(j, :, :)

             alpha_1 = alpha_t1_s
             alpha_3 = alpha_t2_s

             do k = 1, n_ions
               dalpha = alpha_t2(k, :) - alpha_t1(k, :) 
               s1 = s
               s1(k) = -s(k)

               alpha_2 = matmul(transpose(alpha_t1 - alpha_t2), s1)

               do m1 = 1, n_modes
                 do m2 = 1, n_modes
                   if (m1 == m2) then
                     delta_m1_m2 = 1
                   else
                     delta_m1_m2 = 0
                   end if
!                   integrand_vals(i, j) = integrand_vals(i, j) + exp(-20*dot_product(dalpha, dalpha))
                   integrand_vals(i, j) = integrand_vals(i, j) + eta(k, m1)*eta(k, m2)*oomega1*oomega2*cos(mu*t1 + psi)&
                     *cos(mu*t2 + psi)*sin(2*phase1)*sin(2*phase2) &
                     *exp((0, 2)*(s(k)*dot_product(chi_t1(k,:), s) - chi_t1(k,k) &
                                   -s(k)*dot_product(chi_t2(k,:), s) + chi_t2(k,k))) &
                     *exp(-2*dot_product(dalpha, dalpha)) &
                     *exp((0, 4)*aimag(dot_product(alpha_t1(k,:), alpha_t2(k,:))))  &
                     *exp((0, 2)*s(k)*aimag(dot_product(alpha_t1(k, :), alpha_t1_s) - dot_product(alpha_t2(k, :), alpha_t2_s)))&  
                     *(exp(-(0,1)*omegas(m1)*t1 - (0,1)*omegas(m2)*t2)*(alpha_2(m1) + alpha_3(m1))*alpha_3(m2) &
                     + exp( (0,1)*omegas(m1)*t1 - (0,1)*omegas(m2)*t2)*conjg(alpha_1(m1))*alpha_3(m2) &
                     + exp(-(0,1)*omegas(m1)*t1 + (0,1)*omegas(m2)*t2) &
                         *(delta_m1_m2 + conjg(alpha_1(m2) - alpha_2(m2))*(alpha_2(m1) + alpha_3(m1))) &
                     + exp( (0,1)*omegas(m1)*t1 + (0,1)*omegas(m2)*t2)*conjg(alpha_1(m1)*(alpha_1(m2) - alpha_2(m2))))
                 end do
               end do
             end do
           end do
         end do

         val = real(dt*dt*sum(integrand_vals))

      end

