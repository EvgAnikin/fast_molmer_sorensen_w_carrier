      module alpha_chi_expr
          implicit none
          private
          integer, parameter :: kmax = 10
          double precision, parameter :: et_min = 0.4
          public :: alpha_tilde_func, x_tilde_func, ipowexp, ipowexp_a, ipowexp2_a, factorial, et_min

      contains

      function factorial(n)
        implicit none
        integer n, k, f, factorial

        f = 1
        do k = 1, n
          f = f*k
        end do
        factorial = f

      end function

      function ipowexp_a(p, e, t) result(vals)
        implicit none
        integer p, pp, k
        double precision e, t
        double complex vals(p+1), val


        if (abs(e*t) < et_min) then
          do pp = 0, p
            val = 0
            do k = kmax, 0, -1
              val = val + ((0,1)*e)**k*t**(pp + k + 1)/(pp + k + 1)/factorial(k)
            end do
            vals(pp + 1) = val
          end do
        else
          vals(1) = (exp((0,1)*e*t) - 1)/((0,1)*e)

          do pp = 1, p
            vals(pp+1) = t**pp*exp((0,1)*e*t)/((0,1)*e) - pp/((0,1)*e)*vals(pp)
          end do
        end if
      end function

      function ipowexp(p, e, t) result(val)
        implicit none
        integer p, j
        double precision e, t
        double complex vals(p+1)
        double complex val

        vals = ipowexp_a(p, e, t)
        val  = vals(p+1)

      end function

      function ipowexp2_a(p1, p2, e1, e2, t) result(vals)
        implicit none
        integer p1, p2, pp1, pp2, k
        double precision e1, e2, t
        double complex vals(p1+1, p2+1), ipowexp_vals_12(p1+p2+1), ipowexp_vals_1(p1+1), ipowexp_vals_taylor(p1+p2+kmax+2), val

        if (abs(e2*t) < et_min) then
          ipowexp_vals_taylor = ipowexp_a(p1 + p2 + kmax + 1, e1, t)
          do pp1 = 0, p1
            do pp2 = 0, p2
              val = 0
              do k = kmax, 0, -1
                val = val + ((0,1)*e2)**k/factorial(k)/(pp2 + k  + 1)*ipowexp_vals_taylor(pp1 + pp2 + k + 2)
!                val = val + ((0,1)*e2)**k/factorial(k)/(pp2 + k  + 1)*ipowexp(pp1 + pp2 + k + 1, e1, t)
              end do
              vals(pp1+1, pp2+1) = val
!              vals(pp1 + 1, pp2 + 1) = 1/(pp2 + 1)*ipowexp(pp1 + pp2 + 1, e1, t)
            end do
          end do
        else
          ipowexp_vals_1 = ipowexp_a(p1, e1, t)
          ipowexp_vals_12 = ipowexp_a(p1 + p2, e1 + e2, t)
          do pp1 = 0, p1 
            vals(pp1 + 1, 1) = 1/((0,1)*e2) *(ipowexp_vals_12(pp1 + 1) - ipowexp_vals_1(pp1 + 1))
!            vals(pp1 + 1, 1) = 1/((0,1)*e2) *(ipowexp(pp1, e1 + e2, t) - ipowexp(pp1, e1, t))
            do pp2 = 1, p2
!              vals(pp1 + 1, pp2 + 1) = 1/((0,1)*e2)*ipowexp(pp1 + pp2, e1 + e2, t) - pp2/((0,1)*e2)*vals(pp1 + 1, pp2)
              vals(pp1 + 1, pp2 + 1) = 1/((0,1)*e2)*ipowexp_vals_12(pp1 + pp2 + 1) - pp2/((0,1)*e2)*vals(pp1 + 1, pp2)
            end do
          end do
        end if
      end function

      function alpha_tilde_func(t2, t1, mu, omega, psi, pmax)
          implicit none
          integer :: pmax, j
          double precision :: t2, t1, mu, omega, psi
          double complex :: I
          double complex :: alpha_tilde_func_r(pmax+1), alpha_tilde_func(pmax+1)
          I = (0, 1)

          alpha_tilde_func_r = -I/2*exp(I*(omega + mu)*t1 + I*psi)*ipowexp_a(pmax, mu + omega, t2 - t1) - I/2*exp(I*(omega - mu)*t1 - I*psi)*ipowexp_a(pmax, omega - mu, t2 - t1) 

          do j = 1, pmax+1
            alpha_tilde_func(j) = alpha_tilde_func_r(pmax + 2 - j)
          end do

       end function
      
      function x_tilde_func(t2, t1, mu, omega, psi, psi_prime, pmax)
          implicit none
          integer :: pmax, j, k
          double precision t2, t1, mu, omega, psi, psi_prime
          double complex :: I
          double complex :: x_tilde_func(pmax+1, pmax+1), x_tilde_func_rt(pmax+1, pmax+1)

          I = (0, 1)

          x_tilde_func_rt = 0.25*(exp(2*I*mu*t1 + I*(psi + psi_prime))*ipowexp2_a(pmax, pmax, mu - omega, omega + mu, t2 - t1) + &
                                  exp(I*(psi_prime - psi))*ipowexp2_a(pmax, pmax, -mu - omega, omega + mu, t2 - t1) + &
                                  exp(-I*(psi_prime - psi))*ipowexp2_a(pmax, pmax, mu - omega, omega - mu, t2 - t1) + &
                                  exp(-2*I*mu*t1 - I*(psi + psi_prime))*ipowexp2_a(pmax, pmax, - mu - omega, omega - mu, t2 - t1))

          do j = 1, pmax + 1
            do k = 1, pmax + 1
              x_tilde_func(pmax + 2 - k, pmax + 2 - j) = x_tilde_func_rt(j, k)
            end do
          end do

      end function

      end module
