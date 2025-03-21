!
      subroutine alpha_seg_matrix(a_seg_mat, t_bp, n_seg, mu, psi, n_modes, omegas, pmax)
!
         use alpha_chi_expr
         implicit none
         double precision t_bp(n_seg+1)
         double precision mu
         integer pmax, n_seg, n_modes
         double precision psi
         double precision omegas(n_modes)
         double complex a_seg_mat(n_modes, pmax+1, n_seg)
         double complex alpha_tilde(pmax+1)
         integer j, m

!F2PY  intent(out)  :: a_seg_mat
!F2PY  intent(hide) :: n_modes, n_seg

         do j = 1, n_seg
            do m = 1, n_modes
                 alpha_tilde = alpha_tilde_func(t_bp(j+1), t_bp(j), mu, omegas(m), psi, pmax)
                 a_seg_mat(m, :, j) = alpha_tilde
             end do
         end do

         return
      end

!
      subroutine x_seg_matrix(x_seg_mat, t_bp, n_seg, eta1, eta2, mu, psi, omegas, n_modes, pmax)
!
         use alpha_chi_expr
         implicit none
         double precision t_bp(n_seg+1)
         double precision mu
         integer pmax, n_seg, n_modes
         double precision psi
         double precision omegas(n_modes), eta1(n_modes), eta2(n_modes)
         double precision dx
         double precision dx_mat(pmax+1,pmax+1)
         double precision x_seg_mat(pmax+1, n_seg, pmax+1, n_seg)
         double complex alpha_tilde_arr(n_seg, n_modes, pmax+1)
         double complex alpha_tilde_1(pmax+1), alpha_tilde_2(pmax+1)
         integer j, k, m, p2, p1

!F2PY  intent(out)  :: x_seg_mat
!F2PY  intent(hide) :: n_modes, n_seg

         x_seg_mat = 0
         do j = 1, n_seg 
            do m = 1, n_modes
                dx_mat = aimag(eta1(m)*eta2(m)*x_tilde_func(t_bp(j+1), t_bp(j), mu, omegas(m), psi, psi, pmax))
                do p2 = 1, pmax+1
                  do p1 = 1, pmax+1
!                    x_seg_mat(:, j, :, j) = x_seg_mat(:, j, :, j) + dx_mat + transpose(dx_mat)
                    x_seg_mat(p2, j, p1, j) = x_seg_mat(p2, j, p1, j) + dx_mat(p2, p1) + dx_mat(p1, p2)
                  end do
                end do
            end do
         end do

         do j = 1, n_seg
           do m = 1, n_modes
             alpha_tilde_arr(j,m,:) = alpha_tilde_func(t_bp(j+1), t_bp(j), mu, omegas(m), psi, pmax)
           end do
         end do

         ! j < k
         do k = 1, n_seg
            do j = 1, k - 1
                do m = 1, n_modes
                    alpha_tilde_2 = alpha_tilde_arr(j,m,:)
                    alpha_tilde_1 = alpha_tilde_arr(k,m,:)

                    do p2 = 1, pmax+1
                        do p1 = 1, pmax+1
                            dx = aimag(eta1(m)*eta2(m)*alpha_tilde_2(p2)*conjg(alpha_tilde_1(p1)))
                            x_seg_mat(p2, j, p1, k) = x_seg_mat(p2, j, p1, k) + dx
                            x_seg_mat(p1, k, p2, j) = x_seg_mat(p1, k, p2, j) + dx
                        end do
                    end do
                end do
            end do
         end do

         return 
      end
