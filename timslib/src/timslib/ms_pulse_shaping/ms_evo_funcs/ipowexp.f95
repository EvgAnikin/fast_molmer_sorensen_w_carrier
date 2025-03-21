      subroutine ipowexp_arr(val, p, e, t)
        use alpha_chi_expr
        implicit none
        integer p
        double precision e, t
        double complex val(p+1)

!F2PY intent(out) :: val

        val = ipowexp_a(p, e, t)

      end

      subroutine ipowexp2_arr(val, p1, p2, e1, e2, t)
        use alpha_chi_expr
        implicit none
        integer p1, p2
        double precision e1, e2, t
        double complex val(p1+1, p2+1)

!F2PY intent(out) :: val

        val = ipowexp2_a(p1, p2, e1, e2, t)

      end
