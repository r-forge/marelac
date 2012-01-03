       INTEGER FUNCTION Checkit(val, old, ipos)
       DOUBLE PRECISION val, old(*)
       INTEGER I

       Checkit = 0
       if (ipos .eq. 0) return 
       DO I = 1, ipos
         if (val .eq. old(I)) THEN
         Checkit = I
         return
         endif
       ENDDO

       END FUNCTION


       SUBROUTINE crosstab(input, ninput, colnr, rownr, valnr,                   &
     &                    cols, rows, nr, nc, cross)

       IMPLICIT NONE
       INTEGER colnr, rownr, valnr, nr, nc, ninput
       DOUBLE PRECISION input(3,*), rows(nr), cols(nc),cross(nr,nc)

       INTEGER I, J , icol, irow, nrow, ncol 
       DOUBLE PRECISION lrow, lcol
       INTEGER Checkit
  
         DO I = 1, ninput
           lrow = input(rownr,i)
           nrow = Checkit(lrow, rows, nr)
       !    IF (nrow .eq. 0) THEN
       !       call rwarn(" do not find row value")
       !	  ENDIF

           lcol = input(colnr,i)
           ncol = Checkit(lcol, cols, nc)
       !    IF (ncol .eq. 0) THEN
       !       call rwarn(" do not find column value")
       !	  ENDIF

           cross(nrow, ncol) = input(valnr,i)
         ENDDO


       END SUBROUTINE crosstab
