function newton_raphson(f,nf,x0,x1)
IMPLICIT NONE 
REAL(8) :: f, nf, x0, x1, newton_raphson

DO WHILE (i .LE. imax) ! 100 iterations to prevent infinite loop
    IF (ABS(f(x0)) .LT. eps1) THEN
        PRINT '(A32,E25.17)', "The value of x: ", x0/1.e-3
        PRINT '(A32,E25.17)', "The value of f(x): ", f(x0)/1.e-3
        PRINT *, "How many iterations were used in Secant method: ", i
        STOP
    ELSE
        x1 = x0 - f(x0)/nf(f,x0)
        IF (ABS(f(x1)) .LT. eps1 .OR. ABS(x1-x0) .LT. eps2) THEN
            PRINT '(A32,E25.17)', "The value of x: ", x1/1.e-3
            PRINT '(A32,E25.17)', "The value of f(x): ", f(x1)/1.e-3
            PRINT *, "How many iterations were used in Secant method: ", i
            STOP
        ELSE
            x0 = x1
        END IF 
    END IF
    i = i+1
END DO

end