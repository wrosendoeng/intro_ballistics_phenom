function secant(f,x0,x1)
IMPLICIT NONE 
REAL(8) :: f, x0, x1, x2, secant

DO WHILE (i .LE. imax) ! 100 iterations to prevent infinite loop
    IF (ABS(f(x0)) .LT. eps1) THEN
        PRINT '(A32,E25.17)', "The value of x: ", x0/1.e-3
        PRINT '(A32,E25.17)', "The value of f(x): ", f(x0)/1.e-3
        PRINT *, "How many iterations were used in secant method: ", i
        STOP
    ELSE IF (ABS(f(x1)) .LT. eps1 .OR. ABS(x1-x0) .LT. eps2) THEN
        PRINT '(A32,E25.17)', "The value of x: ", x1/1.e-3
        PRINT '(A32,E25.17)', "The value of f(x): ", f(x1)/1.e-3
        PRINT *, "How many iterations were used in secant method: ", i
        STOP
    ELSE
        x2 = x1 - f(x1)/(f(x1)-f(x0))*(x1-x0)
        IF (ABS(f(x2)) .LT. eps1 .OR. ABS(x2-x1) .LT. eps2) THEN
            PRINT '(A32,E25.17)', "The value of x: ", x2/1.e-3
            PRINT '(A32,E25.17)', "The value of f(x): ", f(x2)/1.e-3
            PRINT *, "How many iterations were used in Secant method: ", i
            STOP
        ELSE
            x0 = x1
            x1 = x2
        END IF 
    END IF
    i = i+1
END DO

end