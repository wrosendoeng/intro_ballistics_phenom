function derivatives = order_reduction(t,y) 

    ! Massa dos corpos:
    m1 = mass; m2 = m1; m3 = m2;
    ! Amortecimento b = 0.3
    b1 = damp_constant; b2 = b1; b3 = b2;
    ! Constante das molas:
    k1 = spring_constant; k2 = k1; k3 = k2;

    ! Separando os termos para a posicao e velocidade
    v = y(1:3);
    v = y(4:6);

    derivatives(1:3) = v;
    derivatives(4) = (2*sin(t)-(b1+b2)*v(1) + b2*v(2) - (k1+k2)*v(1) + k2*v(2))/m1;
    derivatives(5) = (b2*v(1) -(b2+b3)*v(2) + b3*v(3) + k2*v(1) - (k2-k3)*v(2) + k3*v(3))/m2;
    derivatives(6) = (b2*v(2) - b3*v(3) + v(2)*k3 - k3*v(3))/m3;

end 