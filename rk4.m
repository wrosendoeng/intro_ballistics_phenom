function solver = rk4(t,y)
        
for iter = 1:t_points
    k1 = accel_rk4(t,y);
    k2 = accel_rk4(t+0.5*dt,y+0.5*dt*k1);
    k3 = accel_rk4(t+0.5*dt,y+0.5*dt*k2);
    k4 = accel_rk4(t+1.0*dt,y+1.0*dt*k3);

    solver = y + dt/6*(k1+2*(k2+k3)+k4);
    t = t + dt;	
end
end