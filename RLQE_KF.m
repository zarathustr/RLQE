function [qk, Sigma_qk, qy] = RLQE_KF(Db, Dr, Sigma_Db, ...
                                      omega, Sigma_gyro, ...
                                      qk_1, Sigma_qk_1, qy, Sigma_qy, dt, j)
                              
    wx = omega(1);        wy = omega(2);        wz = omega(3);

    len = length(Db(1, :));
    
    for i = 1 : len
        b = Db(:, i);
        r = Dr(:, i);
        
        b = b ./ norm(b);
        r = r ./ norm(r);
        
        [qy, W, Sigma_qy] = RLQE(b, r, qy, Sigma_Db, Sigma_qy, j);
    end
    
    omega4 = [     0 , -wx , -wy , -wz ; ...
                  wx ,   0 ,  wz , -wy ; ...
                  wy , -wz ,   0 ,  wx ; ...
                  wz ,  wy , -wx ,   0];
    
    Ak = eye(4) + 0.5 .* omega4 * dt;
    
    q = qk_1;
    
    Ok = [     q(2),     q(3),   q(4);
              -q(1),    -q(4),  -q(3);
               q(3),    -q(1),  -q(2);
              -q(3),     q(2),  -q(1)];
    Qk = dt * dt * 0.25 * Ok * Sigma_gyro * Ok';
    
    qk_ = Ak * qk_1;
    Sigma_qk_ = Ak * Sigma_qk_1 * Ak' + Qk;
    Gk = Sigma_qk_ * inv(Sigma_qk_ + Sigma_qy);
    qk = qk_ + Gk * (qy - qk_);
    Sigma_qk = (eye(4) - Gk) * Sigma_qk_ * (eye(4) - Gk)' + Gk * Sigma_qy * Gk';
      
    qk = qk ./ norm(qk);
    
end
