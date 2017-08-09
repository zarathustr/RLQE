% The implementation of RLQE with no robust law
% 
% author: Jin Wu, Zebo Zhou et al.
% e-mail: jin_wu_uestc@hotmail.com; klinsmann.zhou@gmail.com


function [ qy, W, Sigma_qy] = RLQE(Db, Dr, qk_1, Sigma_Db, Sigma_q_k_1, j)

    W = [ Db(1) * Dr(1) + Db(2) * Dr(2) + Db(3) * Dr(3),                 -Db(3) * Dr(2) + Db(2) * Dr(3),                  Db(3) * Dr(1) - Db(1) * Dr(3),                 -Db(2) * Dr(1) + Db(1) * Dr(2);
                         -Db(3) * Dr(2) + Db(2) * Dr(3),  Db(1) * Dr(1) - Db(2) * Dr(2) - Db(3) * Dr(3),                  Db(2) * Dr(1) + Db(1) * Dr(2),                  Db(3) * Dr(1) + Db(1) * Dr(3);
                          Db(3) * Dr(1) - Db(1) * Dr(3),                  Db(2) * Dr(1) + Db(1) * Dr(2), -Db(1) * Dr(1) + Db(2) * Dr(2) - Db(3) * Dr(3),                  Db(3) * Dr(2) + Db(2) * Dr(3);
                         -Db(2) * Dr(1) + Db(1) * Dr(2),                  Db(3) * Dr(1) + Db(1) * Dr(3),                  Db(3) * Dr(2) + Db(2) * Dr(3), -Db(1) * Dr(1) - Db(2) * Dr(2) + Db(3) * Dr(3)];


    qy = 1.0 / j * (W + (j - 1) * eye(4)) * qk_1;

    K = K_matrix(qk_1, Dr);
    
    Sigma_qy = 1.0 / j^2 * (K' * Sigma_Db * K +  ...
                           (W + (j - 1) * eye(4)) * Sigma_q_k_1 * (W + (j - 1) * eye(4))) / norm(qy)^2;

    qy = qy ./ norm(qy);

end

