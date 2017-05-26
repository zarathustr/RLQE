% The implementation of RLQE with no robust law
% 
% author: Jin Wu, Zebo Zhou et al.
% e-mail: jin_wu_uestc@hotmail.com; klinsmann.zhou@gmail.com


function [ qy, W,  Sigma_qy] = RLQE( Db,Dr, qk_1, Sigma_Db,  Sigma_q_k_1)

    W = [Db(1)*Dr(1) + Db(2)*Dr(2) + Db(3)*Dr(3), -Db(3)*Dr(2) + Db(2)*Dr(3), Db(3)*Dr(1) - Db(1)*Dr(3), -Db(2)*Dr(1) + Db(1)*Dr(2);
        -Db(3)*Dr(2) + Db(2)*Dr(3), Db(1)*Dr(1) - Db(2)*Dr(2) - Db(3)*Dr(3), Db(2)*Dr(1) + Db(1)*Dr(2), Db(3)*Dr(1) + Db(1)*Dr(3);
         Db(3)*Dr(1) - Db(1)*Dr(3), Db(2)*Dr(1) + Db(1)*Dr(2), -Db(1)*Dr(1) + Db(2)*Dr(2) - Db(3)*Dr(3), Db(3)*Dr(2) + Db(2)*Dr(3);
        -Db(2)*Dr(1) + Db(1)*Dr(2), Db(3)*Dr(1) + Db(1)*Dr(3), Db(3)*Dr(2) + Db(2)*Dr(3), -Db(1)*Dr(1) - Db(2)*Dr(2) + Db(3)*Dr(3)];


    qy=0.5*(W+eye(4))*qk_1;

    K=K_matrix(qk_1,Dr);

    Sigma_qy=0.25*(K'*Sigma_Db*K + Sigma_q_k_1 + 2*Sigma_q_k_1*W);

    qy=qy./norm(qy);

end

