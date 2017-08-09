function K = K_matrix(q, Dr)

    Drx = Dr(1);
    Dry = Dr(2);
    Drz = Dr(3);

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    P1 = [  q0,  q1, -q2, -q3;
           -q3,  q2,  q1, -q0;
            q2,  q3,  q0,  q1 ];
    P2 = [  q3,  q2,  q1,  q0;
            q0, -q1,  q2, -q3;
           -q1, -q0,  q3,  q2 ];
    P3 = [ -q2,  q3, -q0,  q1;
            q1,  q0,  q3,  q2;
            q0, -q1, -q2,  q3 ];

    K = Drx * P1 + Dry * P2 + Drz * P3;

end