function P = P_matrix(q)

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    P = [
            q0,     -q1,        -q2,        -q3;
            q1,      q0,        -q3,         q2;
            q2,      q3,         q0,        -q1;
            q3,     -q2,         q1,         q0
    ];

end