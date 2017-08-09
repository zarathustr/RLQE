function W = W_matrix(Db, Dr)

W = [ Db(1) * Dr(1) + Db(2) * Dr(2) + Db(3)*Dr(3),                -Db(3) * Dr(2) + Db(2) * Dr(3),                  Db(3) * Dr(1) - Db(1) * Dr(3),                 -Db(2) * Dr(1) + Db(1) * Dr(2);
                   -Db(3) * Dr(2) + Db(2) * Dr(3), Db(1) * Dr(1) - Db(2) * Dr(2) - Db(3) * Dr(3),                  Db(2) * Dr(1) + Db(1) * Dr(2),                  Db(3) * Dr(1) + Db(1) * Dr(3);
                    Db(3) * Dr(1) - Db(1) * Dr(3),                 Db(2) * Dr(1) + Db(1) * Dr(2), -Db(1) * Dr(1) + Db(2) * Dr(2) - Db(3) * Dr(3),                  Db(3) * Dr(2) + Db(2) * Dr(3);
                   -Db(2) * Dr(1) + Db(1) * Dr(2),                 Db(3) * Dr(1) + Db(1) * Dr(3),                  Db(3) * Dr(2) + Db(2) * Dr(3), -Db(1) * Dr(1) - Db(2) * Dr(2) + Db(3) * Dr(3)];


end