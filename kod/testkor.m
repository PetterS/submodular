E_1234 = [1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016]'
E_1234c = [1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016]'
E_123 = [101, 102, 103, 104, 105, 106, 107,108]'
E_124 = [201, 202, 203, 204, 205, 206, 207, 208]'
E_134 = [301, 302, 303, 304, 305, 306, 307, 308]'
E_234 = [401, 402, 403, 404, 405, 406, 407, 408]'

E_12 = [11, 12, 13, 14]'
E_13 = [21, 22, 23, 24]'
E_14 = [31, 32, 33, 34]'
E_23 = [41, 42, 43, 44]'
E_24 = [51, 52, 53, 54]'
E_34 = [61, 62, 63, 64]'


E_1234 = insert_clique('12', E_1234,E_12)
pause


E_1234 = insert_clique('123', E_1234,E_123);
E_1234 = insert_clique('124', E_1234,E_124);
E_1234 = insert_clique('134', E_1234,E_134);
E_1234 = insert_clique('234', E_1234,E_234);

E_1234 = insert_clique('12', E_1234,E_12);
E_1234 = insert_clique('13', E_1234,E_13);
E_1234 = insert_clique('14', E_1234,E_14);
E_1234 = insert_clique('23', E_1234,E_23);
E_1234 = insert_clique('24', E_1234,E_24);
E_1234 = insert_clique('34', E_1234,E_34);

E_1234c = insert_clique('123', E_1234c,E_123)
E_1234c = insert_clique('124', E_1234c,E_124)
E_1234c = insert_clique('134', E_1234c,E_134)
E_1234c = insert_clique('234', E_1234c,E_234)

E_1234c = insert_clique('12', E_1234c,E_12)
E_1234c = insert_clique('13', E_1234c,E_13)
E_1234c = insert_clique('14', E_1234c,E_14)
E_1234c = insert_clique('23', E_1234c,E_23)
E_1234c = insert_clique('24', E_1234c,E_24)
E_1234c = insert_clique('34', E_1234c,E_34)

E_1234c - E_1234

