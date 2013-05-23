function val3 = insert_clique3(poss, val3, val2)

%val4 vector 1x16
%val2 vector 1x4
%poss = '34'
switch poss
    case '12'
        val3(1:2)   = val3(1:2)   + val2(1);
        val3(3:4)   = val3(3:4)   + val2(2);
        val3(5:6)   = val3(5:6)   + val2(3);
        val3(7:8)   = val3(7:8)   + val2(4);
    case '13'
        val3([1 3]) = val3([1 3]) + val2(1);
        val3([2 4]) = val3([2 4]) + val2(2);
        val3([5 7]) = val3([5 7]) + val2(3);
        val3([6 8]) = val3([6 8]) + val2(4);
    case '23'
        val3([1 5]) = val3([1 5]) + val2(1);
        val3([2 6]) = val3([2 6]) + val2(2);
        val3([3 7]) = val3([3 7]) + val2(3);
        val3([4 8]) = val3([4 8]) + val2(4);
end