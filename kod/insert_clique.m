function val4 = insert_clique(poss, val4, val2)

%val4 vector 1x16
%val2 vector 1x4
%poss = '34'

if length(poss) == 2
    for i = 0:3
        int1 = dec2bin(i,2);
        switch poss
            case '12'
                i1 = [[int1 '00' ]; [int1 '01'] ;[int1 '10'];[int1 '11']];
            case '13'
                i1 = [[int1(1) '0' int1(2) '0']; [int1(1) '0' int1(2) '1'] ;[int1(1) '1' int1(2) '0'];[int1(1) '1' int1(2) '1']];
            case '14'
                i1 = [[int1(1) '00' int1(2)]; [int1(1) '01' int1(2)] ;[int1(1) '10' int1(2)];[int1(1) '11' int1(2)]];
            case '23'
                i1 = [['0' int1 '0']; ['0' int1 '1'] ;['1' int1 '0']; ['1' int1 '1']];
            case '24'
                i1 = [['0' int1(1) '0'  int1(2)]; ['0' int1(1)  '1'  int1(2)] ;['1' int1(1) '0' int1(2) ];['1' int1(1) '1' int1(2) ]];
            case '34'
                i1 = [['00' int1]; ['01' int1] ;['10' int1];['11' int1]];
        end
        for j = 0:15
            i2 = dec2bin(j,4);
            for k = 1:4
                if i1(k,:) == i2
                    val4(j+1) = val4(j+1) + val2(i+1);
                end
            end
        end
    end
end
if length(poss) == 3
    for i = 0:7
        int1 = dec2bin(i,3);
        switch poss
            case '123'
                i1 = [[int1 '0' ]; [int1 '1'] ];
            case '124'
                i1 = [[int1(1:2) '0' int1(3)]; [int1(1:2) '1' int1(3)]];
            case '134'
                i1 = [[int1(1) '0' int1(2:3) ]; [int1(1) '1' int1(2:3) ] ];
            case '234'
                i1 = [['0' int1]; ['1' int1] ];
        end
        for j = 0:15
            i2 = dec2bin(j,4);
            for k = 1:2
                if i1(k,:) == i2
                    [i1(k,:);i2];
                    val4(j+1) = val4(j+1) + val2(i+1);
                end
            end
        end
    end
end

