function col = num2xlcol(num)
%% Convert any number to Excel column name

    div = num;
    base = 26;
    let = char(65:1:90);
    col = "";
    while div > 0
        rem = mod(div,base);
        div = floor(div/base);
        if rem == 0
            rem = 26;
            div = div -1;
        end
        col = strcat(string(let(rem)), col);
    end
end