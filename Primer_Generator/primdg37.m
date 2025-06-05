function dG37 = primdg37(seq,monv,biv,tann,a,d,g,catrat)
    hsfildir = 'G:\Viral Genome Encapsulation Project Exp Data\aPCR_Trial\Unified_HS_1M_NaCl.xlsx';
    hs1MNaCl = readtable(hsfildir,VariableNamingRule="preserve");
    h = 0; s = 0;
    gc = (count(seq,[{'G'},{'C'}]))/strlength(seq);
    for i = 1:16
        pat = hs1MNaCl.Sequence(i);
        nn = length(strfind(seq, pat));
        h = h + nn*(hs1MNaCl.("ΔH° kcal/mol")(i));
        s = s + 0.001*nn*(hs1MNaCl.("ΔS° cal/K·mol")(i));
    end
    % 5' Terminal
    if (seq(1) == 'G')||(seq(1) == 'C')
        h = h + (hs1MNaCl.("ΔH° kcal/mol")(17));
        s = s + 0.001*(hs1MNaCl.("ΔS° cal/K·mol")(17));
    end
    if (seq(1) == 'A')||(seq(1) == 'T')
        h = h + (hs1MNaCl.("ΔH° kcal/mol")(18));
        s = s + 0.001*(hs1MNaCl.("ΔS° cal/K·mol")(18));
    end
    % 3' Terminal
    if (seq(end) == 'G')||(seq(end) == 'C')
        h = h + (hs1MNaCl.("ΔH° kcal/mol")(17));
        s = s + 0.001*(hs1MNaCl.("ΔS° cal/K·mol")(17));
    end
    if (seq(end) == 'A')||(seq(end) == 'T')
        h = h + (hs1MNaCl.("ΔH° kcal/mol")(18));
        s = s + 0.001*(hs1MNaCl.("ΔS° cal/K·mol")(18));
    end
    if seq == char(replace(reverse(string(seq)),["A","C","G","T"],["T","G","C","A"]))
        s = s + 0.001*(hs1MNaCl.("ΔS° cal/K·mol")(19));
    end
    s = 1000*s; h = 1000*h;
    if catrat < 0.22
        s = s + h*((4.29*gc- 3.95)*log(monv)+0.940*(log(monv))^2)*10^(-5);
    elseif catrat >= 0.22
        s = s+h*(a - 0.911*log(biv) + gc*(6.26+d*log(biv))+(1/(2*(strlength(seq)-1))*(-48.2+52.5*log(biv)+g*(log(biv))^2)))*10^(-5);
    end
    tann = tann + 273.15;
    dG37 = h-tann*s; %kcal/mol
end