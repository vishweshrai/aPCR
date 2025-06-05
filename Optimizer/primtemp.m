function [tm,h,s] = primtemp(seq,pcon, monv, biv, catrat, a, d, g)
    hsfildir = 'G:\Viral Genome Encapsulation Project Exp Data\aPCR_Trial\Unified_HS_1M_NaCl.xlsx';
    hs1MNaCl = readtable(hsfildir,VariableNamingRule="preserve");
    h = 0; s = 0;
    gc = (count(seq,[{'G'},{'C'}]))/strlength(seq);
    for i = 1:16
        pat = hs1MNaCl.Sequence(i);
        nn = length(strfind(seq, pat));
        h = h + nn*hs1MNaCl.("ΔH° kcal/mol")(i);
        s = s + 0.001*nn*hs1MNaCl.("ΔS° cal/K·mol")(i);
    end
    % 5' Terminal
    if (seq(1) == 'G')||(seq(1) == 'C')
        h = h + hs1MNaCl.("ΔH° kcal/mol")(17);
        s = s + 0.001*hs1MNaCl.("ΔS° cal/K·mol")(17);
    end
    if (seq(1) == 'A')||(seq(1) == 'T')
        h = h + hs1MNaCl.("ΔH° kcal/mol")(18);
        s = s + 0.001*hs1MNaCl.("ΔS° cal/K·mol")(18);
    end

    % 3' Terminal
    if (seq(end) == 'G')||(seq(end) == 'C')
        h = h + hs1MNaCl.("ΔH° kcal/mol")(17);
        s = s + 0.001*hs1MNaCl.("ΔS° cal/K·mol")(17);
    end
    if (seq(end) == 'A')||(seq(end) == 'T')
        h = h + hs1MNaCl.("ΔH° kcal/mol")(18);
        s = s + 0.001*hs1MNaCl.("ΔS° cal/K·mol")(18);
    end
    if seq == seqrcomplement(seq)
        s = s + 0.001*hs1MNaCl.("ΔS° cal/K·mol")(19);
    end
    R = 1.9872*10^(-3); % R in kcal/(mol.K) : 1.9872*1E-3
    tm = h/(s+R*log(pcon));
    if catrat < 0.22
        tm = ((1/tm) + ((4.29*gc- 3.95)*log(monv)+0.940*(log(monv))^2)*10^(-5))^(-1);
    elseif catrat >= 0.22
        tm = ((1/tm)+(a - 0.911*log(biv) + gc*(6.26+d*log(biv))+(1/(2*(strlength(seq)-1))*(-48.2+52.5*log(biv)+g*(log(biv))^2)))*10^(-5))^(-1);
    end
    tm = tm - 273.15;
end