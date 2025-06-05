%% Primer Generator for aPCR Amplification of Scaffold
clc
clear
clear cache

prompta = {'Enter plasmid file location:','Enter plasmid file name:'};
dlgtitle = 'Input Custom Genome Sequence File';
fieldsize = [1 50; 1 50];
definput = {'G:\Viral Genome Encapsulation Project Exp Data\aPCR_Trial\'...
    ,'D2_E_mut.txt'};
plas = inputdlg(prompta,dlgtitle,fieldsize,definput);
promptb = {'Is it double-stranded (2) or single stranded (1)?',['Enter ' ...
    'plasmid concentration (ng/uL):']};
dlgtitle = 'Input Plasmid Detail';
fieldsize = [1 50; 1 50];
definput = {'2','2'};
dna = inputdlg(promptb,dlgtitle,fieldsize,definput);
promptc = {'Enter the desired length of scaffold:'};
dlgtitle = 'Input Scaffold Length';
fieldsize = [1 50];
definput = {'3700'};
scaf = inputdlg(promptc,dlgtitle,fieldsize,definput);
promptd = {'Enter forward primer concentration (uM):',['Enter forward/' ...
    'reverse primer concentration ratio:'],['Enter monovalent cation' ...
    'concentration (mM):'],['Enter bivalent cation concentration ' ...
    '(mM):'],'Enter dNTPs concentration (mM):'};
dlgtitle = 'Input Reaction Conditions';
fieldsize = [1 50; 1 50; 1 50; 1 50;1 50];
definput = {'1','50','60','2','0.3'};
react = inputdlg(promptd,dlgtitle,fieldsize,definput);
prompte = {'Enter desired minimum length of primers (nt):',['Enter' ...
    ' desired maximum length of primers (nt):'],['Enter desired' ...
    ' minimum melting temperature of primers (oC):'],['Enter desired' ...
    ' maximum melting temperature of primers (oC):'],['Enter desired' ...
    ' minimum temperature difference of primers (oC):'],['Enter ' ...
    'desired maximum temperature difference of primers (oC):']};
dlgtitle = 'Input Primer Design Parameters';
fieldsize = [1 50; 1 50; 1 50; 1 50; 1 50; 1 50];
definput = {'18','22','60','65','2','3'};
descon = inputdlg(prompte,dlgtitle,fieldsize,definput);

%% Processing User Input
dir = plas(1); seqfile = plas(2);
strand = str2double(dna(1)); plaswcon = str2double(dna(2));
Ls = str2double(scaf(1));
monv =  str2double(react(3))*10^(-3);
biv =  str2double(react(4))*10^(-3);
fcon = str2double(react(1))*10^(-6);
dNTP =  str2double(react(5))*10^(-3);
lmin = str2double(descon(1)); lmax = str2double(descon(2));
Tmin = str2double(descon(3)); Tmax = str2double(descon(4));
dTmin = str2double(descon(5)); dTmax = str2double(descon(6));
plasmid = char(extractFileText([dir{1},seqfile{1}]));
if Ls >= strlength(plasmid)
    disp("Insufficient genome length! Choose a longer genome.")
    exit
end
an = count(plasmid,'A'); tn = count(plasmid,'T'); cn = count(plasmid,'C');
gn = count(plasmid,'G');
plasmolwt = (an*313.2)+(tn*304.2)+(cn*289.2)+(gn*329.2)+79.0;
plascon = plaswcon*10^(-3)/(strand*plasmolwt);
ifcon = fcon;
if fcon >= 6*plascon
    fcon = fcon - plascon/2;
end
rcon = ifcon/str2double(react(2));
ircon = rcon;
if rcon >= 6*plascon
    rcon = rcon - plascon/2;
else
    rcon = (rcon + plascon)/4;
end
tbiv = biv;
Ka = 3*10^4;
D = (Ka*dNTP-Ka*biv+1)^2 + 4*Ka*biv;
biv = (-(Ka*dNTP-Ka*biv + 1)+sqrt(D))/(2*Ka);
catrat = sqrt(biv)/monv;
if catrat > 6.0
    a = 3.92;
    d = 1.42;
    g = 8.31;
elseif (catrat >= 0.22)&&(catrat <= 6.0)
    a = 3.92*(0.843-0.352*sqrt(monv)*log(monv));
    d = 1.42*(1.279 - 4.03*10^(-3)*log(monv) - 8.03*10^(-3)*(log(monv))^2);
    g = 8.31*(0.486 - 0.258*log(monv) + 5.25*10^(-3)*(log(monv))^3);
else
    a = 0;
    d = 0;
    g = 0;
end
neighbour = 5;
%% Global Primer Search
T = table();
pat = [{'G'},{'C'}];
fprintf('Launching global primer search...\n');
fprintf('%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n', ...
    'F-Start','F-Length','A/T.3','%F-GC','dlocGC','F-Tm','F-Status', ...
    'R-Start','R-Length','G/C.3','%R-GC','R-Tm','deltaTm','R-Status');
parfor i = 1:(strlength(plasmid))
    % Forward Primer Design
    for j = lmin:lmax
        shiftplasmid = circshift(plasmid,-(i-1));
        primj = shiftplasmid(1:j);
        if ((primj(end) == 'A')||(primj(end) == 'T'))
            GC = count(primj,pat);
            gcf = GC/strlength(primj);
            if (gcf > 0.40)&&(gcf <= 0.60)
                locgc3 = count(primj((j-neighbour+1):j),pat)/neighbour;
                locgc5 = count(primj(1:neighbour),pat)/neighbour;
                if locgc3 < locgc5
                    [tmf,~,~] = primtemp(primj,fcon,monv,biv,catrat,a,d,g);
                    if (tmf >= Tmin)&&(tmf <= Tmax)
                        fprintf(['%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s' ...
                            '%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n'],i,j, ...
                            primj(end),gcf,(locgc5-locgc3),tmf, ...
                            'PASS','NA','NA','NA','NA','NA','NA','NA');
                        % Reverse Primer Design
                        if ((i + Ls)-1) > strlength(plasmid)
                            stop = mod(i + Ls,strlength(plasmid)) - 1;
                        else
                            stop = (i + Ls) - 1;
                        end
                        for k = lmin:lmax
                            shiftplasmid = circshift(plasmid,-stop);
                            primk = seqrcomplement(shiftplasmid(end-k+1:end));
                            if ((primk(end) == 'C')||(primk(end) == 'G'))
                                GC = count(primk,pat);
                                gcr = GC/strlength(primk);
                                if (gcr >= 0.40)&&(gcr <= 0.60)
                                    [tmr,~,~] = primtemp(primk,rcon, monv, biv, catrat, a, d, g);
                                    if (tmr >= Tmin)&&(tmr <= Tmax)
                                        if ((tmr-tmf) >= dTmin)&&((tmr-tmf) <= dTmax)
                                            tabrow = table({primj},i,tmf,j,gcf,{primk},stop,tmr,k,gcr,'VariableNames',{'F_Primer','F_Start','F_Tm','F_Length','F_GC','R_Primer','R_Start','R_Tm','R_Length','R_GC'});
                                            T = [T;tabrow];
                                            fprintf('%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s\n',i,j,primj(end),gcf,(locgc5-locgc3),tmf,'PASS',stop,k,primk(end),gcr,tmr,(tmr-tmf),'PASS');
                                        else
                                            fprintf('%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s\n',i,j,primj(end),gcf,(locgc5-locgc3),tmf,'PASS',stop,k,primk(end),gcr,tmr,(tmr-tmf),'FAIL');
                                        end
                                    else
                                        fprintf('%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s%-9d%-9d%-9s%-9.2f%-9.2f%-9s%-9s\n',i,j,primj(end),gcf,(locgc5-locgc3),tmf,'PASS',stop,k,primk(end),gcr,tmr,'NA','FAIL');
                                    end
                                else
                                    fprintf('%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s%-9d%-9d%-9s%-9.2f%-9s%-9s%-9s\n',i,j,primj(end),gcf,(locgc3-locgc5),tmf,'PASS',stop,k,primk(end),gcr,'NA','NA','FAIL');
                                end
                            else
                                fprintf('%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s%-9d%-9d%-9s%-9s%-9s%-9s%-9s\n',i,j,primj(end),gcf,(locgc3-locgc5),tmf,'PASS',stop,k,primk(end),'NA','NA','NA','FAIL');
                            end
                        end
                    else
                        fprintf('%-9d%-9d%-9s%-9.2f%-9.2f%-9.2f%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n',i,j,primj(end),gcf,(locgc3-locgc5),tmf,'FAIL','NA','NA','NA','NA','NA','NA','NA');
                    end
                else
                    fprintf('%-9d%-9d%-9s%-9.2f%-9.2f%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n',i,j,primj(end),gcf,(locgc3-locgc5),'NA','FAIL','NA','NA','NA','NA','NA','NA','NA');
                end
            else
                fprintf('%-9d%-9d%-9s%-9.2f%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n',i,j,primj(end),gcf,'NA','NA','FAIL','NA','NA','NA','NA','NA','NA','NA');
            end
        else
            fprintf('%-9d%-9d%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n',i,j,primj(end),'NA','NA','NA','FAIL','NA','NA','NA','NA','NA','NA','NA');
        end
    end
end
%% Homodimerization Possibility
q = numel(T.F_Primer);
qualified = false(q,1);
fprintf('\n\nRunning elimination protocol to knock-off primers with high homodimerization possibility.\n\n');
fprintf('%-10s%-10s%-24s%-10s%-10s%-10s%-24s%-10s\n','Primer','F-dG','F-Align-F','F-dG-F','F-Status','R-dG','R-Align-R','R-Status');
F_Primer = parallel.pool.Constant(T.F_Primer);
R_Primer = parallel.pool.Constant(T.R_Primer);
parfor i = 1:q
    % Forward Primer Homodimerization
    primer_f = F_Primer.Value{i,:};
    primer_r = R_Primer.Value{i,:};
    dgf = primdg37(primer_f);
    [~, fAlignf] = swalign(seqcomplement(primer_f),seqrcomplement(primer_f),'GapOpen',1000,'Alphabet','NT');
    dgff = dGmax(fAlignf);
    if dgff > 0.33*dgf
        % Reverse Primer Homodimerization
        dgr = primdg37(primer_r);
        [~, rAlignr] = swalign(seqcomplement(primer_r),seqrcomplement(primer_r),'GapOpen',1000,'Alphabet','NT');
        dgrr = dGmax(rAlignr);
        if dgrr > 0.33*dgr
            qualified(i,1) = true;
            fprintf('%-10d%-10.2f%-24s%-10.2f%-10s%-10.2f%-24s%-10.2f%-10s\n',i,dgf,fAlignf(2,:),dgff,'TRUE',dgr,rAlignr(2,:),dgrr,'TRUE');
        else
            fprintf('%-10d%-10.2f%-24s%-10.2f%-10s%-10.2f%-24s%-10.2f%-10s\n',i,dgf,fAlignf(2,:),dgff,'TRUE',dgr,rAlignr(2,:),dgrr,'FAIL');
        end
    else
        fprintf('%-10d%-10.2f%-24s%-10.2f%-10s%-10s%-24s%-10s%-10s\n',i,dgf,fAlignf(2,:),dgff,'FAIL','NA','NA','NA','NA');
    end
end
Q = T(qualified,:);
%% Homo-/dinucleotide Repeats
q = numel(Q.F_Primer);
skipper = true(q,1);
nucleotides = unique(nchoosek(repmat('ACGT', 1,4), 1), 'rows');
dinucleotides = parallel.pool.Constant(unique(nchoosek(repmat('ACGT', 1,4), 2), 'rows'));
fprintf('\nRunning elimination protocol to further knock-off primers with long homo/dinucleotide runs.\n\n');
fprintf('%-24s%-10s%-10s%-10s%-10s%-10s%-25s%-10s%-10s%-10s%-10s%-10s\n','F-Primer','F-NT','F-#NT','F-dNT','F-#dNT','F_Status','R-Primer','R-NT','R-#NT','R-dNT','R-#dNT','R_Status');
F_Primer = parallel.pool.Constant(Q.F_Primer);
R_Primer = parallel.pool.Constant(Q.R_Primer);
parfor i = 1:q
    % Homonucleotide Repeats
    primer_f = F_Primer.Value{i,:};
    primer_r = R_Primer.Value{i,:};
    j = 1;
    while (j <= length(nucleotides))&&(skipper(i,1) == true)
        pattern = [nucleotides(j) '+'];
        matches = regexp(primer_f, pattern, 'match');
        numocc = cellfun(@length, matches);
        if max(numocc) >= 4
            fprintf('%-24s%-10s%-10d%-10s%-10s%-10s%-25s%-10s%-10s%-10s%-10s%-10s\n',primer_f,nucleotides(j),max(numocc),'NA','NA','FAIL',primer_r,'NA','NA','NA','NA','NA');
            skipper(i,1) = false;
            break
        else
            j = j + 1;
        end
    end
    j = 1;
    while (j <= length(dinucleotides.Value(:,1)))&&(skipper(i,1) == true)
        pattern = ['(' dinucleotides.Value(j,:) ')' '+'];
        matches = regexp(primer_f, pattern, 'match');
        numocc = cellfun(@length, matches)./2;
        if max(numocc) >= 4
            fprintf('%-24s%-10s%-10s%-10s%-10d%-10s%-25s%-10s%-10s%-10s%-10s%-10s\n',primer_f,'NA','NA',dinucleotides.Value(j,:),max(numocc),'FAIL',primer_r,'NA','NA','NA','NA','NA');
            skipper(i,1) = false;
            break
        else
            j = j + 1;
        end
    end
    j = 1;
    while (j <= length(nucleotides))&&(skipper(i,1) == true)
        pattern = [nucleotides(j) '+'];
        matches = regexp(primer_r, pattern, 'match');
        numocc = cellfun(@length, matches);
        if max(numocc) >= 4
            fprintf('%-24s%-10s%-10s%-10s%-10s%-10s%-25s%-10s%-10d%-10s%-10s%-10s\n',primer_f,'NA','NA','NA','NA','NA',primer_r,nucleotides(j),max(numocc),'NA','NA','FAIL');
            skipper(i,1) = false;
            break
        else
            j = j + 1;
        end
    end
    j = 1;
    while (j <= length(dinucleotides.Value(:,1)))&&(skipper(i,1) == true)
        pattern = ['(' dinucleotides.Value(j,:) ')' '+'];
        matches = regexp(primer_r, pattern, 'match');
        numocc = cellfun(@length, matches)./2;
        if max(numocc) >= 4
            fprintf('%-24s%-10s%-10s%-10s%-10s%-10s%-24s%-10s%-10d%-10s%-10s%-10s\n',primer_f,'NA','NA','NA','NA','NA',primer_r,'NA','NA',dinucleotides.Value(j,:),max(numocc),'FAIL');
            skipper(i,1) = false;
            break
        else
            j = j + 1;
        end
    end
    if skipper(i,1) == true
        fprintf('%-24s%-10s%-10s%-10s%-10s%-10s%-25s%-10s%-10s%-10s%-10s%-10s\n',primer_f,'NA','NA','NA','NA','PASS',primer_r,'NA','NA','NA','NA','PASS');
    end
end
V = Q(skipper,:);
%% Heterodimerization Possiblity
q = numel(V.F_Primer);
eliminator = true(q,1);
fprintf('\nRunning elimination protocol to further knock-off primers with heterdimerization possibility.\n\n');
fprintf('%-25s%-25s%-24s%-10s%-10s%-10s%-10s\n','F-Primer','R-Primer','F-Align-R','F-dG','R-dG','F-dG-R','FR-Status');
F_Primer = parallel.pool.Constant(V.F_Primer);
R_Primer = parallel.pool.Constant(V.R_Primer);
parfor i = 1:q
    primer_f = F_Primer.Value{i,:};
    primer_r = R_Primer.Value{i,:};
    [~, fAlignr] = swalign(primer_f,seqrcomplement(primer_r),'GapOpen',1000,'Alphabet','NT');
    dgf = primdg37(primer_f);
    dgr = primdg37(primer_r);
    dgfr = dGmax(fAlignr);
    if (dgfr < 0.33*dgf)||(dgfr < 0.33*dgr)
        fprintf('%-25s%-25s%-24s%-10.2f%-10.2f%-10.2f%-10s\n',primer_f,primer_r,fAlignr(2,:),dgf,dgr,dgfr,'FAIL');
        eliminator(i,1) = false;
    else
        fprintf('%-25s%-25s%-24s%-10.2f%-10.2f%-10.2f%-10s\n',primer_f,primer_r,fAlignr(2,:),dgf,dgr,dgfr,'PASS');
    end
end
Y = V(eliminator,:);
%% Mispriming Possibility
q = numel(Y.F_Primer);
winner = true(q,1);
fprintf('\n\nFurther shortlisting primer-pairs based on mispriming possibility.\n');
fprintf('%-25s%-25s%-10s%-25s%-10s%-10s%-10s%-10s%-10s\n','Primer','p-Align-D','p-dG-G','p-Align-G','p-dG-D','Length','Last5','Stop','Status');
F_Primer = parallel.pool.Constant(Y.F_Primer);
F_Start = parallel.pool.Constant(Y.F_Start);
F_Length = parallel.pool.Constant(Y.F_Length);
R_Primer = parallel.pool.Constant(Y.R_Primer);
R_Start = parallel.pool.Constant(Y.R_Start);
R_Length = parallel.pool.Constant(Y.R_Length);
parfor i = 1:q
    primer_f = F_Primer.Value{i,:};
    primer_r = R_Primer.Value{i,:};
    dgf = primdg37(primer_f);
    dgr = primdg37(primer_r);
    % f:G~fTc> Forward primer alignment with leading strand to check for the
    % availability of alternate binding sites on the lagging strand.
    % If the forward primer binds to the lagging strand, at a location
    % different from the intended one, at that site, the 5'-3' sequence of
    % the primer will be complementary to the 3'-5' sequence of the lagging
    % strand. In short, the 5'-3' sequence of the forward primer will
    % align with the 5'-3' sequence of the leading strand.
    leadplas = [plasmid,plasmid(1:lmax-1)];
    lagplas = seqcomplement(leadplas);
    start_f = F_Start.Value(i);
    length_f = F_Length.Value(i);
    fleadplas = [leadplas(1:start_f-1),leadplas(start_f+length_f:end)];
    fDalstruct = localalign(fleadplas,primer_f,'numaln', 3,'GapOpen',5000,'Alphabet','NT');
    % f:D~fM> Forward primer alignment with lagging strand to check if there are
    % any binding sites on the leading strand. If the forward primer binds
    % to the leading strand, at that site, the 5'-3' sequence will be
    % complementary to the 3'-5' sequence of the forward primer. In short,
    % the 3'-5' sequence of the forward primer will align with the 3'-5'
    % sequence of the lagging strand.
    fGalstruct = localalign(lagplas,seqreverse(primer_f),'numaln', 3,'GapOpen',5000,'Alphabet','NT');
    % r:D~rTc> Reverse primer alignment with lagging strand to check for the
    % availability of alternate binding sites on the leading strand.
    % If the reverse primer binds to the leading strand, at that site,
    % the 3'-5' sequence of the reverse primer will be complementary to the
    % 5'-3' sequence of the leading strand. In short, the 3'-5' sequence of
    % the reverse primer will align with the 3'-5' sequence of the lagging
    % strand.
    start_R = R_Start.Value(i);
    length_R = R_Length.Value(i);
    rlagplas = [lagplas(1:start_R-length_R),lagplas(start_R+1:end)];
    rGalstruct = localalign(rlagplas,seqreverse(primer_r),'numaln', 3,'GapOpen',5000,'Alphabet','NT');
    % r:G~rM> Reverse primer alignment with the leading strand if there are any
    % binding sites on the lagging strand. If the reverse primer binds to
    % the lagging strand, at the site, the 5'-3' sequence of the reverse
    % primer will be complementary to the 3'-5' sequence of the of the
    % lagging strand. In short, the 5'-3' sequence of the reverse primer
    % will align with the 5'-3' sequence of leading strand.
    rDalstruct = localalign(leadplas,primer_r,'numaln', 3,'GapOpen',5000,'Alphabet','NT');
    fDalign = fDalstruct.Alignment;
    fGalign = fGalstruct.Alignment;
    for j = 1:numel(fDalign)
        dgfg = dGmax(fDalign{j});
        matches = regexp(strtrim(replace(fDalign{j}(2,:),'|','X')), ['X' '+'], 'match'); numocc = cellfun(@length, matches);
        if (dgfg < 0.75*dgf)&&(max(numocc) >= 9)            
            if (dgfg < 0.5*dgf)
                fprintf('%-25s%-25s%-10.2f%-25s%-10s%-10d%-10s%-10s%-10s\n',primer_f,fDalign{j}(2,:),dgfg,'NA','NA',max(numocc),'NA','NA','FAIL');
                winner(i,1) = false;
                break;
            elseif (contains(fDalign{j}(2,end-4:end),'') == false)&&(fDalstruct.Stop(j,2)>(strlength(primer_f)-5))
                fprintf('%-25s%-25s%-10.2f%-25s%-10s%-10d%-10s%-10d%-10s\n',primer_f,fDalign{j}(2,:),dgfg,'NA','NA',max(numocc),fDalign{j}(2,end-4:end),fDalstruct.Stop(j,2),'FAIL');
                winner(i,1) = false;
                break;
            else
                fprintf('%-25s%-25s%-10.2f%-25s%-10s%-10d%-10s%-10d%-10s\n',primer_f,fDalign{j}(2,:),dgfg,'NA','NA',max(numocc),fDalign{j}(2,end-4:end),fDalstruct.Stop(j,2),'PASS');
            end
        else
            fprintf('%-25s%-25s%-10.2f%-25s%-10s%-10s%-10s%-10d%-10s\n',primer_f,fDalign{j}(2,:),dgfg,'NA','NA','NA',fDalign{j}(2,end-4:end),fDalstruct.Stop(j,2),'PASS');
        end
    end
    if winner(i,1) == false
        continue;
    end
    for j = 1:numel(fGalign)
        dgfd = dGmax(fGalign{j});
        if (dgfd < 0.25*dgf)
            matches = regexp(strtrim(replace(fGalign{j}(2,:),'|','X')), ['X' '+'], 'match'); numocc = cellfun(@length, matches);
            if (dgfd < 0.5*dgf)||(max(numocc) >= 9)
                fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10d%-10s%-10s%-10s\n',primer_f,'NA','NA',fGalign{j}(2,:),dgfd,max(numocc),'NA','NA','FAIL');
                winner(i,1) = false;
                break;
            elseif (contains(fGalign{j}(2,end-4:end),'') == false)&&(fGalstruct.Stop(j,2)>(strlength(primer_f)-5))
                fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10d%-10s%-10d%-10s\n',primer_f,'NA','NA',fGalign{j}(2,:),dgfd,max(numocc),fGalign{j}(2,end-4:end),fGalstruct.Stop(j,2),'FAIL');
                winner(i,1) = false;
                break;
            else
                fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10d%-10s%-10d%-10s\n',primer_f,'NA','NA',fGalign{j}(2,:),dgfd,max(numocc),fGalign{j}(2,end-4:end),fGalstruct.Stop(j,2),'PASS');
            end
        else
            fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10s%-10s%-10s%-10s\n',primer_f,'NA','NA',fGalign{j}(2,:),dgfd,'NA','NA','NA','PASS');
        end
    end
    if winner(i,1) == false
        continue;
    end
    rGalign = rGalstruct.Alignment;
    rDalign = rDalstruct.Alignment;
    for j = 1:numel(rGalign)
        dgrd = dGmax(rGalign{j});
        if (dgrd < 0.25*dgr)
            matches = regexp(strtrim(replace(rGalign{j}(2,:),'|','X')), ['X' '+'], 'match'); numocc = cellfun(@length, matches);
            if (dgrd < 0.50*dgr)||(max(numocc) >= 9)
                fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10d%-10s%-10s%-10s\n',primer_r,'NA','NA',rGalign{j}(2,:),dgrd,max(numocc),'NA','NA','FAIL');
                winner(i,1) = false;
                break;
            elseif (contains(rGalign{j}(2,end-4:end),'') == false)&&(rGalstruct.Stop(j,2)>(strlength(primer_r)-5))
                fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10d%-10s%-10d%-10s\n',primer_r,'NA','NA',rGalign{j}(2,:),dgrd,max(numocc),rGalign{j}(2,end-4:end),rGalstruct.Stop(j,2),'FAIL');
                winner(i,1) = false;
                break;
            else
                fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10d%-10s%-10d%-10s\n',primer_r,'NA','NA',rGalign{j}(2,:),dgrd,max(numocc),rGalign{j}(2,end-4:end),rGalstruct.Stop(j,2),'PASS');
            end
        else
            fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10s%-10s%-10s%-10s\n',primer_r,'NA','NA',rGalign{j}(2,:),dgrd,'NA','NA','NA','PASS');
        end
    end
    if winner(i,1) == false
        continue;
    end
    for j = 1:numel(rDalign)
        dgrg = dGmax(rDalign{j});
        if (dgrg < 0.25*dgr)
            matches = regexp(strtrim(replace(rDalign{j}(2,:),'|','X')), ['X' '+'], 'match'); numocc = cellfun(@length, matches);
            if (dgrg < 0.50*dgr)||(max(numocc) >= 9)
                fprintf('%-25s%-25s%-10.2f%-25s%-10s%-10d%-10s%-10s%-10s\n',primer_r,rDalign{j}(2,:),dgrg,'NA','NA',max(numocc),'NA','NA','FAIL');
                winner(i,1) = false;
                break;
            elseif (contains(rDalign{j}(2,end-4:end),'') == false)&&(rDalstruct.Stop(j,2)>(strlength(primer_r)-5))
                fprintf('%-25s%-25s%-10s%-25s%-10.2f%-10d%-10s%-10d%-10s\n',primer_r,'NA','NA',rGalign{j}(2,:),dgrg,max(numocc),rGalign{j}(2,end-4:end),rGalstruct.Stop(j,2),'FAIL');
                winner(i,1) = false;
                break;
            else
                fprintf('%-25s%-25s%-10.2f%-25s%-10s%-10d%-10s%-10d%-10s\n',primer_r,rDalign{j}(2,:),dgrg,'NA','NA',max(numocc),rDalign{j}(2,end-4:end),rDalstruct.Stop(j,2),'PASS');
            end
        else
            fprintf('%-25s%-25s%-10.2f%-25s%-10s%-10s%-10s%-10d%-10s\n',primer_r,rDalign{j}(2,:),dgrg,'NA','NA','NA',rDalign{j}(2,end-4:end),rDalstruct.Stop(j,2),'PASS');
        end
    end
end
W = Y(winner,:);
%% Hairpinning/Looping Possibility
q = numel(W.F_Primer);
champion = true(q,1);
F_Primer = parallel.pool.Constant(W.F_Primer);
F_Length = parallel.pool.Constant(W.F_Length);
R_Primer = parallel.pool.Constant(W.R_Primer);
R_Length = parallel.pool.Constant(W.R_Length);
loopmin = 1; stemmin = 5;
fprintf('\n\nFurther shortlisting primer-pairs based on mispriming possibility.\n\n');
fprintf('%-10s%-11s%-10s%-12s%-12s%-15s%-11s%-10s\n','Primer','Stem-Start','Stem-Size','Stem5-3','Stem3-5','Loop','Stem-Align','Status');
parfor i = 1:q
    primer_f = F_Primer.Value{i,:};
    primer_r = R_Primer.Value{i,:};
    for j = 1:F_Length.Value(i) % Stem starts here
        for k = stemmin:floor(0.5*(F_Length.Value(i)-j)) % Length of stem
            stem53 = primer_f(j:j+k-1); % Sequence of stem
            for l = loopmin:(F_Length.Value(i)-j-2*k)
                stem35 = primer_f(j+k+l:j+2*k+l-1); % Sequence of stem
                loop = primer_f(j+k:j+k+l-1);
                [~,stemalign] = swalign(stem53,seqrcomplement(stem35),'GapOpen',100,'Alphabet','NT');
                warning off
                if isempty(stemalign) == false
                    stempat = [{'|||||'},{'|||| ||'},{'|| ||||'},{'||| |||'}];
                    alnstem = contains(stemalign(2,:),stempat);
                    if alnstem == true
                        fprintf('%-10d%-11d%-10d%-12s%-12s%-15s%-11s%-10s\n',i,j,k,stem53,stem35,loop,stemalign(2,:),'FAIL');
                        champion(i,1) = false;
                        break
                    else
                        fprintf('%-10d%-11d%-10d%-12s%-12s%-15s%-11s%-10s\n',i,j,k,stem53,stem35,loop,stemalign(2,:),'PASS');
                    end
                else
                    fprintf('%-10d%-11d%-10d%-12s%-12s%-15s%-11s%-10s\n',i,j,k,stem53,stem35,loop,stemalign(2,:),'PASS');
                end
            end
        end
    end
    for j = 1:R_Length.Value(i) % Stem starts here
        for k = stemmin:floor(0.5*(R_Length.Value(i)-j)) % Length of stem
            stem53 = primer_f(j:j+k-1); % Sequence of stem
            for l = loopmin:(R_Length.Value(i)-j-2*k)
                stem35 = primer_r(j+k+l:j+2*k+l-1); % Sequence of stem
                loop = primer_r(j+k:j+k+l-1);
                [~,stemalign] = swalign(stem53,seqrcomplement(stem35),'GapOpen',100,'Alphabet','NT');
                if isempty(stemalign) == false
                    stempat = [{'|||||'},{'|||| ||'},{'|| ||||'},{'||| |||'}];
                    alnstem = contains(stemalign(2,:),stempat);
                    if alnstem == true
                        fprintf('%-10d%-11d%-10d%-12s%-12s%-15s%-11s%-10s\n',i,j,k,stem53,stem35,loop,stemalign(2,:),'FAIL');
                        champion(i,1) = false;
                        break
                    else
                        fprintf('%-10d%-11d%-10d%-12s%-12s%-15s%-11s%-10s\n',i,j,k,stem53,stem35,loop,stemalign(2,:),'PASS');
                    end
                else
                    fprintf('%-10d%-11d%-10d%-12s%-12s%-15s%-11s%-10s\n',i,j,k,stem53,stem35,loop,stemalign(2,:),'PASS');
                end
            end
        end
    end
end
Z = W(champion,:);
%% Save Results
res = [dir{1},erase(seqfile{1},'.txt'),'_primer\',char(datetime('now','Format','yyyyMMMdd_HHmm')),'\'];
if isfolder(res) ~= 1
    mkdir(res);
end
writetable(Z,[res,erase(seqfile{1},'.txt'),'_',num2str(Ls),'_nt_scaffprimer.xlsx']);
cond = ["Plasmid Name";"Scaffold Length (nt)";"Monovalent Cations (mM)";"Bivalent Cations (mM)";"dNTPs (mM)";"Forward Primer (uM)";"Reverse Primer (uM)";"Tm,min (oC)";"Tm,max (oC)";"dTm,min (oC)";"dTm,max (oC)";"Lmin (nt)";"Lmax (nt)"];
val = [erase(seqfile,".txt");Ls;monv*10^3;tbiv*10^3;dNTP*10^3;ifcon*10^6;ircon*10^6;Tmin;Tmax;dTmin;dTmax;lmin;lmax];
C = table(cond, val,'VariableNames',{'Condition','Design Values'});
range1 = strcat(num2xlcol(numel(Z.Properties.VariableNames)+2),"1");
range2 = strcat(num2xlcol(numel(Z.Properties.VariableNames)+numel(C.Properties.VariableNames)+2),num2str(numel(C.Condition)));
range = char(strcat(range1,":",range2));
writetable(C,[res,erase(seqfile{1},'.txt'),'_',num2str(Ls),'_nt_scaffprimer.xlsx'],"Range",range);
