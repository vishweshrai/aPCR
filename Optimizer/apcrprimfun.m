function [primer_f, primer_r] = apcrprimfun(plas,react,descon,plaswcon,fcon,ftor,Tmin,dTm,Tann)
dir = plas(1); seqfile = plas(2);
strand = str2double(plas(3));
Ls = str2double(plas(4));
monv =  str2double(react(1))*10^(-3);
biv =  str2double(react(2))*10^(-3);
dNTP =  str2double(react(3))*10^(-3);
lmin = str2double(descon(1)); lmax = str2double(descon(2));
plasmid = char(extractFileText([dir{1},seqfile{1}]));
if Ls >= strlength(plasmid)
    disp("Insufficient genome length! Choose a longer genome.")
    exit
end
an = count(plasmid,'A'); tn = count(plasmid,'T'); cn = count(plasmid,'C');
gn = count(plasmid,'G');
plasmolwt = (an*313.2)+(tn*304.2)+(cn*289.2)+(gn*329.2)+79.0;
plascon = plaswcon*10^(-3)/(strand*plasmolwt);
if fcon >= 6*plascon
    fcon = fcon - plascon/2;
end
rcon = fcon/ftor;
if rcon >= 6*plascon
    rcon = rcon - plascon/2;
else
    rcon = (rcon + plascon)/4;
end
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
neighbour = 5; % Definition of neighbourhood for calculating local %GC
% missf = 0.02; % Acceptable mispriming fraction
Tmrange = 0.25;
%% Global Primer Search
pat = [{'A'},{'T'}];
for i = 1:(strlength(plasmid))
    RT = 8.314*1E-3*(273.15+Tann)/4.182;
    % Forward Primer Design
    for j = lmin:lmax
        nucleotides = unique(nchoosek(repmat('ACGT', 1,4), 1), 'rows');
        dinucleotides = unique(nchoosek(repmat('ACGT', 1,4), 2), 'rows');
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
                    if (tmf >= Tmin)&&(tmf <= (Tmin + Tmrange))
                        %% Checking homodimerization possiblity in the forward primer
                        [~, fAlignf] = swalign(seqcomplement(primj),seqrcomplement(primj),'GapOpen',1000,'Alphabet','NT');
                        dgf = primdg(primj,monv,biv,tann,a,d,g,catrat);
                        dgffmin = -RT*log(missf*(plascon/fcon)*exp(-dgf/RT));
                        dgff = dGmax(fAlignf,monv,biv,tann,a,d,g,catrat);
                        if dgff > dgffmin
                            %% Checking homo/multinucleotide repeats in the forward primer
                            skipper = true;
                            for k = 1:length(nucleotides)
                                pattern = [nucleotides(k) '+'];
                                matches = regexp(primj, pattern, 'match');
                                numocc = cellfun(@length, matches);
                                if max(numocc) >= 4
                                    skipper = false;
                                    break
                                end
                            end
                            if skipper == true
                                for k = 1:length(dinucleotides.Value(:,1))
                                    pattern = ['(' dinucleotides.Value(k,:) ')' '+'];
                                    matches = regexp(primj, pattern, 'match');
                                    numocc = cellfun(@length, matches)./2;
                                    if max(numocc) >= 4
                                        skipper = false;
                                        break;
                                    end
                                end
                            end
                            if skipper == true
                                %% Checking mispriming possibility by the forward primer
                                winner = true;
                                leadplas = circshift(plasmid,-(i-1));
                                fleadplas = [repmat('X', 1, j),leadplas(j+1:end)];
                                fDalstruct = localalign(fleadplas,primj,'numaln', 3,'GapOpen',1000,'Alphabet','NT');
                                fGalstruct = localalign([repmat('X',1,j),seqcomplement(leadplas)],seqreverse(primj),'numaln', 3,'GapOpen',5000,'Alphabet','NT');
                                fDalign = fDalstruct.Alignment;
                                fGalign = fGalstruct.Alignment;
                                dgfgmin = -RT*log(missf*exp(-dgf/RT)/(1+(fcon-plascon)*exp(-dgf/RT)));
                                for l = 1:numel(fDalign)
                                    dgfg = dGmax(fDalign{l});
                                    matches = regexp(strtrim(replace(fDalign{l}(2,:),'|','X')), ['X' '+'], 'match'); numocc = cellfun(@length, matches);
                                    if (dgfg < dgfgmin)&&(max(numocc) >= 7)&&(fDalstruct.Stop(l,2)>=(strlength(primj)-3))
                                        winner = false;
                                        break;
                                    end
                                end
                                if winner == true
                                    dgfdmin = -RT*log(missf*exp(-dgf/RT)/(1+(fcon-plascon)*exp(-dgf/RT)));
                                    matches = regexp(strtrim(replace(fGalign{l}(2,:),'|','X')), ['X' '+'], 'match'); numocc = cellfun(@length, matches);
                                    for l = 1:numel(fGalign)
                                        dgfd = dGmax(fGalign{l});
                                        if (dgfd < dgfdmin)&&(max(numocc) >= 7)&&(fGalstruct.Stop(l,2)>=(strlength(primj)-3))
                                            winner = false;
                                            break;
                                        end
                                    end
                                end
                                if winner == true
                                    champion = true;
                                    loopmin = 1; stemmin = 5;
                                    for m = 1:j % Stem starts here
                                        for n = stemmin:floor(0.5*(j-m)) % Length of stem
                                            stem53 = primj(m:m+n-1); % Sequence of stem
                                            for o = loopmin:(j-m-2*n)
                                                stem35 = primj(m+n+o:m+2*n+o-1); % Sequence of stem
                                                [~,stemalign] = swalign(stem53,seqrcomplement(stem35),'GapOpen',100,'Alphabet','NT');
                                                warning off
                                                if isempty(stemalign) == false
                                                    stempat = [{'||||||'},{'|||| |||'},{'||| ||||'},{'|||| ||||'}];
                                                    alnstem = contains(stemalign(2,:),stempat);
                                                    if alnstem == true
                                                        champion = false;
                                                        break;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    if champion == true
                                        if ((i + Ls)-1) > strlength(plasmid)
                                            stop = mod(i + Ls,strlength(plasmid)) - 1;
                                        else
                                            stop = (i + Ls) - 1;
                                        end
                                        for p = lmin:lmax
                                            shiftplasmid = circshift(plasmid,-(i-1));
                                            primk = seqrcomplement(shiftplasmid(end-p+1:end));
                                            if ((primk(end) == 'C')||(primk(end) == 'G'))
                                                GC = count(primk,pat);
                                                gcr = GC/strlength(primk);
                                                if (gcr >= 0.40)&&(gcr <= 0.60)
                                                    [tmr,~,~] = primtemp(primk,rcon, monv, biv, catrat, a, d, g);
                                                    if ((tmr-tmf) >= dTm)&&((tmr-tmf) <= 1.25*dTm)
                                                        %% Checking homodimerization possiblity in the reverse primer
                                                        dgr = primdg(primk,monv,biv,tann,a,d,g,catrat);
                                                        dgrrmin = -RT*log(missf*(plascon/rcon)*exp(-dgr/RT));
                                                        [~, rAlignr] = swalign(seqcomplement(primk),seqrcomplement(primk),'GapOpen',1000,'Alphabet','NT');
                                                        dgrr = dGmax(rAlignr);
                                                        if dgrr > dgrrmin
                                                            %% Checking homo/multinucleotide repeats in the reverse primer
                                                            skipper = true;
                                                            for q = 1:length(nucleotides)
                                                                pattern = [nucleotides(q) '+'];
                                                                matches = regexp(primj, pattern, 'match');
                                                                numocc = cellfun(@length, matches);
                                                                if max(numocc) >= 4
                                                                    skipper = false;
                                                                    break
                                                                end
                                                            end
                                                            if skipper == true
                                                                for q = 1:length(dinucleotides.Value(:,1))
                                                                    pattern = ['(' dinucleotides.Value(q,:) ')' '+'];
                                                                    matches = regexp(primj, pattern, 'match');
                                                                    numocc = cellfun(@length, matches)./2;
                                                                    if max(numocc) >= 4
                                                                        skipper = false;
                                                                        break
                                                                    end
                                                                end
                                                            end
                                                            if skipper == true
                                                                winner = true;
                                                                %% Checking mispriming possibility by the reverse primer
                                                                lagplas = seqcomplement(circshift(plasmid,-(i-1)));
                                                                rlagplas = [lagplas(1:stop-p),repmat('X', 1, p),lagplas(stop+1:end)];
                                                                rGalstruct = localalign(rlagplas,seqreverse(primk),'numaln', 3,'GapOpen',1000,'Alphabet','NT');
                                                                rDalstruct = localalign(leadplas,primk,'numaln', 3,'GapOpen',1000,'Alphabet','NT');
                                                                rGalign = rGalstruct.Alignment;
                                                                rDalign = rDalstruct.Alignment;
                                                                dgrdmin = -RT*log(missf*exp(-dgr/RT)/(1+(rcon-plascon)*exp(-dgr/RT)));
                                                                for r = 1:numel(rGalign)
                                                                    dgrd = dGmax(rGalign{r});
                                                                    matches = regexp(strtrim(replace(rGalign{r}(2,:),'|','X')), ['X' '+'], 'match');
                                                                    numocc = cellfun(@length, matches);
                                                                    if (dgrd < dgrdmin)&&(max(numocc) >= 7)&&(rGalstruct.Stop(l,2)>=(strlength(primk)-3))
                                                                        winner = false;
                                                                        break;
                                                                    end
                                                                end
                                                                if winner == true
                                                                    dgrgmin = -RT*log(missf*exp(-dgr/RT)/(1+(rcon-plascon)*exp(-dgr/RT)));
                                                                    for s = 1:numel(rDalign)
                                                                        dgrg = dGmax(rDalign{s});
                                                                        matches = regexp(strtrim(replace(rGalign{s}(2,:),'|','X')), ['X' '+'], 'match');
                                                                        numocc = cellfun(@length, matches);
                                                                        if (dgrg < dgrgmin)&&(max(numocc) >= 7)&&(rDalstruct.Stop(l,2)>=(strlength(primk)-3))
                                                                            winner = false;
                                                                            break;
                                                                        end
                                                                    end
                                                                end
                                                                if winner == true
                                                                    legend = true;
                                                                    for t = 1:j % Stem starts here
                                                                        for u = stemmin:floor(0.5*(j-t)) % Length of stem
                                                                            stem53 = primj(t:t+u-1); % Sequence of stem
                                                                            for v = loopmin:(j-t-2*u)
                                                                                stem35 = primj(t+u+o:t+2*u+v-1); % Sequence of stem
                                                                                [~,stemalign] = swalign(stem53,seqrcomplement(stem35),'GapOpen',100,'Alphabet','NT');
                                                                                warning off
                                                                                if isempty(stemalign) == false
                                                                                    stempat = [{'||||||'},{'|||| |||'},{'||| ||||'},{'|||| ||||'}];
                                                                                    alnstem = contains(stemalign(2,:),stempat);
                                                                                    if alnstem == true
                                                                                        legend = false;
                                                                                        break;
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                    if legend == true
                                                                        %% Checking primer pair for heterodimerization possiblity
                                                                        [~, fAlignr] = swalign(primj,seqrcomplement(primk),'GapOpen',1000,'Alphabet','NT');
                                                                        dgfr = dGmax(fAlignr);
                                                                        dgfrmin = -RT*log(missf*(plascon/fcon)*exp(-dgr/RT));
                                                                        if dgfr > dgfrmin
                                                                            primer_f = primj;
                                                                            primer_r = primk;
                                                                            break;
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end
