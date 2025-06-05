function dGm = dGmax(Alignment,monv,biv,tann,a,d,g,catrat)
            alignchar = Alignment(2,:);
            alignchar = replace(alignchar,'|','X');
            alignbas =  Alignment(1,:);
            matchchar = 'X';
            pattern = [matchchar '+'];
            trimchar = strtrim(alignchar);
            matches = regexp(trimchar, pattern, 'match');
            startpos = regexp(trimchar, pattern, 'start');
            offset = length(alignchar) - length(trimchar);
            startpos = startpos + offset;
            numocc = cellfun(@length, matches);
            [sorted, ind] = sort(numocc,'descend');
            u = 1;
            dGm = 0;
            if sorted(1) >= 3
                while (u <= numel(sorted))&&(sorted(u) >= 3)
                    dG = primdg37(string(alignbas(1,startpos(ind(u)):startpos(ind(u))+numocc(ind(u))-1)),monv,biv,tann,a,d,g,catrat);
                    u = u + 1;
                    if dG < dGm
                        dGm = dG;
                    end
                end
            end
end