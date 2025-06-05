%% Script for Optimization of aPCR Operating Conditions
clc
clear
clear cache

prompta = {'Enter plasmid file location:','Enter plasmid file name:',...
    'Is it double-stranded (2) or single stranded (1)?',['Enter the' ...
    ' desired length of scaffold:']};
dlgtitle = 'Input Custom Genome Sequence & Scaffold Details';
fieldsize = [1 50; 1 50; 1 50; 1 50];
definput = {'G:\Viral Genome Encapsulation Project Exp Data\aPCR_Trial\'...
    , 'D2_E_mut.txt', '2', '3700'};
plas = inputdlg(prompta,dlgtitle,fieldsize,definput);
promptb = {'Enter monovalent cation concentration (mM):',['Enter ' ...
    'bivalent cation concentration(mM):'],['Enter dNTPs concentration' ...
    ' (mM):']};
dlgtitle = 'Input Reaction Conditions';
fieldsize = [1 50; 1 50; 1 50];
definput = {'60','2','0.3'};
react = inputdlg(promptb,dlgtitle,fieldsize,definput);
promptc = {'Enter desired minimum length of primers (nt):',['Enter' ...
    ' desired maximum length of primers (nt):']};
dlgtitle = 'Input Primer Design Parameters';
fieldsize = [1 50; 1 50];
definput = {'18','22'};
descon = inputdlg(promptc,dlgtitle,fieldsize,definput);

Tmin = [60,61,62,63,64];
dTm = [0.001,1,2,3,5];
Tann = [57,59,61,63,65];
fcon = [0.5,0.75,1,1.5,2];
plaswcon = [0.5,1,1.5,2,3];
ftor = [35,50,65,80,100];
cycle = [25,30,35,40,45];

BigT = table();

parfor i = 1:numel(plaswcon)
    for j = 1:numel(fcon)
        for k = 1:numel(ftor)
            for l = 1:numel(Tmin)
                for m = numel(dTm)
                    for n = 1:numel(Tann)
                        [primer_f, primer_r] = apcrprimfun(plas,react,descon,plaswcon(i),fcon(j),ftor(k),Tmin(l),dTm(m),Tann(n));
                        for o = 1:numel(cycle)
                            [samplicon,damplicon] = apcrcycfun(plas,react,descon,plaswcon(i),fcon(j),ftor(k),primer_f, primer_r,Tann(n),cycle(o));
                            tabrow = table(plaswcon(i),fcon(j),ftor(k),primer_f.Tm,primer_r.Tm,Tann(n),cycle(o),'VariableNames',{'Plasmid_WtCon','F_Prime_Conc','F_by_R','F_Tm','R_Tm','T_ann','Cycles'});
                            BigT = [BigT;tabrow];
                        end
                    end
                end
            end
        end
    end
end
