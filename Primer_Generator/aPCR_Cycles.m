%% Collect User Input
clear
clc
close all
%% 
prompta = {'Enter plasmid file location:','Enter plasmid file name:'};
dlgtitle = 'Input Custom Genome Sequence File';
fieldsize = [1 50; 1 50];
definput = {'G:\Viral Genome Encapsulation Project Exp Data\aPCR_Trial\','D2_E_mut.txt'};
plas = inputdlg(prompta,dlgtitle,fieldsize,definput);
%% 
promptb = {'Is it double-stranded (2) or single stranded (1)?','Enter plasmid concentration (ng/uL):'};
dlgtitle = 'Input Plasmid Detail';
fieldsize = [1 50; 1 50];
definput = {'2','2'};
dna = inputdlg(promptb,dlgtitle,fieldsize,definput);
%% 
promptc = {'Enter the forward primer sequence:','Enter the reverse primer sequence:'};
dlgtitle = 'Input Primer Sequences';
fieldsize = [1 50;1 50];
definput = {'GCAGTTGGAAATGACACA','CCACCTTGTTTCAACGACCTCAC'};
primerseq = inputdlg(promptc,dlgtitle,fieldsize,definput);
%% 
promptd = {'Enter forward primer concentration (uM):','Enter forward/reverse primer concentration ratio:','Enter monovalent cation concentration (mM):',...
    'Enter bivalent cation concentration (mM):','Enter dNTPs concentration (mM):','Enter annealing temperature (oC):',...
    'Enter extension temperature (oC):','Enter the number of amplification cycles:'};
dlgtitle = 'Input Reaction Conditions';
fieldsize = [1 50; 1 50; 1 50; 1 50;1 50;1 50;1 50;1 50];
definput = {'1','65','60','2','0.3','60','65','30'};
react = inputdlg(promptd,dlgtitle,fieldsize,definput);
%% Input Processing
dir = plas(1); seqfile = plas(2);
strand = str2double(dna(1)); plaswcon = str2double(dna(2));
monv =  str2double(react(3))*10^(-3);
biv =  str2double(react(4))*10^(-3);
fcon = str2double(react(1))*10^(-6);
dNTP =  str2double(react(5))*10^(-3);
plasmid = char(extractFileText([dir{1},seqfile{1}]));
an = count(plasmid,'A'); tn = count(plasmid,'T'); cn = count(plasmid,'C'); gn = count(plasmid,'G');
plasmolwt = (an*313.2)+(tn*304.2)+(cn*289.2)+(gn*329.2)+79.0;
plascon = plaswcon*10^(-3)/(strand*plasmolwt);
rcon = fcon/str2double(react(2));
fprime = primerseq{1}; rprime = primerseq{2};
tann = str2double(react(6));
textn = str2double(react(7));
ncycle = str2double(react(8));
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
%% Setting-up the Cycle
leadcon(1,1) = plascon;
if strand == 2
    lagcon(1,1) = plascon;
else
    lagcon(1,1) = 0;
end
forcon(1,1) = fcon;
revcon(1,1) = rcon;
tann = (textn+tann)/2;
for i = 2:ncycle
    % fprintf('Cycle #%d.\n',i-1)
    thetaflag  = boundtheta(fprime,forcon(i-1,1),lagcon(i-1,1),monv,biv,tann,a,d,g,catrat);
    thetaflag(isnan(thetaflag)) = 0;
    thetarlead = boundtheta(rprime,revcon(i-1,1),leadcon(i-1,1),monv,biv,tann,a,d,g,catrat);
    pflag      = forcon(i-1,1)/(forcon(i-1,1)+leadcon(i-1,1));
    pflag(isnan(pflag)) = 0;
    prlead     = revcon(i-1,1)/(revcon(i-1,1)+leadcon(i-1,1));
    
    leadcon(i,1) = leadcon(i-1,1)  + pflag*lagcon(i-1,1)*thetaflag;
    forcon(i,1)  = forcon(i-1,1) - pflag*lagcon(i-1,1)*thetaflag;
    lagcon(i,1)  = lagcon(i-1,1)   + prlead*leadcon(i-1)*thetarlead;
    revcon(i,1)  = revcon(i-1,1) - prlead*leadcon(i-1)*thetarlead;
end
dcolor = distinguishable_colors(4);
figure
ax1 = axes;
semilogy(ax1,(10^9)*leadcon,'LineWidth',5,'Color',dcolor(1,:));
hold on
semilogy(ax1,(10^9)*lagcon,'LineWidth',2,'Color',dcolor(2,:));
hold on
semilogy(ax1,(10^9)*forcon,'LineWidth',2,'Color',dcolor(3,:));
hold on
semilogy(ax1,(10^9)*revcon,'LineWidth',2,'Color',dcolor(4,:));
hold off
ax1.YColor = 'k';
ax1.YLim = [0.1 10000];
yticks = ax1.YTick;
ax1.YTickLabel = log10(yticks);
ylabel(ax1,'log_{10}(Concentration (nM))');
grid on;
ampl = {'Leading Strand','Lagging Strand','Forward Primer','Reverse Primer'};
lgd = legend(ampl);
lgd.Location = 'southeast';
ax2 = axes('Position', ax1.Position, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', [], 'YScale', 'log');
ax2.YColor = 'k';
ylabel(ax2,'log_{10}(Concentration (nM))');
ax2.YLim = ax1.YLim;
yticks = ax2.YTick;
ax2.YTickLabel = log10(yticks);
set(ax1,"FontWeight","bold","FontSize",18);
set(ax2,"FontWeight","bold","FontSize",18);
title("Asymmetric PCR Amplification of ssDNA","FontSize",24);
subtitle("Reverse Primer = "+num2str((10^9)*rcon)+" nM, Forward Primer = "+num2str((10^6)*fcon)+" \muM, ds Template = "+num2str((10^9)*plascon,2)+" nM","FontSize",16);
xlabel('Cycle', 'Position', [mean(xlim), min(ylim) - 0.000005*range(ylim), 0]);
annotation('rectangle',[0 0 1 1],'Color','w');
