%   ######################################################
%   #     swiftshape - Henrik Krauss/Tobias Grabe 2021   #
%   ######################################################
%   #      CONTROL - w/ Optimisation                     #
%   ######################################################
%        Setup with example collimator
%%
clc; clear; close all;

%%   Optimisation function setup

%   Variable optimisation parameters
E_in=[0,pi/4,pi/3];

f=@opti;
opti(E_in,1,0);
options = optimset('PlotFcns',@optimplotfval,'TolX',0.001);

%%   Optimisation loop

[E_in_opt,fval,exitflag,output]=fminsearch(f,E_in,options);
disp('DONE');
% Plot optimised results
figure
opti(E_in_opt,1,1)

function [q] = opti(E_in,plott,export)
if ~exist('plott','var')
    plott=0;
end
if ~exist('export','var')
    export=1;
end

n1=1;
n2=1.45; %n_pmma(785); %=n_objet()

%   Shape generation

%   >>>>>    Shape generation input    <<<<<

%   Iteration parameters
%       n1, n2, t, lmax, p
i_data=[n1, n2,	0.001,	500,	plott];

%   Surface sequence 1

%   Light source data
%       phi1, phi2, Tx1, Ty1
l_data1=[pi/2,  E_in(3),	0,	0];

%   Element construction data
%       etype, cline, wto, sdef, rtype, Px, Py, theta, bQ, bT
E_data1=[0,	1,	0,	0,	-1,	0,	0,          pi/2,	0,	0;
         1,	1,	0,	10,	-1,	0,  E_in(1),	pi/2,	0,	0;
         1,	2,	0,	3,	0,	0,	0,          pi/2,	0,	0;
        ];

%   Execution of shape generation function 1
[Kx1,Ky1,alphat1,material1,Bqx1,Bqy1,Btx1,Bty1,Qx1,Qy1,Tx1,Ty1]=shapegen(l_data1,E_data1,i_data);

%   Light source data
%       phi1, phi2, Tx1, Ty1
l_data2=[E_in(3),   pi/6,	4,	0];

%   Element construction data
%       etype, cline, wto, sdef, rtype, Px, Py, theta, bQ, bT
E_data2=[0,	1,	0,	0,      -1,	0,	0,	pi/2,	0,	0;
         1,	1,	1,	Tx1(2),	0,	0,	0,  E_in(2),0,	0;
        -1,	2,	0,	4,      0,	0,	0,	pi/2,   0,	0;
        1,	1,	0,	1,      0,	0,	0,	pi/2,   0,	0;
        ];

%   Execution of shape generation function 2
[Kx2,Ky2,alphat2,material2,Bqx2,Bqy2,Btx2,Bty2,Qx2,Qy2,Tx2,Ty2]=shapegen(l_data2,E_data2,i_data);

%   >>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<

%   Raytracing
%   >>>>>      Raytracing input      <<<<<

%   Raytracing parameters
%       n1, n2, nray, ndeltax, ndeltaz, deltax, deltaz, p
rdata=[n1,  n2,	100,	11,	11,	4,	6,	0];

%   Light source and detector data
%       Sx, Sy, psi1, psi2, dlambert, eI, Dx1, Dz1, Dx2, Dz2, Dl
ddata=[0,	0,	pi/2,	0,	1,	2,	-16,	50,	16,	50,	20,	10];

%   >>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<

%   Executions of raytracing function
[~,r_hit1,e_hitw1,r_hitw1,diffx1,diffz1]=raytrace(rdata,ddata,Kx1,Ky1,alphat1,material1,Bqx1,Bqy1,Btx1,Bty1,E_data1(:,1));
[~,r_hit2,e_hitw2,r_hitw2,diffx2,diffz2]=raytrace(rdata,ddata,Kx2,Ky2,alphat2,material2,Bqx2,Bqy2,Btx2,Bty2,E_data2(:,1));


if plott==1
    pwidth=8;
    pheight=6;
    figure
    contourf(diffx1,diffz1,(r_hitw1+r_hitw2)'*100);
    colormap(parula(15));
    cb=colorbar('TickLabelInterpreter','latex','FontSize',10);
    set(gcf, 'units', 'centimeters');
    set(gcf, 'Position', [0, 0, pwidth, pheight]);
    xlabel('Light source radius $r_s$');ylabel('Light source position $z_s$');
    ylabel(cb,'Rays hit in \%','Interpreter','latex');
    zlim([0,inf]);
    figure
    contourf(diffx1,diffz1,(e_hitw1+e_hitw2)');
    colormap(parula(15));
    cb=colorbar('TickLabelInterpreter','latex','FontSize',10);
    set(gcf, 'units', 'centimeters');
    set(gcf, 'Position', [0, 0, pwidth, pheight]);
    xlabel('Light source radius $r_s$');ylabel('Light source position $z_s$');
    ylabel(cb,'Surface seq. efficiency $_{\mathrm{2D}}v_{\mathrm{eff}}$','Interpreter','latex');
end

% Optical data export
if export==1
    recycle on % Send to recycle bin instead of permanently deleting.
    n_s=10; % Save every n_s value
    for j=2:length(E_data1(:,1))
        filename=['sequence1_',num2str(j-1), '.xls'];
        delete(filename);
        
        outputArray=[Kx1{j}(1:n_s:end).', Ky1{j}(1:n_s:end).'];
        writematrix(outputArray,filename);
    end
    for j=2:length(E_data2(:,1))
        filename=['sequence12_',num2str(j-1), '.xls'];
        delete(filename);
        
        outputArray=[Kx2{j}(1:n_s:end).', Ky2{j}(1:n_s:end).'];
        writematrix(outputArray,filename);
    end
end

% Optimisation value calculation
q=1-mean(e_hitw1+e_hitw2,'all');

end


%coded by Hernik Krauss