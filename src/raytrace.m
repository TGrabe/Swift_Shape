function [e_hitm, r_hitm,e_hitw,r_hitw,diffx,diffz, qval,trf1,trf2] = raytrace(rdata,ddata,Kx,Ky,alphat,material,Bqx,Bqy,Btx,Bty,etype)

% rdata: raytrace data/ parameters (plot, n_ray, n_trace, diff_size)
% ddata: light source and detection data (SPx, SPy, Psi1, Psi2, Dx1, Dy1,
% Dx2, Dy2)

% ## Define parameters
n1=rdata(1);
n2=rdata(2);
n_ray=rdata(3);
n_tracex=rdata(4);
n_tracez=rdata(5);
dsx=rdata(6);   
dsz=rdata(7);
plott=rdata(8);
n_tracez=n_tracez+1*(mod(n_tracez,2)==0);
n_tracex=n_tracex+1*(mod(n_tracex,2)==0);

u=1:n_tracex;
if n_tracex>1
    diffx=(-dsx/2+dsx/(n_tracex-1)*(u-1));
else
    diffx=0;
end

w=1:n_tracez;
if n_tracez>1
    diffz=(-dsz/2+dsz/(n_tracez-1)*(w-1));
else
    diffz=0;
end



%Define point source and detector lines
psi1=ddata(3);
psi2=ddata(4);
lambert=ddata(5);
exp_I=abs(ddata(6));
Dx1=ddata(7);
Dy1=ddata(8);
Dx2=ddata(9);
Dy2=ddata(10);
Dl=ddata(11);
if length(ddata)<12
DP=2;
else
    if ddata(12)<2
        DP=2;
    else
    DP=ddata(12);
    end
end
D_ang=atan2(Dy2-Dy1,Dx2-Dx1);
Dx3=Dx1-sin(D_ang)*Dl;
Dy3=Dy1+cos(D_ang)*Dl;
Dx4=Dx2-sin(D_ang)*Dl;
Dy4=Dy2+cos(D_ang)*Dl;

m=length(Kx(:,1));
r_hit=zeros(n_tracex,n_tracez);
e_hit=zeros(n_tracex,n_tracez);
trf1=zeros(n_tracex,n_tracez);
trf2=zeros(n_tracex,n_tracez);
qval=zeros(n_tracex,n_tracez,m);
Dn_hit=zeros(1,DP);

if Dx3==Dx4
    Dx=Dx3*ones(1,DP);
else
Dx=Dx3:(Dx4-Dx3)/(DP-1):Dx4;
end
if Dy3==Dy4
Dy=Dy3*ones(1,DP);
else
Dy=Dy3:(Dy4-Dy3)/(DP-1):Dy4;
end
for u=1:n_tracex
    
    for w=1:n_tracez
        
        % Define point source
        SPx=ddata(1)+diffx(u);
        SPy=ddata(2)+diffz(w);
      
        
        % Generate rays
        Sx=zeros(m+1,n_ray);
        Sy=zeros(m+1,n_ray);
        delta=zeros(m,n_ray);
        energy=zeros(m+1,n_ray);
        
        for k=1:n_ray
            delta(1,k)=psi2+(psi1-psi2)/(2*n_ray)+(psi1-psi2)/(n_ray)*(k-1);
        end
        Sx(1,:)=SPx;
        Sy(1,:)=SPy;
        if lambert==1
        energy(1,:)=abs((sin(delta(1,:))).*(cos(delta(1,:))));
        else
        energy(1,:)=abs(cos(delta(1,:)));
        end
        
        % Determine intersection points
        
        
        for j=2:m+1
            
            for k=1:n_ray
                if energy(j-1,k)==0
                    continue;
                end
                if j<=m
                    
                    
                    ihit=intersection(Kx{j},Ky{j},Sx(j-1,k),Sy(j-1,k),delta(j-1,k));
                    
                    
                    
                    if ihit>0
                        Kxhere=Kx{j}(ihit);
                        Kyhere=Ky{j}(ihit);
                        alphathere=alphat{j}(ihit);
                    else
                        if ~isempty(Bqx{j})
                            ihit=intersection(Bqx{j},Bqy{j},Sx(j-1,k),Sy(j-1,k),delta(j-1,k));
                            if ihit>0
                                Kxhere=Bqx{j}(ihit);
                                Kyhere=Bqy{j}(ihit);
                                alphathere=atan2(Bqy{j}(ihit+1)-Bqy{j}(ihit),Bqx{j}(ihit+1)-Bqx{j}(ihit));
                            end
                        end
                        if ihit==0
                            if ~isempty(Btx{j})
                                ihit=intersection(Btx{j},Bty{j},Sx(j-1,k),Sy(j-1,k),delta(j-1,k));
                                if ihit>0
                                    Kxhere=Btx{j}(ihit);
                                    Kyhere=Bty{j}(ihit);
                                    alphathere=atan2(Bty{j}(ihit+1)-Bty{j}(ihit),Btx{j}(ihit+1)-Btx{j}(ihit));
                                end
                            end
                        end
                    end
                    
                    
                    if ihit>0
                        
                        Sx(j,k)=Kxhere;
                        Sy(j,k)=Kyhere;
                        
                        if material(j-1)==2
                            nc1=n2; nc2=n1;
                        else
                            nc1=n1; nc2=n2;
                        end
                        
                        % outgoing angle for reflection
                        if etype(j)==-1
                            delta(j,k)=2*alphathere-delta(j-1,k);
                            esuv=energysurvived(delta(j-1,k),alphathere,nc1,nc2,-1);
                            energy(j,k)=energy(j-1,k)*esuv;
                            if material(j)==2 && esuv < 1
                                trf1(u,w)=1;
                            end
                            qval(u,w,j)=qval(u,w,j)+esuv;
                        else
                            % outgoing angle for refraction
                            delta(j,k)=refrangle(delta(j-1,k),alphathere,nc1,nc2);
                            
                            esuv=energysurvived(delta(j-1,k),alphathere,nc1,nc2,1);
                            energy(j,k)=energy(j-1,k)*esuv;
                            qval(u,w,j)=qval(u,w,j)+esuv;
                            if ~isreal(delta(j,k))
                                energy(j,k)=0;
                                trf2(u,w)=1;
                            end
                        end
                    else
                        energy(j,k)=0;
                    end
                    
                else
                    ihit1=intersection([Dx1, Dx2],[Dy1, Dy2],Sx(j-1,k),Sy(j-1,k),delta(j-1,k));
                    ihit2=intersection(Dx,Dy,Sx(j-1,k),Sy(j-1,k),delta(j-1,k));
                    if (ihit1>0) && (ihit2>0)
                        energy(j,k)=energy(j-1,k);
                        Sx(j,k)=Dx(ihit2);
                        Sy(j,k)=Dy(ihit2);
                        Dn_hit(ihit2)=Dn_hit(ihit2)+0.5;
                        Dn_hit(DP-(ihit2-1))=Dn_hit(DP-(ihit2-1))+0.5;
                    end
                    
                end
            end
            
            if j<=m
                qval(u,w,j)=qval(u,w,j)/length(energy(j,energy(j,:)>0));
            end
        end
        %{f        
        e_hit(u,w)=sum(energy(end,:))/sum(energy(1,:));
        for k=1:n_ray
        yeshit(k)=1*(energy(end,k)>0);
        end
        r_hit(u,w)=sum(yeshit.*energy(1,:))/sum(energy(1,:));
        %}
        if plott==1
        if n_tracex==1  && n_ray<=50 && (w==1 || w==n_tracez/2+0.5 ||w==n_tracez)
    figure(1)
    raycolor=[50,50,50]/255;
    red=[223,5,5]/255;
    blue=[104,165,253]/255;
    if n_tracez==1
    figure
    else
       if w==1
       subplot(1,3,1);
       end
       if w==n_tracez/2+0.5
       subplot(1,3,2);
       end
       if w==n_tracez
       subplot(1,3,3);
       end
    end
    hold on
    grid on
    for j=1:m
        for k=1:n_ray
            if energy(j,k)>1 && j>0
                plot(Sx(j,k),Sy(j,k),'*');
            end
            if energy(j+1,k)>0 %&& j<m
                plot([Sx(j,k), Sx(j+1,k)],[Sy(j,k) Sy(j+1,k)],'Color',raycolor,'LineWidth',1);
            else
                hll=40; % help line length (temporary)
                if energy(j+1,k)>0
                    plot([Sx(j,k), Sx(j,k)+hll*cos(delta(j,k))],[Sy(j,k), Sy(j,k)+hll*sin(delta(j,k))],'Color',raycolor,'LineWidth',0.5);
                else
                    if energy(j+1,k)==0 && energy(j,k)>0
                        plot([Sx(j,k), Sx(j,k)+hll*cos(delta(j,k))],[Sy(j,k), Sy(j,k)+hll*sin(delta(j,k))],'Color',red,'LineWidth',0.5);
                    end
                end
            end
        end
    end
    for j=1:m
        % Element curves
        if j>1
            if etype(j)==1
                plot(Kx{j},Ky{j},'Color',blue,'LineWidth',1.5);
            end
            
            if etype(j)==-1
                plot(Kx{j},Ky{j},'Color',red,'LineWidth',1.5);
            end
            
            if ~isempty(Btx(j))
                plot(Btx{j},Bty{j},'g','LineWidth',2);
            end
            if ~isempty(Bqx(j))
                plot(Bqx{j},Bqy{j},'g','LineWidth',2);
            end
        end
    end
    % Rotational symmetry line
    plot([0 0],[-dsz*3/4, Dy4*2],'w','LineWidth',1.2);
    plot([0 0],[-dsz*3/4, Dy4*2],'-.k','LineWidth',1.2);
    plot([Dx1,Dx2],[Dy1,Dy2],'k','LineWidth',2);
    plot([Dx1,Dx2],[Dy1+Dl,Dy2+Dl],'k','LineWidth',2);
    xlabel('$r$','interpreter','latex');
    ylabel('$z$','interpreter','latex');    
    plot(SPx,SPy,'*','Color',red);
    
    ylim([-dsz*3/4,max([Dy1,Dy2,Dy3,Dy4])*1.05]);
    xlim([1.5*min([Dx1,Dx2,Dx3,Dx4]),-1.5*min([Dx1,Dx2,Dx3,Dx4])+max(cell2mat(cellfun(@max,Kx,'UniformOutput',0)))]);
    set(findall(gcf,'-property','FontSize'),'FontSize',10)
    set(gcf, 'units', 'centimeters');
    set(gcf, 'Position', [0, 0, 15.5, 8.5]);    
    if n_tracez==1 || w==n_tracez     
    h = zeros(7, 1);
        h(1) = plot(NaN,NaN,'Color',blue,'LineWidth',1.5);
        h(2) = plot(NaN,NaN,'Color',red,'LineWidth',1.5);
        h(3) = plot(NaN,NaN,'Color',raycolor,'LineWidth',0.5);
        h(4) = plot(NaN,NaN,'Color',red,'LineWidth',0.5);
        h(5) = plot(NaN,NaN,'-.k','LineWidth',1.2);
        h(6) = plot(NaN,NaN,'*','Color',red);
        h(7) = plot(NaN,NaN,'k','LineWidth',2);


        
        leg=legend(h, 'Refr.','Refl.','Ray','Lost ray','Opt. axis','Source','Detector','interpreter','latex','FontSize',8,'Location', 'Best');
        leg.ItemTokenSize = [15,20];
    end
    %daspect([1 1 1])
    % plot detector intensity
    %figure
    %bar(Dx,Dn_hit);
        end
        end
    end
    
end

e_hitm=e_hit;
r_hitm=r_hit;


imid=(length(e_hit(:,1)))/2+0.5;
for i=1:imid-1
    e_hitm(i,:)=(e_hit(i,:)+e_hit(2*imid-i,:))/2;
    e_hitm(2*imid-i,:)=e_hitm(i,:);
    
    r_hitm(i,:)=(r_hit(i,:)+r_hit(2*imid-i,:))/2;
    r_hitm(2*imid-i,:)=r_hitm(i,:);
end
   e_hitm(imid,:)=e_hit(imid,:);
   r_hitm(imid,:)=r_hit(imid,:);

e_hitw=e_hitm;
r_hitw=r_hitm;

A=zeros(1,imid-1);
wght=zeros(1,imid-1);


for i=1:imid-1
A(i)=pi*((diffx(imid+i)-dsx/(n_tracex-1)/2)^2-(diffx(imid+i)+dsx/(n_tracex-1)/2)^2);
    %A(i)=pi*(diffx(imid+i)^2-diffx(imid+i-1)^2);   
I=(1-(exp_I>0)*(1:i).^exp_I/(i+1)^exp_I); %quadratic
    wght(1:i)=A(1:i).*I(1:i);
c=1/sum(wght(1:i));
for j=1:n_tracez
r_hitw(imid+i,j)=c*sum(wght(1:i).*r_hitm(imid+1:imid+i,j)');
r_hitw(imid-i,j)=r_hitw(imid+i,j);

e_hitw(imid+i,j)=c*sum(wght(1:i).*e_hitm(imid+1:imid+i,j)');
e_hitw(imid-i,j)=e_hitw(imid+i,j);
end
end
%{s
%% plot

% Give one value
if n_tracez==1 && n_tracex==1
    disp('Ideal rays hitting the detector');
    disp([num2str(r_hit*100),' %']);
    
    disp('Energy hitting the detector');
    disp([num2str(e_hit*100),' %']);
end

if plott==1
pwidth=8.5;
pheight=6;
% 3D-plot
if n_tracez>1 && n_tracex>1
edgealpha=0.3;
if n_tracex>80 || n_tracez>80
edgealpha=0.1;
end
figure
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
contourf(diffx,diffz,r_hitm'*100);
colormap(parula(15));
cb=colorbar('TickLabelInterpreter','latex','FontSize',10);
%surf(diffx,diffz,r_hitm'*100,'EdgeColor','k','FaceColor','interp','EdgeAlpha',edgealpha);
xlabel('Light source position $x_s$','interpreter','latex');ylabel('Light source position $z_s$','interpreter','latex');
ylabel(cb,'Rays hit in \%','Interpreter','latex');

figure
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
contourf(diffx,diffz,e_hitm');
colormap(parula(15));
cb=colorbar('TickLabelInterpreter','latex','FontSize',10);
%surf(diffx,diffz,e_hitm','EdgeColor','k','FaceColor','interp','EdgeAlpha',edgealpha);
xlabel('Light source position $x_s$','interpreter','latex');ylabel('Light source position $z_s$','interpreter','latex');
ylabel(cb,'Surface seq. efficiency $_{\mathrm{2D}}v_{\mathrm{eff}}$','Interpreter','latex');

figure
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
contourf(diffx,diffz,r_hitw'*100);
colormap(parula(15));
cb=colorbar('TickLabelInterpreter','latex','FontSize',10);
%surf(diffx,diffz,r_hitw'*100,'EdgeColor','k','FaceColor','interp','EdgeAlpha',edgealpha);
xlabel('Light source radius $r_s$','interpreter','latex');ylabel('Light source position $z_s$','interpreter','latex');
ylabel(cb,'Rays hit in \%','Interpreter','latex');

figure
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
contourf(diffx,diffz,e_hitw');
colormap(parula(15));
cb=colorbar('TickLabelInterpreter','latex','FontSize',10);
%surf(diffx,diffz,e_hitw','EdgeColor','k','FaceColor','interp','EdgeAlpha',edgealpha);
xlabel('Light source radius $r_s$');ylabel('Light source position $z_s$','interpreter','latex');
ylabel(cb,'Surface seq. efficiency $_{\mathrm{2D}}v_{\mathrm{eff}}$','Interpreter','latex','interpreter','latex');

end
% 2D-z-plot
if n_tracez>1 && n_tracex==1

figure
plot(diffz, e_hitm,'k');
grid on
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
xlabel('Light source position $z_s$','interpreter','latex');
ylabel('Surface seq. efficiency $_{\mathrm{2D}}v_{\mathrm{eff}}$','interpreter','latex');
ylim([0,1]);

figure
plot(diffz, r_hitm*100,'k');
grid on
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
xlabel('Light source position $z_s$','interpreter','latex');
ylabel('Rays hit in $r_\mathrm{hit}$ \%','interpreter','latex');
ylim([0,105]);

end
if n_tracez==1 && n_tracex>1
figure
plot(diffx, e_hitw,'k');

set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
xlabel('Light source radius $r_s$','interpreter','latex');
ylabel('Surface seq. efficiency $_{\mathrm{2D}}v_{\mathrm{eff}}$','interpreter','latex');
ylim([0,1]);
figure
plot(diffx, r_hitw*100,'k');
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0, 0, pwidth, pheight]);
xlabel('Light source radius $r_s$ ','interpreter','latex');
ylabel('Rays hit $r_\mathrm{hit}$ in \%','interpreter','latex');
ylim([0,105]);

end



end

end


