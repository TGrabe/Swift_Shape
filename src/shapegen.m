%Henrik Krauß,Tobias Grabe,2019

function [Kx,Ky,alphat,material,Bqx,Bqy,Btx,Bty,Qx,Qy,Tx,Ty] = shapegen(l_data,E_data,i_data)

%% #### Iteration and material parameters ####

%% Initialisation and declaration

% Iteration data
n1=i_data(1);
n2=i_data(2);
t=i_data(3);
l_max=i_data(4);
plott=i_data(5);
% Elemente data
etype=E_data(:,1);
chln=E_data(:,2);
whto=E_data(:,3);
sxy=E_data(:,4);
rtype=E_data(:,5);
Px=E_data(:,6);
Py=E_data(:,7);
theta=E_data(:,8);

if length(E_data(1,:))>8
Bql=E_data(:,9);
else
Bql=zeros(1,length(E_data(:,1)));   
end
if length(E_data(1,:))>9
Btl=E_data(:,10);
else
Btl=zeros(1,length(E_data(:,1)));
end

n=ceil((l_max+max(Bql)+max(Btl))/t);
m=length(etype); % m: Number of shapes that are calculated

% Construction parameters
Qx=zeros(m,1);Qy=zeros(m,1);
Tx=zeros(m,1);Ty=zeros(m,1);
s=zeros(m,1);
material=zeros(m,1);
phi1=zeros(m,1);
phi2=zeros(m,1);
b=zeros(m,1);
Kxm=nan(m,n); Kym=nan(m,n);
Kx=cell(m,1); Ky=cell(m,1);
Bqx=cell(m,1); Bqy=cell(m,1);
Btx=cell(m,1); Bty=cell(m,1);
nBq=zeros(m);
nBt=zeros(m);
alpha1m=nan(m,n); alpha2m=nan(m,n);
alpha1=cell(m,1); alpha2=cell(m,1);
alphatm=nan(m,n);
alphat=cell(m,1);

% Light source specific construction parameters
material(1)=(etype(1)==1)+1; % material: defines refractive index at the moment after passing element j
Qx(1)=Px(1);
Qy(1)=Py(1);
phi1(1)=l_data(1);
phi2(1)=l_data(2);
if rtype(1)==0
Tx(1)=l_data(3);
Ty(1)=l_data(4);
else
Tx(1)=Px(1);
Ty(1)=Py(1);    
end


%% #### Calculation ####


for j=1:m
    nBq(j)=ceil(max(Bql(j))/t);
    nBt(j)=ceil(max(Btl(j))/t);
    
    startx=0;starty=0;phihere=0;
    if j>1
        
        if chln(j)==1
            phihere=phi1(j-1);
            phithere=phi2(j-1);
            
            if chln(j-1)==1
                startx=Qx(j-1); starty=Qy(j-1);
                oppositex=Tx(j-1); oppositey=Ty(j-1);
            else
                startx=Tx(j-1); starty=Ty(j-1);
                oppositex=Qx(j-1);oppositey=Qy(j-1);
            end
        end
        if chln(j)==2
            phihere=phi2(j-1);
            phithere=phi1(j-1);
            
            if chln(j-1)==1
                startx=Tx(j-1); starty=Ty(j-1);
                oppositex=Qx(j-1); oppositey=Qy(j-1);
                
            else
                startx=Qx(j-1); starty=Qy(j-1);
                oppositex=Tx(j-1); oppositey=Ty(j-1);
            end
        end
        
        if whto(j)==0
            s(j)=sxy(j);
        else
            if whto(j)==1
                s(j)=(sxy(j)-startx)/cos(phihere);
            else
                s(j)=(sxy(j)-starty)/sin(phihere);
            end
            if s(j)>l_max
                s(j)=sxy(j);
            end
        end
        Qx(j)=startx+cos(phihere)*s(j);
        Qy(j)=starty+sin(phihere)*s(j);
        
        % define material of next element
        if etype(j)==1
            if material(j-1)==1
                material(j)=2;
            else
                material(j)=1;
            end
        else
            material(j)=material(j-1);
        end
        
        Kxm(j,1)=Qx(j);
        Kym(j,1)=Qy(j);
    end
    
    % #### Curve calculation
    constr=1; % construct 1: Kx, Ky 2: Btx, Bty 3: Bqx, Bqy
    skip=0;
    if j>1
        for i=1:n
            
            % Calculate alpha1m (incident angle)
            if rtype(j-1)==-1
                alpha1m(j,i)=atan2(Kym(j,i)-Py(j-1),Kxm(j,i)-Px(j-1));
            end
            if rtype(j-1)==0
                alpha1m(j,i)=theta(j-1);
            end
            if rtype(j-1)==1
                alpha1m(j,i)=atan2(Py(j-1)-Kym(j,i),Px(j-1)-Kxm(j,i));
            end
            
            %{
if alpha2m(j,i)<0
    alpha2m(j,i)=alpha2m(j,i)+2*pi;
end
            %}
            
            %{
while alpha1m(j,i)<0
    alpha1m(j,i)=alpha1m(j,i)+2*pi;
end
            %}
            
            % Calculate alpha2m (outgoing angle)
            if rtype(j)==0
                alpha2m(j,i)=theta(j);
            end
            if rtype(j)==1
                alpha2m(j,i)=atan2((Py(j)-Kym(j,i)),(Px(j)-Kxm(j,i)));
            end
            if rtype(j)==-1
                alpha2m(j,i)=atan2((Kym(j,i)-Py(j)),(Kxm(j,i)-Px(j)));
            end
            
            %{
while alpha2m(j,i)<0
    alpha2m(j,i)=alpha2m(j,i)+2*pi;
end
            %}
            
            % Changing refractive indices depending on material
            if material(j-1)==2
                nc1=n2; nc2=n1;
                vz=-1;
            else
                nc1=n1; nc2=n2;
                vz=1;
            end
            
            % Transmission surface angle
            if etype(j)==1
                alphatm(j,i)=atan2(vz*(nc1*cos(alpha1m(j,i))-nc2*cos(alpha2m(j,i))),vz*(nc2*sin(alpha2m(j,i))-nc1*sin(alpha1m(j,i))));
                %hier war mal mod(...,2pi)
                
                % Reflection surface angle
            else
                alphatm(j,i)=(alpha2m(j,i)+alpha1m(j,i))/2;
            end
            
            %Vector to next iteration coordinates (does not save all data)
            if i>1
                
            end
            
            if i==1
                corr=0;
                lcheck=s(j)+5*t;
                if ~inpolygon(Kxm(j,i)+t*cos(alphatm(j,i)),Kym(j,i)+t*sin(alphatm(j,i)),[startx,startx+cos(phihere)*lcheck,oppositex+cos(phithere)*lcheck,oppositex],[starty,starty+sin(phihere)*lcheck,oppositey+sin(phithere)*lcheck,oppositey])
                    corr=pi;
                end
            end
            rtx=cos(alphatm(j,i)+corr);
            rty=sin(alphatm(j,i)+corr);
            
            
            
            % Calculation of next point on
            Kxm(j,i+1)=Kxm(j,i)+t*rtx;
            Kym(j,i+1)=Kym(j,i)+t*rty;
            
            if constr==3
                if i==iKdone+nBt(j)+nBq(j)
                    iBqdone=i;
                    break;
                end
                
            else
                if constr==2
                    if i==iKdone+nBt(j)
                        iBtdone=i;
                        if Bql(j)==0
                            break;
                        end
                        constr=3;
                        Kxm(j,i+1)=Qx(j,1);
                        Kym(j,i+1)=Qy(j,1);
                        corr=corr-pi;
                    end
                end
            end
            
            if constr==1
                % ## Criteria to end Iteration
                % point source
           if skip == 0
                if (rtype(j-1)~=0)
                    ang_cov1=atan2(Kym(j,1)-Py(j-1),Kxm(j,1)-Px(j-1));
                    ang_cov2=atan2(Kym(j,i)-Py(j-1),Kxm(j,i)-Px(j-1));
                    ang_cov=abs(ang_cov1-ang_cov2); %Angle covered so far
                    
                    ang_dist=angbtwn(phi1(j-1),phi2(j-1));
                    %abs(2*pi*(phi1(j-1)<phi2(j-1))+phi1(j-1)-phi2(j-1)) old one
                    
                    if ang_cov>=ang_dist
                        iKdone=i;
                        if Btl(j)==0 && Bql(j)==0
                            break;
                        end
                        if Btl(j)>0
                            constr=2;
                            Kxm(j,i+1)=Kxm(j,i);
                            Kym(j,i+1)=Kym(j,i);
                        else
                            iBtdone=i;
                            constr=3;
                            Kxm(j,i+1)=Qx(j,1);
                            Kym(j,i+1)=Qy(j,1);
                            corr=corr-pi;
                        end
                    end
                    skip=floor(orthogonal_distance(Kxm(j,i),Kym(j,i),oppositex,oppositey,phithere)/t);
                else
                
                %parallel
                    w=orthogonal_distance(Kxm(j,i),Kym(j,i),Qx(j),Qy(j),theta(j-1));
                    if w>=b(j-1)
                        iKdone=i;
                        if Btl(j)==0 && Bql(j)==0
                            break;
                        end
                        if Btl(j)>0
                            constr=2;
                            Kxm(j,i+1)=Kxm(j,i);
                            Kym(j,i+1)=Kym(j,i);
                        else
                            iBtdone=i;
                            constr=3;
                            Kxm(j,i+1)=Qx(j,1);
                            Kym(j,i+1)=Qy(j,1);
                            corr=corr-pi;
                        end
                    end
                    skip=floor((b(j-1)-w)/t);
                end
             end
                        skip=skip-1*(skip>0);
            end
            
            if i==n
                iKdone=i;
            end
        end
        
        % Save T-Point and set parameters (end of calculated curve)
        Tx(j)=Kxm(j,iKdone);
        Ty(j)=Kym(j,iKdone);
        if chln(j)==1
            phi1(j)=alpha2m(j,1);
            phi2(j)=alpha2m(j,iKdone);
        end
        if chln(j)==2
            phi1(j)=alpha2m(j,iKdone);
            phi2(j)=alpha2m(j,1);
        end
        
        Kx{j}=Kxm(j,1:iKdone);
        Ky{j}=Kym(j,1:iKdone);
        if Btl(j)>0
            Btx{j}=Kxm(j,iKdone+1:iBtdone);
            Bty{j}=Kym(j,iKdone+1:iBtdone);
        end
        if Bql(j)>0
            Bqx{j}=Kxm(j,iBtdone+1:iBqdone);
            Bqy{j}=Kym(j,iBtdone+1:iBqdone);
        end
        alpha1{j}=alpha1m(j,1:iKdone);
        alpha2{j}=alpha2m(j,1:iKdone);
        alphat{j}=alphatm(j,1:iKdone);
    end
    b(j)=orthogonal_distance(Qx(j),Qy(j),Tx(j),Ty(j),theta(j));
    
end

%% #### Plot ####
if plott==1
blue=[104,165,253]/255;
red=[223,5,5]/255;
marksize=6;

    hold on;
    grid on;
    lmax_x=max([Qx;Tx])*1.2;
    lmax_y=max([Qy;Ty])*1.2;
    %axis([-lmax_x/1.2*0.1,lmax_x/1.2*1.1,-lmax_y/1.2*0.05,lmax_y/1.2*1.15]);
    axis([-lmax_y/1.2*0.05,lmax_y/1.2*1.15,-lmax_y/1.2*0.05,lmax_y/1.2*1.15]);
    daspect([1 1 1])
    set(gcf, 'units', 'centimeters');
    set(gcf, 'Position', [0, 0, 7.5, 7.5]);
    
    for j=1:m
        
        % Element curves
        if j>1
            if etype(j)==1
                %fajsdfsad gleich hier hin schreiben
                plot(Kx{j}(1:end),Ky{j}(1:end),'Color',blue,'LineWidth',1.5);
            end
            
            if etype(j)==-1
                plot(Kx{j}(1:end),Ky{j}(1:end),'Color',red,'LineWidth',1.5);
            end
            
            if Btl(j)>0
                plot(Btx{j}(1:end),Bty{j}(1:end),'g','LineWidth',1.5);
            end
            if Bql(j)>0
                plot(Bqx{j}(1:end),Bqy{j}(1:end),'g','LineWidth',1.5);
            end
            
        end
        
        % P,Q,T-points
        if j>1
            plot(Tx(j),Ty(j),'xk','MarkerSize',marksize,'LineWidth',1);
            plot(Qx(j),Qy(j),'ok','MarkerSize',marksize,'LineWidth',1);
            plot(Qx(j),Qy(j),'.k','MarkerSize',marksize,'LineWidth',1);

        end

        
        % ## Angle lines
        if j<m
            if (chln(j)==1)
                lt1='--';lt2=':';
            else
                lt1=':';lt2='--';
            end
            if chln(j)==chln(j+1)
                plot([Qx(j),Qx(j+1)],[Qy(j),Qy(j+1)],[lt1,'k'],'LineWidth',0.7);
                plot([Tx(j),Tx(j+1)],[Ty(j),Ty(j+1)],[lt2,'k'],'LineWidth',0.7);
            else
                plot([Qx(j),Tx(j+1)],[Qy(j),Ty(j+1)],[lt1,'k'],'LineWidth',0.7);
                plot([Qx(j+1),Tx(j)],[Qy(j+1),Ty(j)],[lt2,'k'],'LineWidth',0.7);
           end
        end
        
        
        if j==m
            hll=lmax_y; %max(Ty(end),Qy(end))/2; %Helping line length
            if chln(j)==1
                plot([Qx(j) Qx(j)+hll*cos(phi1(j))],[Qy(j) Qy(j)+hll*sin(phi1(j))],'--k','LineWidth',0.7);
                plot([Tx(j) Tx(j)+hll*cos(phi2(j))],[Ty(j) Ty(j)+hll*sin(phi2(j))],':k','LineWidth',0.7);
            end
            if chln(j)==2
                plot([Qx(j) Qx(j)+hll*cos(phi2(j))],[Qy(j) Qy(j)+hll*sin(phi2(j))],':k','LineWidth',0.7);
                plot([Tx(j) Tx(j)+hll*cos(phi1(j))],[Ty(j) Ty(j)+hll*sin(phi1(j))],'--k','LineWidth',0.7);
            end
        end
        
        
     
        
       
    end
       % Rotational symmetry line
        plot([0 0],[0, lmax_y],'w','LineWidth',1.2);
        plot([0 0],[0, lmax_y],'-.k','LineWidth',1.2);
        xlabel('$r$','interpreter','latex');
        ylabel('$z$','interpreter','latex');
     if rtype(1)~=0
            plot(Px(1),Py(1),'*','Color',red);
     end
     % Legend
        h = zeros(8, 1);
        h(1) = plot(NaN,NaN,'Color',blue,'LineWidth',1.5);
        h(2) = plot(NaN,NaN,'Color',red,'LineWidth',1.5);
        h(3) = plot(NaN,NaN,'--k','LineWidth',1);
        h(4) = plot(NaN,NaN,':k','LineWidth',1);
        h(5) = plot(NaN,NaN,'-.k','LineWidth',1.2);
        h(6) = plot(NaN,NaN,'*','Color',red);
        h(7) = plot(NaN,NaN,'ok','MarkerSize',marksize,'LineWidth',1);
        h(8) = plot(NaN,NaN,'xk','MarkerSize',marksize,'LineWidth',1);

        
        leg=legend(h, 'Refract. surf.','Reflect. surf.','$L_1$','$L_2$','Opt. ax.','Source','$_{j}Q$','$_{j}T$','interpreter','latex','FontSize',8,'Location', 'Best');
        leg.ItemTokenSize = [15,20];

end
end