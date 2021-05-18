% confirm that rotation works fine, N=2, N=3 and N=100
% rotate in a cartesian plane

clear all
N = 100; 
opt_slerp =0;  % 1 for slerp and ~=1 for rotors
opt_map =1;    % 1 for mapping between initial and final vecs. theta deteremined by initial angle
col = [204 194 193]/255; 

% initial plots
if N <=3; plt_opt = 1; else; plt_opt =0; end

% Create the matrix whose eigenvector are to be rotated
 Mat0 = rand(N,N);  
 Mat0 = 0.5*(Mat0 + Mat0');
 
[eigV0,DV0] = eig(Mat0);
[lam1,I]= sort(diag(DV0));
eigV0 = eigV0(:,I);
V1  = eigV0(:,N);     V2 = eigV0(:,N-1);

% Specify axes of rotation-plane
if N==3
Ax1 = [1; zeros(N-1,1)]; 
Ax2 = [0; 1; zeros(N-2,1)];
else 
Ax1 = V1; 
Ax2 = V2;
end

%%%  calculate the components of bivector
B = zeros(N,N); 
for i = 1:N
    for j = i+1:N
        B(i,j) = Ax1(i)*Ax2(j) - Ax2(i)*Ax1(j); 
    end
end

% make the skew-symmetric matrix
B = B - B';

% define rotation matrix
if opt_map ==1
    th = acos(Ax2'*Ax1);
else
    th = pi/2;
end
R = expm(th*B);

sv = linspace(0,1,51); 


% generate an plot initial axes
ax_mat = eye(N); 
origin = zeros(N,1); 

if plt_opt ==1
figure(1); subplot(1,2,1)
    if N==2
        quiver(origin,origin,ax_mat(:,1), ax_mat(:,2),1,'k','LineWidth',2)
        hold on  
        quiver(origin(1),origin(2),Ax1(1), Ax1(2),1,'b','LineWidth',2)
        if opt_slerp ==1
            quiver(origin(1),origin(2),V2(1), V2(2),1,'b','LineWidth',2)  
        end
        xlabel('x')
        ylabel('y')
        set(gca,'fontsize',14)
    else
        quiver3(origin',origin',origin',ax_mat(1,:), ax_mat(2,:),ax_mat(3,:),1.5,'k','LineWidth',2)
        hold on  
        % plot initial vector that is to be rotated
        quiver3(0,0,0, V1(1), V1(2),V1(3),1.2,'b','LineWidth',2)
%         if opt_map ==1
%             quiver3(0,0,0, Ax2(1), Ax2(2),Ax2(3),1.5,'b','LineWidth',2)
%         end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        set(gca,'fontsize',14)
    %   set view angle
	th_view = deg2rad(129.3); 	ph_view = deg2rad(34.37+180); 
	view([cos(th_view)*cos(ph_view), cos(th_view)*sin(ph_view), sin(th_view)]); 
    hold on 
    end
end
%%%
achk = zeros(N,size(sv,2));

for js = 1:size(sv,2)
    s = sv(js);
    
        % using rotation matrix
        Rs = expm(-s*th*B);        
        MatR = Rs*Mat0*Rs';
        
        [eigV0, DV0] = eig(Mat0);
        [lam0,I0] = sort(diag(DV0));
        eigV0 = eigV0(:,I0);
        
        [eigV1,DV1] = eig(MatR); 
        [lam2,I2]= sort(diag(DV1));
        eigV1 = eigV1(:,I2);
%       norm(abs([eigV1 Rs*eigV0]))
        
        V0 = eigV0(:,N);
        Vt = Rs*V0;  
        norm_chk(js) = norm(Vt);
        detmnt(js) = det(Rs);
        alignmnt(js) = Vt'*V0; % calculate alignment with components parallel to plane 
    
        if plt_opt==1 
        figure(1)
        subplot(1,2,1)
        if N==2
            quiver(origin,origin,ax_mat(:,1), ax_mat(:,2),1,'k','LineWidth',2)
            hold on  
            quiver(origin(1),origin(2), Vt(1), Vt(2),1,'LineWidth',2)
            drawnow
            axis equal
        elseif N==3
             quiver3(0,0,0,Vt(1), Vt(2), Vt(3),1,'LineWidth',2, 'color',col)
%             quiver3(0,0,0,eigV1(1,3), eigV1(2,3), eigV1(3,3),1,'LineWidth',2)
             hold on     
             % set view angle
             th_view = deg2rad(129.3); 	ph_view = deg2rad(34.37+180); 
             
             drawnow
             axis equal
        end
        subplot(1,2,2)
        end
        
end

subplot(1,2,1)
if N==3
quiver3(0,0,0, Vt(1), Vt(2),Vt(3),1.2,'r','LineWidth',2)
end

subplot(1,2,2)
plot(sv,alignmnt)
hold on
axis square








