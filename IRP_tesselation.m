clear all;
close all;
clc;

number=20; % number of particles
dt=0.01; % time interval

grid_size = 5; % n*n axis for velocity plotting
grid_space2 = .5; % resolution

%% Fluid parameters
%Uniform Flow Characteristics
U = 5.0; %flow speed
a_deg = 30.0; %angle of attack (in degrees)
a_rad = -a_deg*pi/180.0; %angle of attack (in radians)


%Source/Sink Characteristics
qs = 0.0; %strength of source

%Doublet Characteristics
R = 2.0; %radius of doublet

%Vortex Characteristics
%  circ = 0.0; %circulation
circ = -4.0; %circulation
%  circ = -10.0; %circulation

%Locations
 % zc = 0.0 + 0.0 * 1i; %centre of the source/doublet/vortex
zc = -0.5 + 0.5 * 1i; %centre of the source/doublet/vortex
%  zc = -1.0 + 1.0 * 1i; %centre of the source/doublet/vortex

%Grid Sizing
reso = -1*grid_size:grid_space2:grid_size;
[x1,y1] = meshgrid(reso,reso);

%% Velocity calculation

pt1 = randi([-grid_size+1,grid_size-1],number,2)+rand(number,2); % generate random coodinate matrix
%pt1=importdata('pt1_50.mat');
pt2 = zeros(max(size(pt1)),2);
z1 = pt1(:,1)+1i*pt1(:,2);
zc1=zc*ones(number,1);
z1_new = z1-zc1;

dF_uniform_dz1 = U * exp(-1i*a_rad); % velocity from uniform flow at angle of attack
dF_source_dz1 = (qs/(2.0*pi)) ./ z1_new; % velocity from source/sink
dF_doublet_dz1 = - U * (R^2 * exp(1i*a_rad)) ./ (z1_new.^2); % velocity from doublet (at angle of attack)
dF_vortex_dz1 = -1i * circ ./ z1_new; % velocity from free vortex

vel1 = (dF_uniform_dz1 + dF_source_dz1 + dF_doublet_dz1 + dF_vortex_dz1); % linear superposition of basic potential flows

dim1 = max(size(pt1));
for temp1=1:dim1
    for temp2=1:2
        if (((pt1(temp1,temp2)-real(zc))^2 + (pt1(temp1,temp2)-imag(zc))^2) < R^2) 
            vel1(temp1) = 0.0;
        end
    end
end

pt2(:,1)=pt1(:,1)+real(vel1)*dt;
pt2(:,2)=pt1(:,2)+imag(vel1)*dt;

%% Joukowski transform

z2 = x1 + 1i*y1; % complex coordinates (coarse grid)
z2_new = z2-zc; % coordinates relative to centre of source/doublet/vortex (coarse grid)

dF_uniform_dz = U * exp(-1i*a_rad); % velocity from uniform flow at angle of attack
dF_source_dz = (qs/(2.0*pi)) ./ z2_new; % velocity from source/sink
dF_doublet_dz = - U * (R^2 * exp(1i*a_rad)) ./ (z2_new.^2); % velocity from doublet (at angle of attack)
dF_vortex_dz = -1i * circ ./ z2_new; % velocity from free vortex
vel = (dF_uniform_dz + dF_source_dz + dF_doublet_dz + dF_vortex_dz); % linear superposition of basic potential flows

a = 1.5; % coefficient in the Joukowski conformal map
z_bar = z2 + (a^2) ./ z2; % the Joukowski conformal map
vel_m = vel ./ (1.0 - ((a^2)./(z2.^2))); % the velocity in the mapped domain

vel_normalised = vel_m./abs(vel_m); % normalising the vector length (easier to visualise)

dim = max(size(reso));
for temp1=1:dim
    for temp2=1:dim
        if (((x1(temp1,temp2)-real(zc))^2 + (y1(temp1,temp2)-imag(zc))^2) < R^2) % removing points inside the body, hence visualising only the flow around the doublet
            vel_normalised(temp1,temp2) = 0.0;
        end
    end
end
figure;
quiver(real(z_bar),imag(z_bar),real(vel_normalised),imag(vel_normalised),0.5); %plotting the velocity direction as arrows
title('Velocity field')
axis square;
    
%% Delaunay
pt1;
pt2;

DT1 = delaunayTriangulation(pt1);
DT2 = delaunayTriangulation(pt2);

hold on
triplot(DT1,'-r');
triplot(DT2,'-o');

NTriangles1 = size(DT1.ConnectivityList,1);
Areas1 = zeros(NTriangles1,1);
NTriangles2 = size(DT2.ConnectivityList,1);
Areas2 = zeros(NTriangles2,1);
for i = 1:NTriangles1
    PointIndexes = DT1.ConnectivityList(i,:);
    Areas1(i) = polyarea(DT1.Points(PointIndexes,1),DT1.Points(PointIndexes,2));
end
for i = 1:NTriangles2
    PointIndexes = DT2.ConnectivityList(i,:);
    Areas2(i) = polyarea(DT2.Points(PointIndexes,1),DT2.Points(PointIndexes,2));
end

%% Soft thresholding & Egozi
h=0.7;
affinity_s=[]; %single match matrix
for i=1:length(Areas1)
    for j=1:length(Areas2)
        affinity_s=[affinity_s;exp(-abs(Areas2(j)-Areas1(i))/h)]; %does norm2 need sqrt?
    end
end

affinity_p=affinity_s*transpose(affinity_s);
B=affinity_p;
IterNo=15;      %iteration number
[r,c]=size(B);
delta=10^-3;     %threshold
p=transpose(1/c*ones(1,c)); %initialize p     
p_matrix=(p'); %threshold storage for step 7
for t = 0:(IterNo-1)                         
    q = B*p_matrix(t+1,:)' ;                                          
    %p = normalize(q);
    p=0.5*normalize(q)+0.5*p;
    p_matrix = [p_matrix;p']; %store p
    B = B .* ( p_matrix(t+2,1 ) / p_matrix( t+1,1) ); 
end 

%% Matching result
C=B*p;
C_re=reshape(affinity_s,max(size(Areas1)),max(size(Areas2)))';  
[max_a,index]=max(C_re,[],2); % find location of the most probable matching

p_index=[];

for i =1: max(size(C_re))
    Position_index1=DT1.ConnectivityList(index(i),:);
    Position_index2=DT2.ConnectivityList(index(i),:);
    p_index=[ p_index; Position_index1 , Position_index2];
end

temp3=size(p_index);
figure;
for i=1: 3
   for j=1:temp3(1)
    x_start=pt1( p_index(j,i),1 );
    x_finish= pt2(p_index(j,i+3),1);
    y_start=pt1( p_index(j,i),2 );
    y_finish=pt2(p_index(j,i+3),2);
    line([x_start,x_finish],[y_start,y_finish],'color','red','linestyle','-');
   end
end
hold on
quiver(real(z_bar),imag(z_bar),real(vel_normalised),imag(vel_normalised),0.5); %plotting the velocity direction as arrows
legend('velocity direction')
hold on

plot(pt1(:,1),pt1(:,2),'s','DisplayName','old position');
hold on
plot(pt2(:,1),pt2(:,2),'+','DisplayName','new position');
title('Velocity field')
xlabel('x');
ylabel('y');
axis square;