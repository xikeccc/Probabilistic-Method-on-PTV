clear all;
close all;
clc;

number=100; % number of particles
dt=0.01; % time interval

pt1=[0.503840081846848,-0.859744640730556;0.612809588835559,2.26013019411244;3.81942224223188,4.08681510086583;-0.468110830507074,3.42939733708581;0.202075096548432,0.257282784769860;4.45389346568052,2.29755538415112;-3.57208908017983,-0.575141588295374;4.96605274961200,4.11920725942129;-2.37994494586653,4.49506692380046;2.69538995191800,3.70640722753756;1.72016461201798,-0.756426627319049;2.34689519028761,0.785070081934006;-0.483009579474690,-1.92591042313907;1.55669463201210,3.39388342698170;3.15649521911404,3.00339412296431;-3.43794395443823,4.22067689193302;-3.30519671089154,1.00130056870060;4.42645553571024,1.18917967099474;1.83627042067631,-2.85751594517903;-1.26861293383807,4.26807599912338;-0.639968957914228,0.174892065335783;-2.54578762814660,-2.86135102840741;-1.61361010140783,4.59888561036787;-1.22444535774946,2.90105790568273;-1.26572889431406,3.93937976588424;-2.56972215492584,-1.77881554406092;-0.306247425866087,2.48267137539641;-2.05478651283154,1.37601111592475;3.78423259823174,-2.47621996303413;-3.29442814186893,-0.735127406451410;4.10933423952435,-1.93164277952917;-0.610069342843515,2.43632707748010;-3.40909526959466,-1.82614696263500;-0.540619952033675,4.02610710815491;2.05033998673895,3.95467827408045;3.22868758402022,-0.569403480140583;0.834189061056900,0.961558573103663;2.01564469269653,2.76241448400299;4.86371086505587,3.00734866110285;-3.92193094696757,1.68003864046014;-1.33095740991844,1.70595075424888;-3.49978867571890,-1.35487121339262;-2.78200620124985,0.552309844507730;2.57161572540409,2.21810872978328;2.12218915092342,3.77236621361304;3.67116623298865,2.22802832360045;1.59958554812281,-3.62913528102659;-3.94402384261221,2.89092882102570;4.05634301854336,-0.143623092271739;3.15250063700247,-0.597566450186556;-1.98037893323442,-2.68198088250327;0.435175545696952,3.60863543063499;4.83222147528282,-1.08980477330552;2.61739017145408,-1.09090181219912;3.52012941531074,-0.408405591092562;-0.136131777066389,-0.667428592664449;0.0976979181614431,0.853063629209919;1.90805220318677,1.44239789304117;-1.89198330586324,-0.0956445217820563;2.51699675809695,-0.966820594047873;0.143156022083576,0.532426482500133;1.55937057240300,1.71649734650438;-1.99542037605268,4.17930184382179;-2.23331800137851,2.33653292580106;0.848709226458282,-0.812287051424733;-0.0831787297462620,3.32192718318292;3.98696827478366,-2.59614328876979;3.50513310179882,-3.45143370014304;2.27142162441752,-3.95126140721773;-2.89924948807876,-2.44726786682106;-3.49215116917046,-1.72518859516250;1.58560912570188,-1.75849825826589;-1.23711290408926,-3.75685482059345;4.08296264911054,0.154159449170074;4.66159619308272,-3.04358363899198;-1.48302098529379,-2.06433862448340;3.17104801752545,1.81871443720824;4.93855786433184,3.72826184782283;1.59048317714257,4.17581172827553;3.44063468076084,1.36037097704245;4.94191893031128,4.18878997264161;0.655913820257141,0.00119839684004341;2.45194570926039,0.316419513732543;1.83969742071967,-1.30038301360294;-3.46737649752565,-0.374744819820960;0.553887065791275,0.543062175343852;1.68006553008336,-3.56096279661233;0.367189905317367,3.28742726748254;-0.760709393806455,-3.49834089325331;4.57892349245909,-0.238453814330353;3.86688705467251,3.76240804878871;3.40677676021523,-0.423944098503443;-0.887384858974953,1.74766283764172;1.44384583672696,3.64553450598132;3.30018440121390,3.12321951832451;4.40138685381449,4.50439786009277;2.83336356345313,-2.65273868727837;-2.59637133722639,-1.90785231519490;1.39017593813061,4.14784946803252;-3.63955110662106,1.19816970106640];
% this is original coordinates data used for report, I keep it the same for 
% the purpose of controling variables when testing the algorithm in different time
% interval dt, but it can change to random in 'Velocity calculation' section.

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
grid_size = 5; % n*n axis for velocity plotting
grid_space2 = .5; % resolution
reso = -1*grid_size:grid_space2:grid_size;
[x1,y1] = meshgrid(reso,reso);

%% Velocity calculation

%pt1 = randi([-grid_size+1,grid_size-1],number,2)+rand(number,2); % generate random coodinate matrix
%pt1=importdata('pt1_100.mat');
pt2 = zeros(max(size(pt1)),2);
z1 = pt1(:,1)+1i*pt1(:,2);
zc1=zc*ones(number,1);
z1_new = z1-zc1;

dF_uniform_dz1 = U * exp(-1i*a_rad); % velocity from uniform flow at angle of attack
dF_source_dz1 = (qs/(2.0*pi)) ./ z1_new; % velocity from source/sink
dF_doublet_dz1 = - U * (R^2 * exp(1i*a_rad)) ./ (z1_new.^2); % velocity from doublet (at angle of attack)
dF_vortex_dz1 = -1i * circ ./ z1_new; % velocity from free vortex

vel1 = (dF_uniform_dz1 + dF_source_dz1 + dF_doublet_dz1 + dF_vortex_dz1); % linear superposition of basic potential flows

pt2(:,1)=pt1(:,1)+real(vel1)*dt;
pt2(:,2)=pt1(:,2)+imag(vel1)*dt;

% Plot
figure;
plot(pt1(:,1),pt1(:,2),'o');
hold on
plot(pt2(:,1),pt2(:,2),'s');
title('Particle Distribution')

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

%% Computing affinity matrix

pt1; %original position matirx
pt2; %final position matrix
h=0.7; %gaussian kernal bandwidth ?
affinity_s=[]; %single match matrix

for i=1:length(pt1)
    for j=1:length(pt2)
        affinity_s=[affinity_s;exp(-((pt2(j,1)-pt1(i,1))^2 + (pt2(j,2)-pt1(i,2))^2)/h)]; %does norm2 need sqrt?
    end
end

affinity_p=affinity_s*transpose(affinity_s); %pair match matrix

%% Soft thresholding

A=affinity_p;
u1=0.01; % u
u2=0.5; % threshold 
countB=0; % iteration number

sum_matrix=zeros();
A_1_1=zeros();
A_1_5=zeros();
A_1_102=zeros();
A_1_104=zeros();
A_1_203=zeros();
A_1_215=zeros();

for ix = 1:number^2  % go through all row
    while sum(A(ix,:))>(number-u2)
        for iy = 1:number^2 % go through all row
            for iz = 1:number^2 % go through all column
            A(iy,iz) = max( 0,abs(A(iy,iz))-u1/2 ); 
            end
        end

        countB=countB+1;
      
        %sum_matrix = [ sum_matrix ; sum( A(:) ) ];
        %A_1_1 = [ A_1_1; A(1,1)];
        %A_1_5 = [ A_1_5; A(1,5)];
        %A_1_102 = [ A_1_102; A(1,102)];
        %A_1_104 = [ A_1_104; A(1,104)];
        %A_1_203 = [ A_1_203; A(1,203)];
        %A_1_215 = [ A_1_215; A(1,215)];
    end 
end

%plot
%figure;
%plot([1:10:countB],A_1_1(2:10:max(size(A_1_1)),:),'s-');
%hold on
%plot([1:10:countB],A_1_5(2:10:max(size(A_1_5)),:),'o-');
%hold on
%plot([1:10:countB],A_1_102(2:10:max(size(A_1_102)),:),'d-');
%hold on
%plot([1:10:countB],A_1_104(2:10:max(size(A_1_104)),:),'<-');
%hold on
%plot([1:10:countB],A_1_203(2:10:max(size(A_1_203)),:),'>-');
%hold on
%plot([1:10:countB],A_1_215(2:10:max(size(A_1_215)),:),'^-');
%legend('1-1','1-2','2-2','2-3','3-3','3-4');
%title('Pairwise Probability');
%xlabel('Iteration');
%ylabel('Probability');
%grid on
%set(gca,'yscale','log')          
                
 %% Egozi Algorithm 1

B=A; %given affinity matrix 
%B=affinity_p;
IterNo=15;      %iteration number
[r,c]=size(B);
delta=10^-3;     %threshold
p=transpose(1/c*ones(1,c)); %initialize p     
p_matrix=(p'); %threshold storage for step 7

result_1_1=zeros();
result_1_2=zeros();
result_2_2=zeros(); 
result_2_3=zeros();
result_3_1=zeros(); 
result_3_3=zeros();% storage matrix for plotting

for t = 0:(IterNo-1)                         
    q = B*p_matrix(t+1,:)' ;                                          
    %p = normalize(q);
    p=0.5*normalize(q)+0.5*p;
    p_matrix = [p_matrix;p']; %store p
   
    B = B .* ( p_matrix(t+2,1 ) / p_matrix( t+1,1) );

    %if norm ( p_matrix(t+2,1 ) - p_matrix( t+1,1 ),2 ) / (r*c) < delta
   % break
    
    %end 
   
    result=B*p_matrix(t+1,:)';
    result_1_1=[result_1_1;result(1)/10^8];
    result_1_2=[result_1_2;result(5)/10^8];
    result_2_2=[result_2_2;result(102)/10^8];
    result_2_3=[result_2_3;result(104)/10^8];
    result_3_1=[result_3_1;result(203)/10^8];
    result_3_3=[result_3_3;result(215)/10^8];

end 

%figure;
%plot(1:max(size(p)),p(:));

figure;
plot((1:t),result_1_1(3:max(size(result_1_1)),:),'s-')
hold on
plot((1:t),result_1_2(3:max(size(result_1_2)),:),'o-')
hold on
plot((1:t),result_2_2(3:max(size(result_2_2)),:),'d-')
hold on
plot((1:t),result_2_3(3:max(size(result_2_3)),:),'<-')
hold on
plot((1:t),result_3_1(3:max(size(result_3_1)),:),'>-')
hold on
plot((1:t+1),result_3_3(2:max(size(result_3_3)),:),'^-')
legend('1-1','1-2','2-2','2-3','3-1','3-3')
xlabel('iteration');
ylabel('matching index');
title('Matching index')

%C=B*affinity_s; % adjacency matrix
C=B*p;
C_re=reshape(C,number,number)'; % reshape C to a n*n matrix
[max_a,index]=max(C_re,[],2); % find location of the most probable matching
%[max_a,index]=max(abs(C_re),[],2);

figure;
for i=1:number
    %line( [pt1(i,1) , pt2(index(i),1)] , [pt1(i,2),pt2(index(i),2)],'color','red','linestyle','-');
    quiver( pt1(i,1), pt1(i,2), pt2(index(i),1)-pt1(i,1), pt2(index(i),2)-pt1(i,2),'MaxHeadSize',3, 'linewidth',1,'color','r')
    hold on
    
end


quiver(real(z_bar),imag(z_bar),real(vel_normalised),imag(vel_normalised),0.5); %plotting the velocity direction as arrows
legend('velocity direction')
hold on

plot(pt1(:,1),pt1(:,2),'s','DisplayName','old position');
hold on
plot(pt2(:,1),pt2(:,2),'+','DisplayName','new position');
title('Velocity field')
axis square;