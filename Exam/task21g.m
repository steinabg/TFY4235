% Task 2.1g. Write the code for solving (1) numerically
% 
close all
clear all
clc
tic; % Initialize timer of program

jmax = 201; % Number of discrete spatial points
tend = 5e5; % SET ENDTIME OF SIMULATION HERE! 5e5 virker som en god tid
xstart = 0; % Left boundary. a.
xend = 1;   % Right boundary. b.
dx = (xend-xstart)/(jmax-1); % Calculate spatial stepsize dx
dt = 1e-1; % Set timestep size
r = dt/dx^2;    % Calculate value of much used constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set diffusivity constants
Da = 4e-7; % Diffusivity
Db = 2/3*Da;
Dc = 8/15*Da;


% Comment on Diffusivity values:
% The size of Dc relative to (Da,Db) decides whether or not sheets of s
% will form. If Dc is very large in relation to (Da,Db), the gas c, will
% diffuse in a manner that makes a high enough density of c for s to occur,
% hard. If Dc is very small (i.e. Dc -> 0), the gas c will have "a hard
% time diffusing", so a large concentration of c can easily occur. This
% means that s will be generated. When s is generated, the density of c
% decreases in that spot. Scaling all values Da, Db, Dc, results in the
% same formation of s sheets.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate constants used in FDM
ra = Da*r;
rb = Db*r;
rc = Dc*r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Nucleation Threshold %%%%%%%%%
c0 = 3e-2; % Given in task
% c0 = 3e-1;
% c0 = 3e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Control variable diffusivity
s0 = 1e0; % Variable for varying DBarG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare and initialize constants for the PDEs
R = 1; % Reaction speed
N1 = R/10; % 
N2 = N1;
gc = 1;
L = xend-xstart;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xspace = linspace(xstart,xend,jmax); % Vector for spatial coordinates

% Vectors for storing the solutions
a = zeros(length(xspace),1);
b = zeros(length(xspace),1);
c = zeros(length(xspace),1);
s = zeros(length(xspace),1);

% Vectors for storing the spatially varying diffusivity values
DBarA = zeros(length(xspace),1);
DBarB = DBarA;
DBarC = DBarA;


% Set Boundary conditions for the reservoirs:
a(1) = 1;
a(jmax) = 0;
b(1) = 0;
b(jmax) = a(1)*10;


figure(1); % Prepare for plotting!
for t = 1:length(dt:dt:tend)
    
    % Vectors for storing solution vectors for the previous timestep for
    % use in the FTCS scheme.
    tempa = a;
    tempb = b;
    tempc = c;
    temps = s;
    
    % Calculate the values of the spatially varying diffusivity
    DBarA = Da./(1+s/s0);
    DBarB = Db./(1+s/s0);
    DBarC = Dc./(1+s/s0);
    
    for jj = 2:jmax-1
        a(jj) = tempa(jj) + r/4*( (DBarA(jj+1)-DBarA(jj-1))* ...
            (tempa(jj+1)-tempa(jj-1)) ) + DBarA(jj)*r*(tempa(jj+1)-2* ...
            tempa(jj)+tempa(jj-1)) - R*dt*tempa(jj)*tempb(jj);
        
        b(jj) = tempb(jj) + r/4*( (DBarB(jj+1)-DBarB(jj-1))* ...
            (tempb(jj+1)-tempb(jj-1)) ) + DBarB(jj)*r*(tempb(jj+1)-2*...
            tempb(jj)+tempb(jj-1)) - R*dt*tempa(jj)*tempb(jj);
        
        c(jj) = tempc(jj) + r/4*( (DBarC(jj+1)-DBarC(jj-1))*(tempc(jj+1)-tempc(jj-1)) ) + ...
            DBarC(jj)*r*(tempc(jj+1)-2*tempc(jj)+tempc(jj-1)) + R*dt*tempa(jj)*tempb(jj)- ...
            N1*dt*(tempc(jj)>=c0)*(tempc(jj))^2 - N2*dt*tempc(jj)*temps(jj);
              
        s(jj) = temps(jj) + N1*dt*(tempc(jj)>=c0)*(tempc(jj)^2) + ...
            N2*dt*tempc(jj)*temps(jj) ;
        
    end
     

   if(~mod(t,10000)) % Limit live-plotting
    figure(1);
    subplot(2,2,1:2);
    plot(xspace,a,xspace,b)
    title(sprintf('jmax = %d, dt= %0.2f, t = %d, c0 = %0.2f, ', jmax, dt, t*dt, c0));
    subplot(2,2,3)
    plot(xspace,c,'k')
    title('Concentration of C')
    subplot(2,2,4);
    plot(xspace, s,'c')
    title('Concentration of S')
    drawnow();
   end
    
end

% Plot the final result
h = figure(1);
subplot(2,2,1:2);
plot(xspace,a,xspace,b)
title(sprintf('jmax = %d, dt= %0.2f, t = %d, c0 = %0.2f, s_0 = %1.0e ', jmax, dt, t*dt, c0,s0));
subplot(2,2,3)
plot(xspace,c,'k')
title('Concentration of C')
subplot(2,2,4);
plot(xspace, s,'c')
title('Concentration of S')

% Save the figure and the s-vector
figureName = sprintf('NonConstD_jmax%d_dt%0.2f_t%1.0e_c0__%0.2f_Da%1.2e_Db%1.2e_Dc%1.2e_s0_%1.0e.fig', ...
    jmax, dt, t*dt, c0, Da, Db, Dc, s0);

savefig(h,figureName);


savefile = sprintf('sVector_NonConstD_jmax%d_dt%0.2f_t%1.0e_c0__%0.2f_Da%1.2e_Db%1.2e_Dc%1.2e_s0_%1.0e.mat', ...
    jmax, dt, t*dt, c0, Da, Db, Dc, s0);
save(savefile,'s')
toc
% Make a noise when the program is done
load handel
sound(y,4*Fs)