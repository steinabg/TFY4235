% Task 2.1b. Write the code for solving (1) numerically
% 
close all
clear all
clc
tic; % time the program
jmax = 401; % Number of discrete points
tend = 5e5; % Endtime. 5e5 virker som en god tid
xstart = 0; % Left boundary. a.
xend = 1; % Right boundary. b.
dx = (xend-xstart)/(jmax-1); % spatial stepsize
dt = 1e-1; % timestep
r = dt/dx^2; % define much used constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
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

ra = Da*r;
rb = Db*r;
rc = Dc*r;

c0 = 3e-2; % Given in task
% c0 = 1.5e-2;
% c0 = 6e-2;

R = 1;
N1 = R/10;
N2 = N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



xspace = linspace(xstart,xend,jmax); % For plotting

% Initialize solution vectors
a = zeros(length(xspace),1);
b = zeros(length(xspace),1);
c = zeros(length(xspace),1);
s = zeros(length(xspace),1);

% set ICs:
a(1) = 1;
a(jmax) = 0;
b(1) = 0;
b(jmax) = a(1)*10;




figure(1);

for t = 1:length(dt:dt:tend)
% for t = 1:noOfTimeSteps
    
    % Solve for a and b with FTCS scheme
    tempa = a;
    tempb = b;
    tempc = c;
    temps = s;
    
    for jj = 2:jmax-1
        a(jj) = tempa(jj) + ra*(tempa(jj+1)-2*tempa(jj)+tempa(jj-1)) ...
            - tempa(jj)*tempb(jj)*dt*R;
        
        b(jj) = tempb(jj) + rb*(tempb(jj+1)-2*tempb(jj)+tempb(jj-1)) ...
            - tempa(jj)*tempb(jj)*dt*R;
        
        c(jj) = tempc(jj) + rc*(tempc(jj+1)-2*tempc(jj)+tempc(jj-1)) + ...
            tempa(jj)*tempb(jj)*dt*R - (tempc(jj)>=c0)* ...
            N1*dt*(tempc(jj))^2 - N2*tempc(jj)*temps(jj)*dt;
        
        s(jj) = temps(jj) + N1*dt*(tempc(jj)>=c0)*(tempc(jj)^2) + ...
            N2*dt*tempc(jj)*temps(jj) ;
    end
     

   if(~mod(t,10000)) %Limit live-plotting
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
title(sprintf('jmax = %d, dt= %0.2f, t = %d, c0 = %0.2f, ', jmax, dt, t*dt, c0));
subplot(2,2,3)
plot(xspace,c,'k')
title('Concentration of C')
subplot(2,2,4);
plot(xspace, s,'c')
title('Concentration of S')

% Save the figure and the s-vector
figureName = sprintf('jmax%d_dt%0.2f_t%1.0e_c0__%0.2f_Da%1.2e_Db%1.2e_Dc%1.2e.fig', ...
    jmax, dt, t*dt, c0, Da, Db, Dc);

savefig(h,figureName);


savefile = sprintf('sVector_jmax%d_dt%0.2f_t%1.0e_c0__%0.2f_Da%1.2e_Db%1.2e_Dc%1.2e.mat', ...
    jmax, dt, t*dt, c0, Da, Db, Dc);
save(savefile,'s')
toc
% Make a noise when the program is done
load handel
sound(y,Fs)