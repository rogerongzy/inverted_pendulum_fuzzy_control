function test(varargin)

if mod(length(varargin),2)==1
    error('Error input. input should be (''t_0'',10)')
end
%Generate default values
t_step=0.001;
L=floor(6./t_step);
g0=1;
g1=1;
h=1;
t_0=0;
theta_0=5*pi/180;
dtheta_0=pi/180;
ddtheta_0=0;
x_0=0;
dx_0=-0.5;
ddx_0=0;
F_0=0;
reference_theta=0;
reference_dtheta=0;
rulebase=[6,6,6,5,4,3; % 1
          6,6,5,4,3,3;
          6,5,4,3,3,2;
          5,4,4,3,2,1;
          4,4,3,2,1,1;
          4,3,2,1,1,1];
centerpoint=[-5*pi/12 -pi/4 -pi/12 pi/12 pi/4 5*pi/12; -5*pi/36 -3*pi/36 -pi/36 pi/36 3*pi/36 5*pi/36; -50/6 -30/6 -10/6 10/6 30/6 50/6]; % 1
width=[pi/3 pi/3 pi/3 pi/3 pi/3 pi/3; pi/9 pi/9 pi/9 pi/9 pi/9 pi/9; 40/6 40/6 40/6 40/6 40/6 40/6]; % 1
% rulebase=[0,0,5,0,0; % 2
%           0,4,0,3,0;
%           0,0,3,0,0;
%           0,3,0,2,0;
%           0,0,1,0,0];
% centerpoint=[-pi/3 -pi/6 0 pi/6 pi/3; -pi/9 -pi/18 0 pi/18 pi/9; -20/3 -10/3 0 10/3 20/3]; % 2
% width=[pi/3 pi/3 pi/3 pi/3 pi/3; pi/9 pi/9 pi/9 pi/9 pi/9; 40/6 40/6 40/6 40/6 40/6]; % 2
FigureNumber=1;
functiontype='triangle';
COGtype='min';
ForceInput=0;
StartTime=0;
EndTime=0;
% Get input
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'t_step')
        t_step=varargin{i+1};
    elseif strcmpi(varargin{i},'L')
        L=varargin{i+1};
    elseif strcmpi(varargin{i},'g0')
        g0=varargin{i+1};
    elseif strcmpi(varargin{i},'g1')
        g1=varargin{i+1};
    elseif strcmpi(varargin{i},'h')
        h=varargin{i+1};
    elseif strcmpi(varargin{i},'t_0')
        t_0=varargin{i+1};
    elseif strcmpi(varargin{i},'theta')
        theta_0=varargin{i+1}(1);
        dtheta_0=varargin{i+1}(2);
        ddtheta_0=varargin{i+1}(3);
        reference_theta=varargin{i+1}(4);
        reference_dtheta=varargin{i+1}(5);
    elseif strcmpi(varargin{i},'x')
        x_0=varargin{i+1}(1);
        dx_0=varargin{i+1}(2);
        ddx_0=varargin{i+1}(3);
    elseif strcmpi(varargin{i},'F_0')
        F_0=varargin{i+1};
    elseif strcmpi(varargin{i},'rulebase')
        rulebase=varargin{i+1};
    elseif strcmpi(varargin{i},'centerpoint')
        centerpoint=varargin{i+1};
    elseif strcmpi(varargin{i},'width')
        width=varargin{i+1};
    elseif strcmpi(varargin{i},'functiontype')
        functiontype=varargin{i+1};
    elseif strcmpi(varargin{i},'COGtype')
        COGtype=varargin{i+1};
    elseif strcmpi(varargin{i},'ForceInput')
        ForceInput=varargin{i+1};
    elseif strcmpi(varargin{i},'StartTime')
        StartTime=varargin{i+1};
    elseif strcmpi(varargin{i},'EndTime')
        EndTime=varargin{i+1};
    elseif strcmpi(varargin{i},'FigureNumber')
        FigureNumber=varargin{i+1};
    else
        error(['Unknown inputs: ' varargin{i}]);
    end
end
% Initial variables
t=zeros(1,L);
theta=zeros(1,L);
dtheta=zeros(1,L);
ddtheta=zeros(1,L);
x=zeros(1,L);
dx=zeros(1,L);
ddx=zeros(1,L);
F=zeros(1,L);
inputF=zeros(1,L);
t(1)=t_0;
theta(1)=theta_0;
dtheta(1)=dtheta_0;
ddtheta(1)=ddtheta_0;
x(1)=x_0;
dx(1)=dx_0;
ddx(1)=ddx_0;
F(1)=F_0;
inputF(1)=F_0;
if ForceInput==0 || StartTime==EndTime
    isForceInput=0;
else
    isForceInput=1;
end
% begin to test
for i=2:L
    % Calculate next input force according to previous situation
    inputF(i)=FuzzyController(reference_theta-theta(i-1),...
        reference_dtheta-dtheta(i-1),g0,g1,h,rulebase,...
        centerpoint,width,functiontype,COGtype);
    % Calculate next situation according to next input force
    if ~isForceInput
        [t(i),theta(i),dtheta(i),ddtheta(i),x(i),dx(i),ddx(i),F(i)]=...
            InvertedPendulum(t(i-1),theta(i-1),dtheta(i-1),ddtheta(i-1),...
            x(i-1),dx(i-1),ddx(i-1),F(i-1),inputF(i),t_step);
    else
        [t(i),theta(i),dtheta(i),ddtheta(i),x(i),dx(i),ddx(i),F(i)]=...
            InvertedPendulum(t(i-1),theta(i-1),dtheta(i-1),ddtheta(i-1),...
            x(i-1),dx(i-1),ddx(i-1),F(i-1),inputF(i),t_step,...
            ForceInput,StartTime,EndTime);
    end
end
% plot results
fontsize=10;
linewidth=2;
marksize=10;
% figure
figure(FigureNumber(1));
subplot(5,1,1)
% angle
plot(t,theta,'LineWidth',linewidth,'MarkerSize',marksize);
axis([min(t) max(t) min(theta) max(theta)])
grid on;
xlabel('Time (s)','FontSize',fontsize);
ylabel('Angle (rad)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
% angular velocity
subplot(5,1,2)
plot(t,dtheta,'LineWidth',linewidth,'MarkerSize',marksize);
axis([min(t) max(t) min(dtheta) max(dtheta)])
grid on;
xlabel('Time (s)','FontSize',fontsize);
ylabel('Angular Velocity (rad/s)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
% position
subplot(5,1,3)
plot(t,x,'LineWidth',linewidth,'MarkerSize',marksize);
axis([min(t) max(t) min([min(x) min(dx)]) max([max(x) max(dx)])])
grid on;
xlabel('Time (s)','FontSize',fontsize);
ylabel('Distance (m)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
% velocity
subplot(5,1,4)
plot(t,dx,'LineWidth',linewidth,'MarkerSize',marksize);
axis([min(t) max(t) min([min(x) min(dx)]) max([max(x) max(dx)])])
grid on;
xlabel('Time (s)','FontSize',fontsize);
ylabel('Velocity (m/s)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
% controller output
subplot(5,1,5)
plot(t,F,'LineWidth',linewidth,'MarkerSize',marksize);
axis([min(t) max(t) min(F) max(F)])
grid on;
xlabel('Time (s)','FontSize',fontsize);
ylabel('Controller Output (N)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
end