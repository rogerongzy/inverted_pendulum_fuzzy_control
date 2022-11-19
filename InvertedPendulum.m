function [t,theta,dtheta,ddtheta,x,dx,ddx,F]=InvertedPendulum(t0,theta0,...
    dtheta0,ddtheta0,x0,dx0,ddx0,F0,inputF,t_step,varargin)

% Get force
if ~isempty(varargin)
    if t0>varargin{2} && t0<varargin{3}
        inputF=inputF+varargin{1};
    end
end
fx_F=@(t,F,inputF) -100.*F+100.*inputF;
fx=@(t,F) fx_F(t,F,inputF);
[t,F]=ODE_RK(t0,F0,t_step,fx,1);
t=t(end);
F=F(end);
% Get theta
theta=theta0+t_step.*dtheta0+0.5.*ddtheta0.*t_step.^2;
% Get dtheta
fx_dtheta=@(t,theta,dtheta,F) (9.81.*sin(theta)+cos(theta).*...
    ((-F-0.25.*dtheta.^2.*sin(theta))./(1.5)))./...
    (0.5.*(4/3-1/3.*cos(theta).^2));
fx=@(t,dtheta) fx_dtheta(t,theta,dtheta,F);
[~,dtheta]=ODE_RK(t0,dtheta0,t_step,fx,1);
dtheta=dtheta(end);
% Get ddtheta
ddtheta=(9.81.*sin(theta)+cos(theta).*((-F-0.25.*dtheta.^2.*sin(theta))./...
    (1.5)))./(0.5.*(4/3-1/3.*cos(theta).^2));
% Get car position
x=x0+t_step.*dx0+0.5.*ddx0.*t_step.^2;
% Get car speed
dx=dx0+t_step.*ddx0;
% Get car acceleration
% l=0.5; m=0.5; M=1;
l=0.3; m=0.3; M=1.5;
ddx=(F-m*l*ddtheta*cos(theta)+m*l*dtheta^2*sin(theta))./(m+M);

end