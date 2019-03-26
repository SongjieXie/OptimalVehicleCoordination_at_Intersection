addpath( '/Users/xiesongjie/Desktop/casadi-osx-matlabR2015a-v3.4.5')
import casadi.*
%%
%The initial condition
rat = 1000/(60*60);%The unit switching from km/h to m/s 
v_0 = rat*[70,75,80,85];
p_0 = [-160,-163,-166,-166];

Q= 1; Vref = rat*80;%penalty of quaritic form of v
R= 1; %penalty of quaritic form of a

N = 4;% the number of vehicle
Pin = 0;
Pout = 10;
K = 200;
L = 30;

%%
%The initiate the optimize parameter(x_ii) , constaints(g_ii), cost(f)

x_ii = {};%The input parameters which contain a, v, p and t
lbx_ii = [];
ubx_ii = [];

g_ii = {};%The constraints contain {p(tin) = pin and p(tout) = pout},{tout>tin},{P=...},{v=...}
lbg_ii = [];
ubg_ii = [];

f = 0; %The cost function

%%
for ii = 1:N %The iith vehicle
    Tin = MX.sym(['Tin' int2str(ii)]);
    lbx_ii = [lbx_ii;-Inf];
    ubx_ii = [ubx_ii;Inf];
    x_ii = [x_ii;{Tin}];
    
    Tout = MX.sym(['Tout' int2str(ii)]);
    lbx_ii = [lbx_ii;-Inf];
    ubx_ii = [ubx_ii;Inf];
    x_ii = [x_ii;{Tout}];
    
    
    
%     t0 = MX.sym(['t' int2str(ii) '0']); %Add the initial time t0 = 0
%     x_ii = [x_ii;{t0}];
%     lbx_ii = [lbx_ii ;0];
%     ubx_ii = [ubx_ii ;0];
    
    
    a0 = MX.sym(['a' int2str(ii) '0']);%Add the initial accelerator a0 = 0
    x_ii = [x_ii;{a0}];
    lbx_ii = [lbx_ii;0];
    ubx_ii = [ubx_ii;0];
    
    
    
    v0 = MX.sym(['v' int2str(ii) '0']);%Add the initial velocity as the first(3) parameter
    x_ii = [x_ii;{v0}];
    lbx_ii = [lbx_ii; v_0(ii)];
    ubx_ii = [ubx_ii; v_0(ii)];
    
    
    p0 = MX.sym(['p' int2str(ii) '0']);%Add the initial position as the second (4) parameter
    x_ii = [x_ii;{p0}];
    lbx_ii = [lbx_ii; p_0(ii)];
    ubx_ii = [ubx_ii; p_0(ii)];

    g_ii = [g_ii; {Tout - Tin} ];%%%%%%%  unsettle  !!!
    lbg_ii = [lbg_ii; 0];
    ubg_ii = [ubg_ii; Inf];
    
    T1 = Tin/K;
    

    for i = 1:K% in the time interval 0-15
        
%         t = MX.sym(['t' int2str(ii) '+' int2str(i)]);%define t
%         lbx_ii = [lbx_ii;i*T1];
%         ubx_ii = [ubx_ii;i*T1];
%         x_ii = [x_ii;{t}];
        
        a = MX.sym(['a' int2str(ii) '+' int2str(i)]);%define a
        lbx_ii = [lbx_ii;-2];
        ubx_ii = [ubx_ii;2];
        x_ii = [x_ii;{a}];
        
        v = MX.sym(['v' int2str(ii) '+' int2str(i)]);%define v
        lbx_ii = [lbx_ii;0];
        ubx_ii = [ubx_ii;Inf];
        x_ii = [x_ii;{v}];
        
        
        VV = x_ii(end);
        V = x_ii(end-3);
        A = x_ii(end-4);
        g_ii = [g_ii;{VV{:}-V{:}-A{:}*T1}];%Vk+1 - (Vk + ak*T)
        lbg_ii = [lbg_ii;0];
        ubg_ii = [ubg_ii;0];
        
        p = MX.sym(['p' int2str(ii) '+' int2str(i)]);%define p
        lbx_ii = [lbx_ii;-Inf];
        ubx_ii = [ubx_ii;Inf];
        x_ii = [x_ii;{p}];
        
        
        PP = x_ii(end);
        P = x_ii(end-3);
        V = x_ii(end-4);
        A = x_ii(end-5);
         
        g_ii = [g_ii;{PP{:}-...
            (P{:}+V{:}*T1+0.5*A{:}*T1^2)}];% pk+1 = pk + vkT + 0.5*ak*T^2
        lbg_ii = [lbg_ii;0];
        ubg_ii = [ubg_ii;0];
        
        V = x_ii(end-1);
        A = x_ii(end-2);
        f = f+ Q*(V{:}-Vref)^2 + R*A{:}^2;

    end
%------------------------------------------------------------------------------------------------------------
    f = f + 100*(Tout-Tin)^2;  % Add a term (1)
%-------------------------------------------------------------------------------------------------------------
    pp = x_ii(end);   %P(Tin) = Pin
    g_ii = [g_ii;{pp{:}-Pin}];
    lbg_ii = [lbg_ii;0];
    ubg_ii = [ubg_ii;0];
    
    T1 = (Tout-Tin)/L;
    for i = 1:L
%         t = MX.sym(['t' int2str(ii) '+' int2str(i)]);%define t
%         lbx_ii = [lbx_ii;i*T1];
%         ubx_ii = [ubx_ii;i*T1];
%         x_ii = [x_ii;{t}];
        
        a = MX.sym(['a' int2str(ii) '+' int2str(i)]);%define a
        lbx_ii = [lbx_ii;-Inf];
        ubx_ii = [ubx_ii;Inf];
        x_ii = [x_ii;{a}];
        
        v = MX.sym(['v' int2str(ii) '+' int2str(i)]);%define v
        lbx_ii = [lbx_ii;0];
        ubx_ii = [ubx_ii;Inf];
        x_ii = [x_ii;{v}];
        
        
        VV = x_ii(end);
        V = x_ii(end-3);
        A = x_ii(end-4);
        g_ii = [g_ii;{VV{:}-V{:}-A{:}*T1}];%Vk+1 - (Vk + ak*T)
        lbg_ii = [lbg_ii;0];
        ubg_ii = [ubg_ii;0];
        
        
        p = MX.sym(['p' int2str(ii) '+' int2str(i)]);%define p
        lbx_ii = [lbx_ii;-Inf];
        ubx_ii = [ubx_ii;Inf];
        x_ii = [x_ii;{p}];
        
        
        PP = x_ii(end);
        P = x_ii(end-3);
        V = x_ii(end-4);
        A = x_ii(end-5);
         
        g_ii = [g_ii;{PP{:}-...
            (P{:}+V{:}*T1+0.5*A{:}*T1^2)}];% pk+1 = pk + vkT + 0.5*ak*T^2
        lbg_ii = [lbg_ii;0];
        ubg_ii = [ubg_ii;0];
        
        V = x_ii(end-1);
        A = x_ii(end-2);
        f = f+ Q*(V{:}-Vref)^2 + R*A{:}^2;
    end
    
    pp = x_ii(end);   %P(Tout) = Pout
    g_ii = [g_ii;{pp{:}-Pout}];
    lbg_ii = [lbg_ii;0];
    ubg_ii = [ubg_ii;0];
    
    if ii > 1
        Tinkk_index = 3*(K+L) + 5 -1;
        Tinkk = x_ii(end-Tinkk_index);
        
        Toutk_index = 2*(3*(K+L)+5) -2;
        Toutk = x_ii(end-Toutk_index);
        
        g_ii = [g_ii; {Tinkk{:}-Toutk{:}}];
        lbg_ii = [lbg_ii; 0];
        ubg_ii = [ubg_ii; Inf];
    end
end


%%

%nlp structure
intersection_problem = struct(...
    'f', f,...
    'x', vertcat(x_ii{:}),...
    'g', vertcat(g_ii{:})...
);

S = nlpsol('S', 'ipopt', intersection_problem);

disp(S);

X_initial = ones((3*(K+L) + 5)*N,1);

r = S('x0',X_initial,...
      'lbg',lbg_ii,'ubg',ubg_ii,'lbx',lbx_ii,'ubx',ubx_ii);
x_opt = r.x;

%%
%Result processing

length = 3*(K+L) + 5;

Tinl = [];
Toutl = [];
for i = 0:N-1
    Tinl = [Tinl;x_opt(1+i*length)];
    Toutl = [Toutl;x_opt(2+i*length)];
end



P = [];
for ii = 0:N-1
    ll = ((K+L)*3+5)*ii;
    for i = 0:(K+L)
        P= [P;x_opt(ll+5+3*i)];
    end
end

Vl = [];

for ii = 0:N-1
    ll = ((K+L)*3+5)*ii;
    for i = 0:(K+L)
        Vl= [Vl;x_opt(ll+4+3*i)];
    end
end

Al = [];

for ii = 0:N-1
    ll = ((K+L)*3+5)*ii;
    for i = 0:(K+L)
        Al= [Al;x_opt(ll+3+3*i)];
    end
end

Positionl = full(P);

PositionM = [];
for i =  1:N
    PositionM = [PositionM,Positionl((i-1)*(K+L+1)+1:i*(K+L+1))];
end

Velocity = full(Vl);
VelM = [];
for i =  1:N
    VelM = [VelM,Velocity((i-1)*(K+L+1)+1:i*(K+L+1))];
end


Acceleration = full(Al);
AccM = [];
for i = 1: N
    AccM = [AccM,Acceleration((i-1)*(K+L+1)+1:i*(K+L+1))];
end

Toutlist = full(Toutl);

Tinlist = full(Tinl);

time = [];
for i = 1:N
    ttK = 0:Tinlist(i)/K:Tinlist(i);
    ttL = Tinlist(i)+(Toutlist(i)-Tinlist(i))/L:(Toutlist(i)-Tinlist(i))/L:Toutlist(i);
    ttt = [ttK,ttL];
    time = [time;ttt];
    trick = time';%  The time line!!! 
end


%% Plot

ttt = 0:0.1:9;

title('The original objective fucntion')
subplot(3,1,1);

plot(trick(:,1),PositionM(:,1),'r-',ttt,p_0(1)+v_0(1)*ttt ,'r:');
hold all;
plot(trick(:,2),PositionM(:,2),'b-',ttt,p_0(2)+v_0(2)*ttt ,'b:');
plot(trick(:,3),PositionM(:,3),'g-',ttt,p_0(3)+v_0(3)*ttt ,'g:');
plot(trick(:,4),PositionM(:,4),'k-',ttt,p_0(4)+v_0(4)*ttt ,'k:');

rectangle('position',[Tinlist(1) 0 Toutlist(1)-Tinlist(1) 10] );
rectangle('position',[Tinlist(2) 0 Toutlist(2)-Tinlist(2) 10] );
rectangle('position',[Tinlist(3) 0 Toutlist(3)-Tinlist(3) 10] );
rectangle('position',[Tinlist(4) 0 Toutlist(4)-Tinlist(4) 10] );


axis([6 10 -20 10]);
xlabel('Time(s)');

ylabel('Position(m)');

legend('Vehicle1_{pre}','Vehicle1_{post}','Vehicle2_{pre}','Vehicle2_{post}','Vehicle3_{pre}','Vehicle3_{post}','Vehilce4_{post}','Vehicle4_{pre}');

grid on;


subplot(3,1,2);

plot(trick(:,1),(1/rat)*VelM(:,1),'r-');
hold all;
plot(trick(:,2),(1/rat)*VelM(:,2),'b-');
plot(trick(:,3),(1/rat)*VelM(:,3),'g-');
plot(trick(:,4),(1/rat)*VelM(:,4),'k-');

axis([0 10 70 90]);
xlabel('Time(s)');

ylabel('Velocity(Km/h)');

legend('Vehicle1','Vehicle2','Vehicle3','Vehicle4');

grid on;


subplot(3,1,3);
plot(trick(:,1),AccM(:,1),'r-');
hold all;
plot(trick(:,2),AccM(:,2),'b-');
plot(trick(:,3),AccM(:,3),'g-');
plot(trick(:,4),AccM(:,4),'k-');

axis([0 9 -4 5]);
xlabel('Time(s)');

ylabel('Acceleration(m/s^2)');

legend('Vehicle1','Vehicle2','Vehicle3','Vehicle4');

grid on;

suptitle('The objective function with additional (T_{out}-T_{in})^2')


    