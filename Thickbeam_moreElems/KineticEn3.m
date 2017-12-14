function [Ekin_avL,Epot_av,Epot_avx,Epot_avy,P,Px,Py,Ekin_avLx,Ekin_avLy,ve] = KineticEn3(PhiQa_p1,lpa2,Ma,Ka,w,La,AMP)
[dofa,c] = size(PhiQa_p1);

DD = (Ka-w^2*Ma);

Ek=zeros(c,1);
Ex=zeros(c,1);
Ey=zeros(c,1);
E_pery=zeros(c,1);
E_perx=zeros(c,1);

Ekin_av=zeros(c,1);
Epot_av=zeros(c,1);
Etot_av=zeros(c,1);
P=zeros(c,1);
ve=zeros(c,1);

VV=[PhiQa_p1;PhiQa_p1*diag(lpa2)];
FF = DD*VV;

Vx=zeros(max(size(VV)),c);
Vy=zeros(max(size(VV)),c);
Fx=zeros(max(size(VV)),c);
Fy=zeros(max(size(VV)),c);
for qq=1:2:max(size(VV))
    Vx(qq,:)=(VV(qq,:));
    Fx(qq,:)=FF(qq,:);
    Vy(qq+1,:)=(VV(qq+1,:));
    Fy(qq+1,:)=FF(qq+1,:);
end
Vx;
Vy;
% Total Kinetic Energy
Ek=1/2*(((1i*w).*VV*AMP)'*Ma*((1i*w).*VV*AMP));
% Kinetic Energy per direction
Ex=1/2*(((1i*w).*Vx*AMP)'*Ma*((1i*w).*Vx*AMP));
Ey=1/2*(((1i*w).*Vy*AMP)'*Ma*((1i*w).*Vy*AMP));
E_pery=100*real(Ey)/real(Ek);
E_perx=100*real(Ex)/real(Ek);
% Time-averaged kinetic energy
Ekin_av=1/2*real(Ek);
% Time-averaged kinetic energy per unit length
Ekin_avL=1/2*real(Ek)/La;
Ekin_avLx=1/2*real(Ex)/La;
Ekin_avLy=1/2*real(Ey)/La;
Ekin_pery=100*real(Ekin_avLy)/real(Ekin_avL);
Ekin_perx=100*real(Ekin_avLx)/real(Ekin_avL);
% Time-Averaged Potencial Energy
Epot_av=1/4*real((AMP'*VV'*(Ka)*VV*AMP))/La;
Epot_avx=1/4*real((AMP'*Vx'*(Ka)*Vx*AMP))/La;
Epot_avy=1/4*real((AMP'*Vy'*(Ka)*Vy*AMP))/La;
Epot_perx = Epot_avx/Epot_av;
Epot_pery = Epot_avy/Epot_av;
% Time-averaged total energy per unit length
Etot_av=(Ekin_av+Epot_av)/(La);
% Time-averaged Energy flow
PP=real(1/2*[(AMP'*FF(1:dofa,:)'*1i*w*VV(1:dofa,:)*AMP),(AMP'*FF(dofa+1:2*dofa,:)'*1i*w*VV(dofa+1:2*dofa,:)*AMP)]);
P=PP(1);
PPx=real(1/2*[(AMP'*Fx(1:dofa,:)'*1i*w*Vx(1:dofa,:)*AMP),(AMP'*Fx(dofa+1:2*dofa,:)'*1i*w*Vx(dofa+1:2*dofa,:)*AMP)]);
Px=PPx(1);
PPy=real(1/2*[(AMP'*Fy(1:dofa,:)'*1i*w*Vy(1:dofa,:)*AMP),(AMP'*Fy(dofa+1:2*dofa,:)'*1i*w*Vy(dofa+1:2*dofa,:)*AMP)]);
Py=PPy(1);
P_perx = 100*Px./P;
P_pery = 100*Py/P;
% Energy Velocity
ve=P*La/(Ekin_av+Epot_av);%(P(1)*cos(theta)+P(1)*sin(theta))*L/(Ekin_av+Epot_av); %sqrt(vve(1)^2+vve(2)^2)
