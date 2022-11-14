function [residual, g1, g2, g3] = RP_omega_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 52, 1);

%
% Model equations
%

T22 = y(22)/y(20);
T49 = y(10)*(1-params(3))+y(11)*y(29)-(params(7)*(y(29)-1)+params(8)*0.5*(y(29)-1)^2);
T71 = 1/(1+params(10)*params(6));
T72 = (params(10)*params(5)/(params(10)-1))^T71;
T84 = y(4)^(1+params(6));
T88 = (y(12)*y(13))^(params(10)*(1+params(6)));
T90 = params(2)*params(11)*y(38)^(1+params(10)*params(6))+T84*T88;
T96 = y(13)^params(10);
T98 = y(12)^(params(10)-1);
T100 = params(2)*params(11)*y(39)^(1+params(10)*params(6))+y(18)*y(4)*T96*T98;
T104 = y(12)^(1-params(10));
T105 = y(13)^(1-params(10));
T132 = y(12)^params(16);
T138 = y(12)^(params(16)-1);
T143 = y(12)^(1-params(16));
T157 = y(30)^(params(19)-1);
T158 = y(23)*exp(y(40))*params(19)*T157;
T160 = y(4)^params(20);
T164 = params(14)^params(21);
T168 = params(15)^params(22);
T172 = y(30)^params(19);
T175 = y(4)^(params(20)-1);
T185 = params(14)^(params(21)-1);
T195 = params(15)^(params(22)-1);
T220 = params(23)-y(9)/params(24)-y(19)/(y(10)*params(24));
T221 = y(10)*y(46)*T220;
T383 = y(10)*y(17)/y(16);
T387 = params(24)*params(23)-y(9)-y(19)/y(10);
lhs =1/(y(1)-y(1)*params(1))-params(1)*params(2)/(y(1)-y(1)*params(1));
rhs =y(18);
residual(1)= lhs-rhs;
lhs =params(4)/T22;
rhs =y(18)*(y(16)-1);
residual(2)= lhs-rhs;
lhs =y(11);
rhs =params(7)+params(8)*(y(29)-1);
residual(3)= lhs-rhs;
lhs =y(18)*y(10);
rhs =params(2)*y(18)*T49;
residual(4)= lhs-rhs;
lhs =y(18);
rhs =y(16)*params(2)*y(18)/y(12);
residual(5)= lhs-rhs;
lhs =y(12);
rhs =y(20)*y(24)/y(20);
residual(6)= lhs-rhs;
lhs =y(37);
rhs =T72*y(38)/y(39);
residual(7)= lhs-rhs;
lhs =y(38);
rhs =T90^T71;
residual(8)= lhs-rhs;
lhs =y(39);
rhs =T100^T71;
residual(9)= lhs-rhs;
lhs =T104*T105;
rhs =params(13)*(1-params(11))*y(37)^(1-params(10))+T105*params(11)*T104;
residual(10)= lhs-rhs;
lhs =y(34);
rhs =params(16)*y(35)/((params(16)-1)*y(36));
residual(11)= lhs-rhs;
lhs =y(35);
rhs =y(35)*params(2)*params(17)+y(18)*y(23)*y(5)*T132;
residual(12)= lhs-rhs;
lhs =y(36);
rhs =y(36)*params(2)*params(17)+y(18)*y(5)*T138;
residual(13)= lhs-rhs;
lhs =T143;
rhs =(1-params(17))*y(34)^(1-params(16))+params(17)*T143;
residual(14)= lhs-rhs;
lhs =y(11);
rhs =T158*T160*T164*T168;
residual(15)= lhs-rhs;
lhs =y(13);
rhs =T168*T164*y(23)*exp(y(40))*params(20)*T172*T175;
residual(16)= lhs-rhs;
lhs =y(14);
rhs =T168*T160*T172*y(23)*exp(y(40))*params(21)*T185;
residual(17)= lhs-rhs;
lhs =y(15);
rhs =T164*T160*T172*y(23)*exp(y(40))*params(22)*T195;
residual(18)= lhs-rhs;
lhs =y(8)*y(17);
rhs =y(46)*y(19)*y(6)/params(24);
residual(19)= lhs-rhs;
lhs =y(16)*(1-y(22)*params(13))/y(20);
rhs =y(6)*T221;
residual(20)= lhs-rhs;
lhs =y(32);
rhs =y(6)*y(19)*y(46)*params(25)/(y(10)*params(24));
residual(21)= lhs-rhs;
lhs =y(33);
rhs =y(6)*y(9)*y(46)*params(26)/params(24);
residual(22)= lhs-rhs;
lhs =y(8);
rhs =y(32)*exp(y(42))*(y(11)+y(10)*(1-params(3)))+params(15)*y(15);
residual(23)= lhs-rhs;
lhs =y(7);
rhs =y(33)*(y(11)+y(10)*(1-params(3)))+params(14)*y(14);
residual(24)= lhs-rhs;
lhs =y(3);
rhs =y(6)*y(19)*y(46)*(1-params(25))/params(24);
residual(25)= lhs-rhs;
lhs =y(2);
rhs =y(6)*y(9)*y(10)*y(46)*(1-params(26))/params(24);
residual(26)= lhs-rhs;
lhs =y(6)+y(3)+y(2)+y(1)*params(13)+y(19)*y(6);
rhs =y(5);
residual(27)= lhs-rhs;
lhs =y(31);
rhs =(1-params(3))*y(31)+y(6)*y(46)*params(23);
residual(28)= lhs-rhs;
lhs =log(y(47)/(y(47)));
rhs =log(y(47)/(y(47)))*params(27)+(1-params(27))*(params(28)*log(y(12)/(y(12)))+params(29)*log(y(5)/(y(5))))+y(41);
residual(29)= lhs-rhs;
lhs =y(16);
rhs =y(47);
residual(30)= lhs-rhs;
lhs =y(5);
rhs =T168*T164*T160*exp(y(40))*T172-params(30);
residual(31)= lhs-rhs;
lhs =params(15)*y(26);
rhs =(1-y(22)*params(13))/y(20);
residual(32)= lhs-rhs;
lhs =y(30);
rhs =y(32)+y(33)+y(29)*(y(31)-y(33)-y(32));
residual(33)= lhs-rhs;
lhs =y(25);
rhs =y(12);
residual(34)= lhs-rhs;
lhs =y(27);
rhs =1;
residual(35)= lhs-rhs;
lhs =y(21);
rhs =y(6)-y(7);
residual(36)= lhs-rhs;
lhs =y(28);
rhs =y(3)+y(2)+y(1)*params(13);
residual(37)= lhs-rhs;
lhs =y(40);
rhs =y(40)*params(33)-x(1);
residual(38)= lhs-rhs;
lhs =y(41);
rhs =y(41)*params(34)+x(2);
residual(39)= lhs-rhs;
lhs =y(42);
rhs =y(42)*params(35)-x(3);
residual(40)= lhs-rhs;
lhs =y(43);
rhs =log(y(5));
residual(41)= lhs-rhs;
lhs =y(44);
rhs =log(y(6));
residual(42)= lhs-rhs;
lhs =y(9);
rhs =params(44)*(1+y(19)*params(46))^(-params(43));
residual(43)= lhs-rhs;
lhs =y(46);
rhs =params(39);
residual(44)= lhs-rhs;
lhs =y(45);
rhs =(y(45))+params(41)*((y(6)-y(7)-y(5)*y(48))/y(5)-((y(6))-(y(7)))/(y(5)))-x(5);
residual(45)= lhs-rhs;
lhs =y(6);
rhs =y(5)*y(48)+y(7)+y(8)*y(45);
residual(46)= lhs-rhs;
lhs =y(45);
rhs =1+T383*T387/y(19)-y(17)*params(24)/y(46);
residual(47)= lhs-rhs;
lhs =y(48);
rhs =y(49)/(y(5));
residual(48)= lhs-rhs;
lhs =y(49);
rhs =y(49)*params(49)+x(4);
residual(49)= lhs-rhs;
lhs =y(50);
rhs =log(y(1)-y(1)*params(1))-params(20)*params(5)*(1-params(1)*params(2))/(params(5)*y(1)*(1-params(1))*params(12))/(1+params(6))+params(4)*log(T22);
residual(50)= lhs-rhs;
lhs =y(51);
rhs =y(50)+params(2)*y(50);
residual(51)= lhs-rhs;
lhs =y(52);
rhs =y(19)/(y(10)*0.264276472435772);
residual(52)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(52, 52);

  %
  % Jacobian matrix
  %

T449 = getPowerDeriv(T90,T71,1);
T454 = getPowerDeriv(T100,T71,1);
T457 = getPowerDeriv(y(4),params(20),1);
T595 = getPowerDeriv(y(12)*y(13),params(10)*(1+params(6)),1);
T604 = getPowerDeriv(y(12),1-params(10),1);
T615 = getPowerDeriv(y(12),1-params(16),1);
T634 = getPowerDeriv(y(13),1-params(10),1);
T765 = getPowerDeriv(y(30),params(19),1);
T861 = ((y(47))-y(47))/((y(47))*(y(47)))/(y(47)/(y(47)));
  g1(1,1)=(-(1-params(1)))/((y(1)-y(1)*params(1))*(y(1)-y(1)*params(1)))-(-(params(1)*params(2)*(1-params(1))))/((y(1)-y(1)*params(1))*(y(1)-y(1)*params(1)));
  g1(1,18)=(-1);
  g1(2,16)=(-y(18));
  g1(2,18)=(-(y(16)-1));
  g1(2,20)=(-(params(4)*(-y(22))/(y(20)*y(20))))/(T22*T22);
  g1(2,22)=(-(params(4)*1/y(20)))/(T22*T22);
  g1(3,11)=1;
  g1(3,29)=(-params(8));
  g1(4,10)=y(18)-params(2)*y(18)*(1-params(3));
  g1(4,11)=(-(params(2)*y(18)*y(29)));
  g1(4,18)=y(10)-params(2)*T49;
  g1(4,29)=(-(params(2)*y(18)*(y(11)-(params(7)+params(8)*0.5*2*(y(29)-1)))));
  g1(5,12)=(-((-(y(16)*params(2)*y(18)))/(y(12)*y(12))));
  g1(5,16)=(-(params(2)*y(18)/y(12)));
  g1(5,18)=1-params(2)*y(16)/y(12);
  g1(6,12)=1;
  g1(6,24)=(-1);
  g1(7,37)=1;
  g1(7,38)=(-(T72/y(39)));
  g1(7,39)=(-((-(T72*y(38)))/(y(39)*y(39))));
  g1(8,4)=(-(T88*getPowerDeriv(y(4),1+params(6),1)*T449));
  g1(8,12)=(-(T449*T84*y(13)*T595));
  g1(8,13)=(-(T449*T84*y(12)*T595));
  g1(8,38)=1-T449*params(2)*params(11)*getPowerDeriv(y(38),1+params(10)*params(6),1);
  g1(9,4)=(-(T98*y(18)*T96*T454));
  g1(9,12)=(-(T454*y(18)*y(4)*T96*getPowerDeriv(y(12),params(10)-1,1)));
  g1(9,13)=(-(T454*T98*y(18)*y(4)*getPowerDeriv(y(13),params(10),1)));
  g1(9,18)=(-(T454*T98*y(4)*T96));
  g1(9,39)=1-T454*params(2)*params(11)*getPowerDeriv(y(39),1+params(10)*params(6),1);
  g1(10,12)=T105*T604-T105*params(11)*T604;
  g1(10,13)=T104*T634-params(11)*T104*T634;
  g1(10,37)=(-(params(13)*(1-params(11))*getPowerDeriv(y(37),1-params(10),1)));
  g1(11,34)=1;
  g1(11,35)=(-(params(16)/((params(16)-1)*y(36))));
  g1(11,36)=(-((-(params(16)*y(35)*(params(16)-1)))/((params(16)-1)*y(36)*(params(16)-1)*y(36))));
  g1(12,5)=(-(y(18)*y(23)*T132));
  g1(12,12)=(-(y(18)*y(23)*y(5)*getPowerDeriv(y(12),params(16),1)));
  g1(12,18)=(-(T132*y(23)*y(5)));
  g1(12,23)=(-(T132*y(18)*y(5)));
  g1(12,35)=1-params(2)*params(17);
  g1(13,5)=(-(y(18)*T138));
  g1(13,12)=(-(y(18)*y(5)*getPowerDeriv(y(12),params(16)-1,1)));
  g1(13,18)=(-(y(5)*T138));
  g1(13,36)=1-params(2)*params(17);
  g1(14,12)=T615-params(17)*T615;
  g1(14,34)=(-((1-params(17))*getPowerDeriv(y(34),1-params(16),1)));
  g1(15,4)=(-(T168*T164*T158*T457));
  g1(15,11)=1;
  g1(15,23)=(-(T168*T164*T160*T157*exp(y(40))*params(19)));
  g1(15,30)=(-(T168*T164*T160*y(23)*exp(y(40))*params(19)*getPowerDeriv(y(30),params(19)-1,1)));
  g1(15,40)=(-(T158*T160*T164*T168));
  g1(16,4)=(-(T168*T164*y(23)*exp(y(40))*params(20)*T172*getPowerDeriv(y(4),params(20)-1,1)));
  g1(16,13)=1;
  g1(16,23)=(-(T168*T164*T175*T172*exp(y(40))*params(20)));
  g1(16,30)=(-(T168*T164*T175*y(23)*exp(y(40))*params(20)*T765));
  g1(16,40)=(-(T168*T164*y(23)*exp(y(40))*params(20)*T172*T175));
  g1(17,4)=(-(T168*T185*T172*y(23)*exp(y(40))*params(21)*T457));
  g1(17,14)=1;
  g1(17,23)=(-(T168*T185*T160*T172*exp(y(40))*params(21)));
  g1(17,30)=(-(T168*T185*T160*y(23)*exp(y(40))*params(21)*T765));
  g1(17,40)=(-(T168*T160*T172*y(23)*exp(y(40))*params(21)*T185));
  g1(18,4)=(-(T195*T164*T172*y(23)*exp(y(40))*params(22)*T457));
  g1(18,15)=1;
  g1(18,23)=(-(T195*T164*T160*T172*exp(y(40))*params(22)));
  g1(18,30)=(-(T195*T164*T160*y(23)*exp(y(40))*params(22)*T765));
  g1(18,40)=(-(T164*T160*T172*y(23)*exp(y(40))*params(22)*T195));
  g1(19,6)=(-(y(46)*y(19)/params(24)));
  g1(19,8)=y(17);
  g1(19,17)=y(8);
  g1(19,19)=(-(y(46)*y(6)/params(24)));
  g1(19,46)=(-(y(19)*y(6)/params(24)));
  g1(20,6)=(-T221);
  g1(20,9)=(-(y(6)*y(10)*y(46)*(-(1/params(24)))));
  g1(20,10)=(-(y(6)*(y(46)*T220+y(10)*y(46)*(-((-(y(19)*params(24)))/(y(10)*params(24)*y(10)*params(24)))))));
  g1(20,16)=(1-y(22)*params(13))/y(20);
  g1(20,19)=(-(y(6)*y(10)*y(46)*(-(1/(y(10)*params(24))))));
  g1(20,20)=(-(y(16)*(1-y(22)*params(13))))/(y(20)*y(20));
  g1(20,22)=y(16)*(-params(13))/y(20);
  g1(20,46)=(-(y(6)*y(10)*T220));
  g1(21,6)=(-(y(19)*y(46)*params(25)/(y(10)*params(24))));
  g1(21,10)=(-((-(params(24)*y(6)*y(19)*y(46)*params(25)))/(y(10)*params(24)*y(10)*params(24))));
  g1(21,19)=(-(y(6)*y(46)*params(25)/(y(10)*params(24))));
  g1(21,32)=1;
  g1(21,46)=(-(y(6)*y(19)*params(25)/(y(10)*params(24))));
  g1(22,6)=(-(y(9)*y(46)*params(26)/params(24)));
  g1(22,9)=(-(y(6)*y(46)*params(26)/params(24)));
  g1(22,33)=1;
  g1(22,46)=(-(y(6)*y(9)*params(26)/params(24)));
  g1(23,8)=1;
  g1(23,10)=(-(y(32)*(1-params(3))*exp(y(42))));
  g1(23,11)=(-(y(32)*exp(y(42))));
  g1(23,15)=(-params(15));
  g1(23,32)=(-(exp(y(42))*(y(11)+y(10)*(1-params(3)))));
  g1(23,42)=(-(y(32)*exp(y(42))*(y(11)+y(10)*(1-params(3)))));
  g1(24,7)=1;
  g1(24,10)=(-((1-params(3))*y(33)));
  g1(24,11)=(-y(33));
  g1(24,14)=(-params(14));
  g1(24,33)=(-(y(11)+y(10)*(1-params(3))));
  g1(25,3)=1;
  g1(25,6)=(-(y(19)*y(46)*(1-params(25))/params(24)));
  g1(25,19)=(-(y(6)*y(46)*(1-params(25))/params(24)));
  g1(25,46)=(-(y(6)*y(19)*(1-params(25))/params(24)));
  g1(26,2)=1;
  g1(26,6)=(-(y(9)*y(10)*y(46)*(1-params(26))/params(24)));
  g1(26,9)=(-(y(6)*y(10)*y(46)*(1-params(26))/params(24)));
  g1(26,10)=(-(y(6)*y(9)*y(46)*(1-params(26))/params(24)));
  g1(26,46)=(-(y(6)*y(9)*y(10)*(1-params(26))/params(24)));
  g1(27,1)=params(13);
  g1(27,2)=1;
  g1(27,3)=1;
  g1(27,5)=(-1);
  g1(27,6)=1+y(19);
  g1(27,19)=y(6);
  g1(28,6)=(-(y(46)*params(23)));
  g1(28,31)=1-(1-params(3));
  g1(28,46)=(-(y(6)*params(23)));
  g1(29,5)=(-((1-params(27))*params(29)*((y(5))-y(5))/((y(5))*(y(5)))/(y(5)/(y(5)))));
  g1(29,12)=(-((1-params(27))*params(28)*((y(12))-y(12))/((y(12))*(y(12)))/(y(12)/(y(12)))));
  g1(29,41)=(-1);
  g1(29,47)=T861-params(27)*T861;
  g1(30,16)=1;
  g1(30,47)=(-1);
  g1(31,4)=(-(T168*T164*exp(y(40))*T172*T457));
  g1(31,5)=1;
  g1(31,30)=(-(T168*T164*T160*exp(y(40))*T765));
  g1(31,40)=(-(T168*T164*T160*exp(y(40))*T172));
  g1(32,20)=(-((-(1-y(22)*params(13)))/(y(20)*y(20))));
  g1(32,22)=(-((-params(13))/y(20)));
  g1(32,26)=params(15);
  g1(33,29)=(-(y(31)-y(33)-y(32)));
  g1(33,30)=1;
  g1(33,31)=(-y(29));
  g1(33,32)=(-(1-y(29)));
  g1(33,33)=(-(1-y(29)));
  g1(34,12)=(-1);
  g1(34,25)=1;
  g1(35,27)=1;
  g1(36,6)=(-1);
  g1(36,7)=1;
  g1(36,21)=1;
  g1(37,1)=(-params(13));
  g1(37,2)=(-1);
  g1(37,3)=(-1);
  g1(37,28)=1;
  g1(38,40)=1-params(33);
  g1(39,41)=1-params(34);
  g1(40,42)=1-params(35);
  g1(41,5)=(-(1/y(5)));
  g1(41,43)=1;
  g1(42,6)=(-(1/y(6)));
  g1(42,44)=1;
  g1(43,9)=1;
  g1(43,19)=(-(params(44)*params(46)*getPowerDeriv(1+y(19)*params(46),(-params(43)),1)));
  g1(44,46)=1;
  g1(45,5)=(-(params(41)*((y(5)*(-y(48))-(y(6)-y(7)-y(5)*y(48)))/(y(5)*y(5))-(-((y(6))-(y(7))))/((y(5))*(y(5))))));
  g1(45,6)=(-(params(41)*(1/y(5)-1/(y(5)))));
  g1(45,7)=(-(params(41)*((-1)/y(5)-(-1)/(y(5)))));
  g1(45,48)=(-(params(41)*(-y(5))/y(5)));
  g1(46,5)=(-y(48));
  g1(46,6)=1;
  g1(46,7)=(-1);
  g1(46,8)=(-y(45));
  g1(46,45)=(-y(8));
  g1(46,48)=(-y(5));
  g1(47,9)=(-((-T383)/y(19)));
  g1(47,10)=(-((T387*y(17)/y(16)+T383*(-((-y(19))/(y(10)*y(10)))))/y(19)));
  g1(47,16)=(-(T387*(-(y(10)*y(17)))/(y(16)*y(16))/y(19)));
  g1(47,17)=(-(T387*y(10)/y(16)/y(19)-params(24)/y(46)));
  g1(47,19)=(-((y(19)*T383*(-(1/y(10)))-T383*T387)/(y(19)*y(19))));
  g1(47,45)=1;
  g1(47,46)=(-(y(17)*params(24)))/(y(46)*y(46));
  g1(48,5)=(-((-y(49))/((y(5))*(y(5)))));
  g1(48,48)=1;
  g1(48,49)=(-(1/(y(5))));
  g1(49,49)=1-params(49);
  g1(50,1)=(-((1-params(1))/(y(1)-y(1)*params(1))-(-(params(20)*params(5)*(1-params(1)*params(2))*params(12)*params(5)*(1-params(1))))/(params(5)*y(1)*(1-params(1))*params(12)*params(5)*y(1)*(1-params(1))*params(12))/(1+params(6))));
  g1(50,20)=(-(params(4)*(-y(22))/(y(20)*y(20))/T22));
  g1(50,22)=(-(params(4)*1/y(20)/T22));
  g1(50,50)=1;
  g1(51,50)=(-(1+params(2)));
  g1(51,51)=1;
  g1(52,10)=(-((-(y(19)*0.264276472435772))/(y(10)*0.264276472435772*y(10)*0.264276472435772)));
  g1(52,19)=(-(1/(y(10)*0.264276472435772)));
  g1(52,52)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],52,2704);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],52,140608);
end
end
end
end
