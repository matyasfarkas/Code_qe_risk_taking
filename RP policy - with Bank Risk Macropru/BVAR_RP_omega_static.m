function [residual, g1, g2, g3] = BVAR_RP_omega_static(y, x, params)
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

residual = zeros( 133, 1);

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
T388 = y(10)*y(17)/y(16);
T392 = params(24)*params(23)-y(9)-y(19)/y(10);
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
rhs =y(40)*params(33)+x(1);
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
rhs =1-(1-params(39))-((y(6)-y(7))/y(5)-((y(6))-(y(7)))/(y(5)))*params(42);
residual(44)= lhs-rhs;
lhs =y(45);
rhs =(y(45))+params(41)*(y(133)-(y(133)));
residual(45)= lhs-rhs;
lhs =y(6);
rhs =y(7)+y(8)*y(45)+y(5)*y(48);
residual(46)= lhs-rhs;
lhs =y(45);
rhs =1+T388*T392/y(19)-y(17)*params(24)/y(46);
residual(47)= lhs-rhs;
lhs =y(48);
rhs =y(49)/(y(5));
residual(48)= lhs-rhs;
lhs =y(133);
rhs =y(19)/(y(10)*0.284025416687741);
residual(49)= lhs-rhs;
lhs =y(49);
rhs =y(49)*params(49)+x(4);
residual(50)= lhs-rhs;
lhs =y(50);
rhs =y(19)/y(10);
residual(51)= lhs-rhs;
lhs =y(51);
rhs =y(12);
residual(52)= lhs-rhs;
lhs =y(52);
rhs =y(51);
residual(53)= lhs-rhs;
lhs =y(53);
rhs =y(52);
residual(54)= lhs-rhs;
lhs =y(54);
rhs =y(53);
residual(55)= lhs-rhs;
lhs =y(55);
rhs =y(54);
residual(56)= lhs-rhs;
lhs =y(56);
rhs =y(55);
residual(57)= lhs-rhs;
lhs =y(57);
rhs =y(56);
residual(58)= lhs-rhs;
lhs =y(58);
rhs =y(57);
residual(59)= lhs-rhs;
lhs =y(59);
rhs =y(58);
residual(60)= lhs-rhs;
lhs =y(60);
rhs =y(59);
residual(61)= lhs-rhs;
lhs =y(61);
rhs =y(60);
residual(62)= lhs-rhs;
lhs =y(62);
rhs =y(61);
residual(63)= lhs-rhs;
lhs =y(63);
rhs =y(62);
residual(64)= lhs-rhs;
lhs =y(64);
rhs =y(63);
residual(65)= lhs-rhs;
lhs =y(65);
rhs =y(64);
residual(66)= lhs-rhs;
lhs =y(66);
rhs =y(65);
residual(67)= lhs-rhs;
lhs =y(67);
rhs =y(66);
residual(68)= lhs-rhs;
lhs =y(68);
rhs =y(67);
residual(69)= lhs-rhs;
lhs =y(69);
rhs =y(68);
residual(70)= lhs-rhs;
lhs =y(70);
rhs =y(69);
residual(71)= lhs-rhs;
lhs =y(71);
rhs =y(70);
residual(72)= lhs-rhs;
lhs =y(72);
rhs =y(71);
residual(73)= lhs-rhs;
lhs =y(73);
rhs =y(72);
residual(74)= lhs-rhs;
lhs =y(74);
rhs =y(73);
residual(75)= lhs-rhs;
lhs =y(75);
rhs =y(74);
residual(76)= lhs-rhs;
lhs =y(76);
rhs =y(75);
residual(77)= lhs-rhs;
lhs =y(77);
rhs =y(76);
residual(78)= lhs-rhs;
lhs =y(78);
rhs =y(77);
residual(79)= lhs-rhs;
lhs =y(79);
rhs =y(78);
residual(80)= lhs-rhs;
lhs =y(80);
rhs =y(79);
residual(81)= lhs-rhs;
lhs =y(81);
rhs =y(80);
residual(82)= lhs-rhs;
lhs =y(82);
rhs =y(81);
residual(83)= lhs-rhs;
lhs =y(83);
rhs =y(82);
residual(84)= lhs-rhs;
lhs =y(84);
rhs =y(83);
residual(85)= lhs-rhs;
lhs =y(85);
rhs =y(84);
residual(86)= lhs-rhs;
lhs =y(86);
rhs =y(85);
residual(87)= lhs-rhs;
lhs =y(87);
rhs =y(86);
residual(88)= lhs-rhs;
lhs =y(88);
rhs =y(87);
residual(89)= lhs-rhs;
lhs =y(89);
rhs =y(88);
residual(90)= lhs-rhs;
lhs =y(90);
rhs =y(89);
residual(91)= lhs-rhs;
lhs =y(91);
rhs =y(16);
residual(92)= lhs-rhs;
lhs =y(92);
rhs =y(91);
residual(93)= lhs-rhs;
lhs =y(93);
rhs =y(92);
residual(94)= lhs-rhs;
lhs =y(94);
rhs =y(93);
residual(95)= lhs-rhs;
lhs =y(95);
rhs =y(94);
residual(96)= lhs-rhs;
lhs =y(96);
rhs =y(95);
residual(97)= lhs-rhs;
lhs =y(97);
rhs =y(96);
residual(98)= lhs-rhs;
lhs =y(98);
rhs =y(97);
residual(99)= lhs-rhs;
lhs =y(99);
rhs =y(98);
residual(100)= lhs-rhs;
lhs =y(100);
rhs =y(99);
residual(101)= lhs-rhs;
lhs =y(101);
rhs =y(100);
residual(102)= lhs-rhs;
lhs =y(102);
rhs =y(101);
residual(103)= lhs-rhs;
lhs =y(103);
rhs =y(102);
residual(104)= lhs-rhs;
lhs =y(104);
rhs =y(103);
residual(105)= lhs-rhs;
lhs =y(105);
rhs =y(104);
residual(106)= lhs-rhs;
lhs =y(106);
rhs =y(105);
residual(107)= lhs-rhs;
lhs =y(107);
rhs =y(106);
residual(108)= lhs-rhs;
lhs =y(108);
rhs =y(107);
residual(109)= lhs-rhs;
lhs =y(109);
rhs =y(108);
residual(110)= lhs-rhs;
lhs =y(110);
rhs =y(109);
residual(111)= lhs-rhs;
lhs =y(111);
rhs =y(110);
residual(112)= lhs-rhs;
lhs =y(112);
rhs =y(111);
residual(113)= lhs-rhs;
lhs =y(113);
rhs =y(112);
residual(114)= lhs-rhs;
lhs =y(114);
rhs =y(113);
residual(115)= lhs-rhs;
lhs =y(115);
rhs =y(114);
residual(116)= lhs-rhs;
lhs =y(116);
rhs =y(115);
residual(117)= lhs-rhs;
lhs =y(117);
rhs =y(116);
residual(118)= lhs-rhs;
lhs =y(118);
rhs =y(117);
residual(119)= lhs-rhs;
lhs =y(119);
rhs =y(118);
residual(120)= lhs-rhs;
lhs =y(120);
rhs =y(119);
residual(121)= lhs-rhs;
lhs =y(121);
rhs =y(120);
residual(122)= lhs-rhs;
lhs =y(122);
rhs =y(121);
residual(123)= lhs-rhs;
lhs =y(123);
rhs =y(122);
residual(124)= lhs-rhs;
lhs =y(124);
rhs =y(123);
residual(125)= lhs-rhs;
lhs =y(125);
rhs =y(124);
residual(126)= lhs-rhs;
lhs =y(126);
rhs =y(125);
residual(127)= lhs-rhs;
lhs =y(127);
rhs =y(126);
residual(128)= lhs-rhs;
lhs =y(128);
rhs =y(127);
residual(129)= lhs-rhs;
lhs =y(129);
rhs =y(128);
residual(130)= lhs-rhs;
lhs =y(130);
rhs =y(129);
residual(131)= lhs-rhs;
lhs =y(131);
rhs =log(y(1)-y(1)*params(1))-params(20)*params(5)*(1-params(1)*params(2))/(params(5)*y(1)*(1-params(1))*params(12))/(1+params(6))+params(4)*log(T22);
residual(132)= lhs-rhs;
lhs =y(132);
rhs =y(131)+params(2)*y(131);
residual(133)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(133, 133);

  %
  % Jacobian matrix
  %

T3 = (-1);
T615 = getPowerDeriv(T90,T71,1);
T620 = getPowerDeriv(T100,T71,1);
T623 = getPowerDeriv(y(4),params(20),1);
T760 = getPowerDeriv(y(12)*y(13),params(10)*(1+params(6)),1);
T769 = getPowerDeriv(y(12),1-params(10),1);
T780 = getPowerDeriv(y(12),1-params(16),1);
T799 = getPowerDeriv(y(13),1-params(10),1);
T930 = getPowerDeriv(y(30),params(19),1);
T1026 = ((y(47))-y(47))/((y(47))*(y(47)))/(y(47)/(y(47)));
  g1(1,1)=(-(1-params(1)))/((y(1)-y(1)*params(1))*(y(1)-y(1)*params(1)))-(-(params(1)*params(2)*(1-params(1))))/((y(1)-y(1)*params(1))*(y(1)-y(1)*params(1)));
  g1(1,18)=T3;
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
  g1(6,24)=T3;
  g1(7,37)=1;
  g1(7,38)=(-(T72/y(39)));
  g1(7,39)=(-((-(T72*y(38)))/(y(39)*y(39))));
  g1(8,4)=(-(T88*getPowerDeriv(y(4),1+params(6),1)*T615));
  g1(8,12)=(-(T615*T84*y(13)*T760));
  g1(8,13)=(-(T615*T84*y(12)*T760));
  g1(8,38)=1-T615*params(2)*params(11)*getPowerDeriv(y(38),1+params(10)*params(6),1);
  g1(9,4)=(-(T98*y(18)*T96*T620));
  g1(9,12)=(-(T620*y(18)*y(4)*T96*getPowerDeriv(y(12),params(10)-1,1)));
  g1(9,13)=(-(T620*T98*y(18)*y(4)*getPowerDeriv(y(13),params(10),1)));
  g1(9,18)=(-(T620*T98*y(4)*T96));
  g1(9,39)=1-T620*params(2)*params(11)*getPowerDeriv(y(39),1+params(10)*params(6),1);
  g1(10,12)=T105*T769-T105*params(11)*T769;
  g1(10,13)=T104*T799-params(11)*T104*T799;
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
  g1(14,12)=T780-params(17)*T780;
  g1(14,34)=(-((1-params(17))*getPowerDeriv(y(34),1-params(16),1)));
  g1(15,4)=(-(T168*T164*T158*T623));
  g1(15,11)=1;
  g1(15,23)=(-(T168*T164*T160*T157*exp(y(40))*params(19)));
  g1(15,30)=(-(T168*T164*T160*y(23)*exp(y(40))*params(19)*getPowerDeriv(y(30),params(19)-1,1)));
  g1(15,40)=(-(T158*T160*T164*T168));
  g1(16,4)=(-(T168*T164*y(23)*exp(y(40))*params(20)*T172*getPowerDeriv(y(4),params(20)-1,1)));
  g1(16,13)=1;
  g1(16,23)=(-(T168*T164*T175*T172*exp(y(40))*params(20)));
  g1(16,30)=(-(T168*T164*T175*y(23)*exp(y(40))*params(20)*T930));
  g1(16,40)=(-(T168*T164*y(23)*exp(y(40))*params(20)*T172*T175));
  g1(17,4)=(-(T168*T185*T172*y(23)*exp(y(40))*params(21)*T623));
  g1(17,14)=1;
  g1(17,23)=(-(T168*T185*T160*T172*exp(y(40))*params(21)));
  g1(17,30)=(-(T168*T185*T160*y(23)*exp(y(40))*params(21)*T930));
  g1(17,40)=(-(T168*T160*T172*y(23)*exp(y(40))*params(21)*T185));
  g1(18,4)=(-(T195*T164*T172*y(23)*exp(y(40))*params(22)*T623));
  g1(18,15)=1;
  g1(18,23)=(-(T195*T164*T160*T172*exp(y(40))*params(22)));
  g1(18,30)=(-(T195*T164*T160*y(23)*exp(y(40))*params(22)*T930));
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
  g1(27,5)=T3;
  g1(27,6)=1+y(19);
  g1(27,19)=y(6);
  g1(28,6)=(-(y(46)*params(23)));
  g1(28,31)=1-(1-params(3));
  g1(28,46)=(-(y(6)*params(23)));
  g1(29,5)=(-((1-params(27))*params(29)*((y(5))-y(5))/((y(5))*(y(5)))/(y(5)/(y(5)))));
  g1(29,12)=(-((1-params(27))*params(28)*((y(12))-y(12))/((y(12))*(y(12)))/(y(12)/(y(12)))));
  g1(29,41)=T3;
  g1(29,47)=T1026-params(27)*T1026;
  g1(30,16)=1;
  g1(30,47)=T3;
  g1(31,4)=(-(T168*T164*exp(y(40))*T172*T623));
  g1(31,5)=1;
  g1(31,30)=(-(T168*T164*T160*exp(y(40))*T930));
  g1(31,40)=(-(T168*T164*T160*exp(y(40))*T172));
  g1(32,20)=(-((-(1-y(22)*params(13)))/(y(20)*y(20))));
  g1(32,22)=(-((-params(13))/y(20)));
  g1(32,26)=params(15);
  g1(33,29)=(-(y(31)-y(33)-y(32)));
  g1(33,30)=1;
  g1(33,31)=(-y(29));
  g1(33,32)=(-(1-y(29)));
  g1(33,33)=(-(1-y(29)));
  g1(34,12)=T3;
  g1(34,25)=1;
  g1(35,27)=1;
  g1(36,6)=T3;
  g1(36,7)=1;
  g1(36,21)=1;
  g1(37,1)=(-params(13));
  g1(37,2)=T3;
  g1(37,3)=T3;
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
  g1(44,5)=params(42)*((-(y(6)-y(7)))/(y(5)*y(5))-(-((y(6))-(y(7))))/((y(5))*(y(5))));
  g1(44,6)=params(42)*(1/y(5)-1/(y(5)));
  g1(44,7)=params(42)*(T3/y(5)-T3/(y(5)));
  g1(44,46)=1;
  g1(46,5)=(-y(48));
  g1(46,6)=1;
  g1(46,7)=T3;
  g1(46,8)=(-y(45));
  g1(46,45)=(-y(8));
  g1(46,48)=(-y(5));
  g1(47,9)=(-((-T388)/y(19)));
  g1(47,10)=(-((T392*y(17)/y(16)+T388*(-((-y(19))/(y(10)*y(10)))))/y(19)));
  g1(47,16)=(-(T392*(-(y(10)*y(17)))/(y(16)*y(16))/y(19)));
  g1(47,17)=(-(T392*y(10)/y(16)/y(19)-params(24)/y(46)));
  g1(47,19)=(-((y(19)*T388*(-(1/y(10)))-T388*T392)/(y(19)*y(19))));
  g1(47,45)=1;
  g1(47,46)=(-(y(17)*params(24)))/(y(46)*y(46));
  g1(48,5)=(-((-y(49))/((y(5))*(y(5)))));
  g1(48,48)=1;
  g1(48,49)=(-(1/(y(5))));
  g1(49,10)=(-((-(y(19)*0.284025416687741))/(y(10)*0.284025416687741*y(10)*0.284025416687741)));
  g1(49,19)=(-(1/(y(10)*0.284025416687741)));
  g1(49,133)=1;
  g1(50,49)=1-params(49);
  g1(51,10)=(-((-y(19))/(y(10)*y(10))));
  g1(51,19)=(-(1/y(10)));
  g1(51,50)=1;
  g1(52,12)=T3;
  g1(52,51)=1;
  g1(53,51)=T3;
  g1(53,52)=1;
  g1(54,52)=T3;
  g1(54,53)=1;
  g1(55,53)=T3;
  g1(55,54)=1;
  g1(56,54)=T3;
  g1(56,55)=1;
  g1(57,55)=T3;
  g1(57,56)=1;
  g1(58,56)=T3;
  g1(58,57)=1;
  g1(59,57)=T3;
  g1(59,58)=1;
  g1(60,58)=T3;
  g1(60,59)=1;
  g1(61,59)=T3;
  g1(61,60)=1;
  g1(62,60)=T3;
  g1(62,61)=1;
  g1(63,61)=T3;
  g1(63,62)=1;
  g1(64,62)=T3;
  g1(64,63)=1;
  g1(65,63)=T3;
  g1(65,64)=1;
  g1(66,64)=T3;
  g1(66,65)=1;
  g1(67,65)=T3;
  g1(67,66)=1;
  g1(68,66)=T3;
  g1(68,67)=1;
  g1(69,67)=T3;
  g1(69,68)=1;
  g1(70,68)=T3;
  g1(70,69)=1;
  g1(71,69)=T3;
  g1(71,70)=1;
  g1(72,70)=T3;
  g1(72,71)=1;
  g1(73,71)=T3;
  g1(73,72)=1;
  g1(74,72)=T3;
  g1(74,73)=1;
  g1(75,73)=T3;
  g1(75,74)=1;
  g1(76,74)=T3;
  g1(76,75)=1;
  g1(77,75)=T3;
  g1(77,76)=1;
  g1(78,76)=T3;
  g1(78,77)=1;
  g1(79,77)=T3;
  g1(79,78)=1;
  g1(80,78)=T3;
  g1(80,79)=1;
  g1(81,79)=T3;
  g1(81,80)=1;
  g1(82,80)=T3;
  g1(82,81)=1;
  g1(83,81)=T3;
  g1(83,82)=1;
  g1(84,82)=T3;
  g1(84,83)=1;
  g1(85,83)=T3;
  g1(85,84)=1;
  g1(86,84)=T3;
  g1(86,85)=1;
  g1(87,85)=T3;
  g1(87,86)=1;
  g1(88,86)=T3;
  g1(88,87)=1;
  g1(89,87)=T3;
  g1(89,88)=1;
  g1(90,88)=T3;
  g1(90,89)=1;
  g1(91,89)=T3;
  g1(91,90)=1;
  g1(92,16)=T3;
  g1(92,91)=1;
  g1(93,91)=T3;
  g1(93,92)=1;
  g1(94,92)=T3;
  g1(94,93)=1;
  g1(95,93)=T3;
  g1(95,94)=1;
  g1(96,94)=T3;
  g1(96,95)=1;
  g1(97,95)=T3;
  g1(97,96)=1;
  g1(98,96)=T3;
  g1(98,97)=1;
  g1(99,97)=T3;
  g1(99,98)=1;
  g1(100,98)=T3;
  g1(100,99)=1;
  g1(101,99)=T3;
  g1(101,100)=1;
  g1(102,100)=T3;
  g1(102,101)=1;
  g1(103,101)=T3;
  g1(103,102)=1;
  g1(104,102)=T3;
  g1(104,103)=1;
  g1(105,103)=T3;
  g1(105,104)=1;
  g1(106,104)=T3;
  g1(106,105)=1;
  g1(107,105)=T3;
  g1(107,106)=1;
  g1(108,106)=T3;
  g1(108,107)=1;
  g1(109,107)=T3;
  g1(109,108)=1;
  g1(110,108)=T3;
  g1(110,109)=1;
  g1(111,109)=T3;
  g1(111,110)=1;
  g1(112,110)=T3;
  g1(112,111)=1;
  g1(113,111)=T3;
  g1(113,112)=1;
  g1(114,112)=T3;
  g1(114,113)=1;
  g1(115,113)=T3;
  g1(115,114)=1;
  g1(116,114)=T3;
  g1(116,115)=1;
  g1(117,115)=T3;
  g1(117,116)=1;
  g1(118,116)=T3;
  g1(118,117)=1;
  g1(119,117)=T3;
  g1(119,118)=1;
  g1(120,118)=T3;
  g1(120,119)=1;
  g1(121,119)=T3;
  g1(121,120)=1;
  g1(122,120)=T3;
  g1(122,121)=1;
  g1(123,121)=T3;
  g1(123,122)=1;
  g1(124,122)=T3;
  g1(124,123)=1;
  g1(125,123)=T3;
  g1(125,124)=1;
  g1(126,124)=T3;
  g1(126,125)=1;
  g1(127,125)=T3;
  g1(127,126)=1;
  g1(128,126)=T3;
  g1(128,127)=1;
  g1(129,127)=T3;
  g1(129,128)=1;
  g1(130,128)=T3;
  g1(130,129)=1;
  g1(131,129)=T3;
  g1(131,130)=1;
  g1(132,1)=(-((1-params(1))/(y(1)-y(1)*params(1))-(-(params(20)*params(5)*(1-params(1)*params(2))*params(12)*params(5)*(1-params(1))))/(params(5)*y(1)*(1-params(1))*params(12)*params(5)*y(1)*(1-params(1))*params(12))/(1+params(6))));
  g1(132,20)=(-(params(4)*(-y(22))/(y(20)*y(20))/T22));
  g1(132,22)=(-(params(4)*1/y(20)/T22));
  g1(132,131)=1;
  g1(133,131)=(-(1+params(2)));
  g1(133,132)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],133,17689);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],133,2352637);
end
end
end
end
