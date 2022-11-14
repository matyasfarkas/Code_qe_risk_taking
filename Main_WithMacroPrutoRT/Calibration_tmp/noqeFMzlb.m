%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'noqeFMzlb';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = 'z_shk';
M_.exo_names_tex = 'z\_shk';
M_.exo_names_long = 'z_shk';
M_.exo_names = char(M_.exo_names, 'mp_shk');
M_.exo_names_tex = char(M_.exo_names_tex, 'mp\_shk');
M_.exo_names_long = char(M_.exo_names_long, 'mp_shk');
M_.exo_names = char(M_.exo_names, 'bk_shk');
M_.exo_names_tex = char(M_.exo_names_tex, 'bk\_shk');
M_.exo_names_long = char(M_.exo_names_long, 'bk_shk');
M_.exo_names = char(M_.exo_names, 'qe_shk');
M_.exo_names_tex = char(M_.exo_names_tex, 'qe\_shk');
M_.exo_names_long = char(M_.exo_names_long, 'qe_shk');
M_.endo_names = 'ch';
M_.endo_names_tex = 'ch';
M_.endo_names_long = 'ch';
M_.endo_names = char(M_.endo_names, 'Ce');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ce');
M_.endo_names_long = char(M_.endo_names_long, 'Ce');
M_.endo_names = char(M_.endo_names, 'Cb');
M_.endo_names_tex = char(M_.endo_names_tex, 'Cb');
M_.endo_names_long = char(M_.endo_names_long, 'Cb');
M_.endo_names = char(M_.endo_names, 'H');
M_.endo_names_tex = char(M_.endo_names_tex, 'H');
M_.endo_names_long = char(M_.endo_names_long, 'H');
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names_long = char(M_.endo_names_long, 'Y');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names_long = char(M_.endo_names_long, 'I');
M_.endo_names = char(M_.endo_names, 'bigN');
M_.endo_names_tex = char(M_.endo_names_tex, 'bigN');
M_.endo_names_long = char(M_.endo_names_long, 'bigN');
M_.endo_names = char(M_.endo_names, 'bigA');
M_.endo_names_tex = char(M_.endo_names_tex, 'bigA');
M_.endo_names_long = char(M_.endo_names_long, 'bigA');
M_.endo_names = char(M_.endo_names, 'smallb');
M_.endo_names_tex = char(M_.endo_names_tex, 'smallb');
M_.endo_names_long = char(M_.endo_names_long, 'smallb');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'rk');
M_.endo_names_tex = char(M_.endo_names_tex, 'rk');
M_.endo_names_long = char(M_.endo_names_long, 'rk');
M_.endo_names = char(M_.endo_names, 'infl');
M_.endo_names_tex = char(M_.endo_names_tex, 'infl');
M_.endo_names_long = char(M_.endo_names_long, 'infl');
M_.endo_names = char(M_.endo_names, 'w_h');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_h');
M_.endo_names_long = char(M_.endo_names_long, 'w_h');
M_.endo_names = char(M_.endo_names, 'w_e');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_e');
M_.endo_names_long = char(M_.endo_names_long, 'w_e');
M_.endo_names = char(M_.endo_names, 'w_b');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_b');
M_.endo_names_long = char(M_.endo_names_long, 'w_b');
M_.endo_names = char(M_.endo_names, 'Rd');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rd');
M_.endo_names_long = char(M_.endo_names_long, 'Rd');
M_.endo_names = char(M_.endo_names, 'Ra');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ra');
M_.endo_names_long = char(M_.endo_names_long, 'Ra');
M_.endo_names = char(M_.endo_names, 'lam');
M_.endo_names_tex = char(M_.endo_names_tex, 'lam');
M_.endo_names_long = char(M_.endo_names_long, 'lam');
M_.endo_names = char(M_.endo_names, 'mu');
M_.endo_names_tex = char(M_.endo_names_tex, 'mu');
M_.endo_names_long = char(M_.endo_names_long, 'mu');
M_.endo_names = char(M_.endo_names, 'p');
M_.endo_names_tex = char(M_.endo_names_tex, 'p');
M_.endo_names_long = char(M_.endo_names_long, 'p');
M_.endo_names = char(M_.endo_names, 'TL');
M_.endo_names_tex = char(M_.endo_names_tex, 'TL');
M_.endo_names_long = char(M_.endo_names_long, 'TL');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 's');
M_.endo_names_tex = char(M_.endo_names_tex, 's');
M_.endo_names_long = char(M_.endo_names_long, 's');
M_.endo_names = char(M_.endo_names, 'mgrowth');
M_.endo_names_tex = char(M_.endo_names_tex, 'mgrowth');
M_.endo_names_long = char(M_.endo_names_long, 'mgrowth');
M_.endo_names = char(M_.endo_names, 'expinfl');
M_.endo_names_tex = char(M_.endo_names_tex, 'expinfl');
M_.endo_names_long = char(M_.endo_names_long, 'expinfl');
M_.endo_names = char(M_.endo_names, 'smalld');
M_.endo_names_tex = char(M_.endo_names_tex, 'smalld');
M_.endo_names_long = char(M_.endo_names_long, 'smalld');
M_.endo_names = char(M_.endo_names, 'gY');
M_.endo_names_tex = char(M_.endo_names_tex, 'gY');
M_.endo_names_long = char(M_.endo_names_long, 'gY');
M_.endo_names = char(M_.endo_names, 'totC');
M_.endo_names_tex = char(M_.endo_names_tex, 'totC');
M_.endo_names_long = char(M_.endo_names_long, 'totC');
M_.endo_names = char(M_.endo_names, 'u');
M_.endo_names_tex = char(M_.endo_names_tex, 'u');
M_.endo_names_long = char(M_.endo_names_long, 'u');
M_.endo_names = char(M_.endo_names, 'keff');
M_.endo_names_tex = char(M_.endo_names_tex, 'keff');
M_.endo_names_long = char(M_.endo_names_long, 'keff');
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, 'K');
M_.endo_names_long = char(M_.endo_names_long, 'K');
M_.endo_names = char(M_.endo_names, 'Kb');
M_.endo_names_tex = char(M_.endo_names_tex, 'Kb');
M_.endo_names_long = char(M_.endo_names_long, 'Kb');
M_.endo_names = char(M_.endo_names, 'Ke');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ke');
M_.endo_names_long = char(M_.endo_names_long, 'Ke');
M_.endo_names = char(M_.endo_names, 'ptilde');
M_.endo_names_tex = char(M_.endo_names_tex, 'ptilde');
M_.endo_names_long = char(M_.endo_names_long, 'ptilde');
M_.endo_names = char(M_.endo_names, 'nump');
M_.endo_names_tex = char(M_.endo_names_tex, 'nump');
M_.endo_names_long = char(M_.endo_names_long, 'nump');
M_.endo_names = char(M_.endo_names, 'denp');
M_.endo_names_tex = char(M_.endo_names_tex, 'denp');
M_.endo_names_long = char(M_.endo_names_long, 'denp');
M_.endo_names = char(M_.endo_names, 'wtilde');
M_.endo_names_tex = char(M_.endo_names_tex, 'wtilde');
M_.endo_names_long = char(M_.endo_names_long, 'wtilde');
M_.endo_names = char(M_.endo_names, 'numw');
M_.endo_names_tex = char(M_.endo_names_tex, 'numw');
M_.endo_names_long = char(M_.endo_names_long, 'numw');
M_.endo_names = char(M_.endo_names, 'denw');
M_.endo_names_tex = char(M_.endo_names_tex, 'denw');
M_.endo_names_long = char(M_.endo_names_long, 'denw');
M_.endo_names = char(M_.endo_names, 'lz');
M_.endo_names_tex = char(M_.endo_names_tex, 'lz');
M_.endo_names_long = char(M_.endo_names_long, 'lz');
M_.endo_names = char(M_.endo_names, 'lmp');
M_.endo_names_tex = char(M_.endo_names_tex, 'lmp');
M_.endo_names_long = char(M_.endo_names_long, 'lmp');
M_.endo_names = char(M_.endo_names, 'lbk');
M_.endo_names_tex = char(M_.endo_names_tex, 'lbk');
M_.endo_names_long = char(M_.endo_names_long, 'lbk');
M_.endo_names = char(M_.endo_names, 'log_y');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_y');
M_.endo_names_long = char(M_.endo_names_long, 'log_y');
M_.endo_names = char(M_.endo_names, 'log_I');
M_.endo_names_tex = char(M_.endo_names_tex, 'log\_I');
M_.endo_names_long = char(M_.endo_names_long, 'log_I');
M_.endo_names = char(M_.endo_names, 'gammag');
M_.endo_names_tex = char(M_.endo_names_tex, 'gammag');
M_.endo_names_long = char(M_.endo_names_long, 'gammag');
M_.endo_names = char(M_.endo_names, 'alpha');
M_.endo_names_tex = char(M_.endo_names_tex, 'alpha');
M_.endo_names_long = char(M_.endo_names_long, 'alpha');
M_.endo_names = char(M_.endo_names, 'rnot');
M_.endo_names_tex = char(M_.endo_names_tex, 'rnot');
M_.endo_names_long = char(M_.endo_names_long, 'rnot');
M_.endo_names = char(M_.endo_names, 'CBBS');
M_.endo_names_tex = char(M_.endo_names_tex, 'CBBS');
M_.endo_names_long = char(M_.endo_names_long, 'CBBS');
M_.endo_names = char(M_.endo_names, 'lqe');
M_.endo_names_tex = char(M_.endo_names_tex, 'lqe');
M_.endo_names_long = char(M_.endo_names_long, 'lqe');
M_.endo_names = char(M_.endo_names, 'inflexp1');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp1');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp1');
M_.endo_names = char(M_.endo_names, 'inflexp2');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp2');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp2');
M_.endo_names = char(M_.endo_names, 'inflexp3');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp3');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp3');
M_.endo_names = char(M_.endo_names, 'inflexp4');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp4');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp4');
M_.endo_names = char(M_.endo_names, 'inflexp5');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp5');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp5');
M_.endo_names = char(M_.endo_names, 'inflexp6');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp6');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp6');
M_.endo_names = char(M_.endo_names, 'inflexp7');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp7');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp7');
M_.endo_names = char(M_.endo_names, 'inflexp8');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp8');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp8');
M_.endo_names = char(M_.endo_names, 'inflexp9');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp9');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp9');
M_.endo_names = char(M_.endo_names, 'inflexp10');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp10');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp10');
M_.endo_names = char(M_.endo_names, 'inflexp11');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp11');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp11');
M_.endo_names = char(M_.endo_names, 'inflexp12');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp12');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp12');
M_.endo_names = char(M_.endo_names, 'inflexp13');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp13');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp13');
M_.endo_names = char(M_.endo_names, 'inflexp14');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp14');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp14');
M_.endo_names = char(M_.endo_names, 'inflexp15');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp15');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp15');
M_.endo_names = char(M_.endo_names, 'inflexp16');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp16');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp16');
M_.endo_names = char(M_.endo_names, 'inflexp17');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp17');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp17');
M_.endo_names = char(M_.endo_names, 'inflexp18');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp18');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp18');
M_.endo_names = char(M_.endo_names, 'inflexp19');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp19');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp19');
M_.endo_names = char(M_.endo_names, 'inflexp20');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp20');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp20');
M_.endo_names = char(M_.endo_names, 'inflexp21');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp21');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp21');
M_.endo_names = char(M_.endo_names, 'inflexp22');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp22');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp22');
M_.endo_names = char(M_.endo_names, 'inflexp23');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp23');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp23');
M_.endo_names = char(M_.endo_names, 'inflexp24');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp24');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp24');
M_.endo_names = char(M_.endo_names, 'inflexp25');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp25');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp25');
M_.endo_names = char(M_.endo_names, 'inflexp26');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp26');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp26');
M_.endo_names = char(M_.endo_names, 'inflexp27');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp27');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp27');
M_.endo_names = char(M_.endo_names, 'inflexp28');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp28');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp28');
M_.endo_names = char(M_.endo_names, 'inflexp29');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp29');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp29');
M_.endo_names = char(M_.endo_names, 'inflexp30');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp30');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp30');
M_.endo_names = char(M_.endo_names, 'inflexp31');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp31');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp31');
M_.endo_names = char(M_.endo_names, 'inflexp32');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp32');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp32');
M_.endo_names = char(M_.endo_names, 'inflexp33');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp33');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp33');
M_.endo_names = char(M_.endo_names, 'inflexp34');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp34');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp34');
M_.endo_names = char(M_.endo_names, 'inflexp35');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp35');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp35');
M_.endo_names = char(M_.endo_names, 'inflexp36');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp36');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp36');
M_.endo_names = char(M_.endo_names, 'inflexp37');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp37');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp37');
M_.endo_names = char(M_.endo_names, 'inflexp38');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp38');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp38');
M_.endo_names = char(M_.endo_names, 'inflexp39');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp39');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp39');
M_.endo_names = char(M_.endo_names, 'inflexp40');
M_.endo_names_tex = char(M_.endo_names_tex, 'inflexp40');
M_.endo_names_long = char(M_.endo_names_long, 'inflexp40');
M_.endo_names = char(M_.endo_names, 'Rexp1');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp1');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp1');
M_.endo_names = char(M_.endo_names, 'Rexp2');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp2');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp2');
M_.endo_names = char(M_.endo_names, 'Rexp3');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp3');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp3');
M_.endo_names = char(M_.endo_names, 'Rexp4');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp4');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp4');
M_.endo_names = char(M_.endo_names, 'Rexp5');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp5');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp5');
M_.endo_names = char(M_.endo_names, 'Rexp6');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp6');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp6');
M_.endo_names = char(M_.endo_names, 'Rexp7');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp7');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp7');
M_.endo_names = char(M_.endo_names, 'Rexp8');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp8');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp8');
M_.endo_names = char(M_.endo_names, 'Rexp9');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp9');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp9');
M_.endo_names = char(M_.endo_names, 'Rexp10');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp10');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp10');
M_.endo_names = char(M_.endo_names, 'Rexp11');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp11');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp11');
M_.endo_names = char(M_.endo_names, 'Rexp12');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp12');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp12');
M_.endo_names = char(M_.endo_names, 'Rexp13');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp13');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp13');
M_.endo_names = char(M_.endo_names, 'Rexp14');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp14');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp14');
M_.endo_names = char(M_.endo_names, 'Rexp15');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp15');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp15');
M_.endo_names = char(M_.endo_names, 'Rexp16');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp16');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp16');
M_.endo_names = char(M_.endo_names, 'Rexp17');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp17');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp17');
M_.endo_names = char(M_.endo_names, 'Rexp18');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp18');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp18');
M_.endo_names = char(M_.endo_names, 'Rexp19');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp19');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp19');
M_.endo_names = char(M_.endo_names, 'Rexp20');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp20');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp20');
M_.endo_names = char(M_.endo_names, 'Rexp21');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp21');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp21');
M_.endo_names = char(M_.endo_names, 'Rexp22');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp22');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp22');
M_.endo_names = char(M_.endo_names, 'Rexp23');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp23');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp23');
M_.endo_names = char(M_.endo_names, 'Rexp24');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp24');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp24');
M_.endo_names = char(M_.endo_names, 'Rexp25');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp25');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp25');
M_.endo_names = char(M_.endo_names, 'Rexp26');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp26');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp26');
M_.endo_names = char(M_.endo_names, 'Rexp27');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp27');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp27');
M_.endo_names = char(M_.endo_names, 'Rexp28');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp28');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp28');
M_.endo_names = char(M_.endo_names, 'Rexp29');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp29');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp29');
M_.endo_names = char(M_.endo_names, 'Rexp30');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp30');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp30');
M_.endo_names = char(M_.endo_names, 'Rexp31');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp31');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp31');
M_.endo_names = char(M_.endo_names, 'Rexp32');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp32');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp32');
M_.endo_names = char(M_.endo_names, 'Rexp33');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp33');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp33');
M_.endo_names = char(M_.endo_names, 'Rexp34');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp34');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp34');
M_.endo_names = char(M_.endo_names, 'Rexp35');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp35');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp35');
M_.endo_names = char(M_.endo_names, 'Rexp36');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp36');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp36');
M_.endo_names = char(M_.endo_names, 'Rexp37');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp37');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp37');
M_.endo_names = char(M_.endo_names, 'Rexp38');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp38');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp38');
M_.endo_names = char(M_.endo_names, 'Rexp39');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp39');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp39');
M_.endo_names = char(M_.endo_names, 'Rexp40');
M_.endo_names_tex = char(M_.endo_names_tex, 'Rexp40');
M_.endo_names_long = char(M_.endo_names_long, 'Rexp40');
M_.endo_names = char(M_.endo_names, 'BVaR');
M_.endo_names_tex = char(M_.endo_names_tex, 'BVaR');
M_.endo_names_long = char(M_.endo_names_long, 'BVaR');
M_.param_names = 'habit';
M_.param_names_tex = 'habit';
M_.param_names_long = 'habit';
M_.param_names = char(M_.param_names, 'bet');
M_.param_names_tex = char(M_.param_names_tex, 'bet');
M_.param_names_long = char(M_.param_names_long, 'bet');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'reta');
M_.param_names_tex = char(M_.param_names_tex, 'reta');
M_.param_names_long = char(M_.param_names_long, 'reta');
M_.param_names = char(M_.param_names, 'psi_l');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_l');
M_.param_names_long = char(M_.param_names_long, 'psi_l');
M_.param_names = char(M_.param_names, 'elas_l');
M_.param_names_tex = char(M_.param_names_tex, 'elas\_l');
M_.param_names_long = char(M_.param_names_long, 'elas_l');
M_.param_names = char(M_.param_names, 'chi1');
M_.param_names_tex = char(M_.param_names_tex, 'chi1');
M_.param_names_long = char(M_.param_names_long, 'chi1');
M_.param_names = char(M_.param_names, 'chi2');
M_.param_names_tex = char(M_.param_names_tex, 'chi2');
M_.param_names_long = char(M_.param_names_long, 'chi2');
M_.param_names = char(M_.param_names, 'pi_ss');
M_.param_names_tex = char(M_.param_names_tex, 'pi\_ss');
M_.param_names_long = char(M_.param_names_long, 'pi_ss');
M_.param_names = char(M_.param_names, 'xi_w');
M_.param_names_tex = char(M_.param_names_tex, 'xi\_w');
M_.param_names_long = char(M_.param_names_long, 'xi_w');
M_.param_names = char(M_.param_names, 'phi_w');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_w');
M_.param_names_long = char(M_.param_names_long, 'phi_w');
M_.param_names = char(M_.param_names, 'mark_w');
M_.param_names_tex = char(M_.param_names_tex, 'mark\_w');
M_.param_names_long = char(M_.param_names_long, 'mark_w');
M_.param_names = char(M_.param_names, 'eta_h');
M_.param_names_tex = char(M_.param_names_tex, 'eta\_h');
M_.param_names_long = char(M_.param_names_long, 'eta_h');
M_.param_names = char(M_.param_names, 'eta_e');
M_.param_names_tex = char(M_.param_names_tex, 'eta\_e');
M_.param_names_long = char(M_.param_names_long, 'eta_e');
M_.param_names = char(M_.param_names, 'eta_b');
M_.param_names_tex = char(M_.param_names_tex, 'eta\_b');
M_.param_names_long = char(M_.param_names_long, 'eta_b');
M_.param_names = char(M_.param_names, 'xi_p');
M_.param_names_tex = char(M_.param_names_tex, 'xi\_p');
M_.param_names_long = char(M_.param_names_long, 'xi_p');
M_.param_names = char(M_.param_names, 'phi_p');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_p');
M_.param_names_long = char(M_.param_names_long, 'phi_p');
M_.param_names = char(M_.param_names, 'mark_p');
M_.param_names_tex = char(M_.param_names_tex, 'mark\_p');
M_.param_names_long = char(M_.param_names_long, 'mark_p');
M_.param_names = char(M_.param_names, 'theta_k');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_k');
M_.param_names_long = char(M_.param_names_long, 'theta_k');
M_.param_names = char(M_.param_names, 'theta_h');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_h');
M_.param_names_long = char(M_.param_names_long, 'theta_h');
M_.param_names = char(M_.param_names, 'theta_e');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_e');
M_.param_names_long = char(M_.param_names_long, 'theta_e');
M_.param_names = char(M_.param_names, 'theta_b');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_b');
M_.param_names_long = char(M_.param_names_long, 'theta_b');
M_.param_names = char(M_.param_names, 'bigR');
M_.param_names_tex = char(M_.param_names_tex, 'bigR');
M_.param_names_long = char(M_.param_names_long, 'bigR');
M_.param_names = char(M_.param_names, 'delalpha');
M_.param_names_tex = char(M_.param_names_tex, 'delalpha');
M_.param_names_long = char(M_.param_names_long, 'delalpha');
M_.param_names = char(M_.param_names, 'tau_b');
M_.param_names_tex = char(M_.param_names_tex, 'tau\_b');
M_.param_names_long = char(M_.param_names_long, 'tau_b');
M_.param_names = char(M_.param_names, 'tau_e');
M_.param_names_tex = char(M_.param_names_tex, 'tau\_e');
M_.param_names_long = char(M_.param_names_long, 'tau_e');
M_.param_names = char(M_.param_names, 'lam_r');
M_.param_names_tex = char(M_.param_names_tex, 'lam\_r');
M_.param_names_long = char(M_.param_names_long, 'lam_r');
M_.param_names = char(M_.param_names, 'lam_pi');
M_.param_names_tex = char(M_.param_names_tex, 'lam\_pi');
M_.param_names_long = char(M_.param_names_long, 'lam_pi');
M_.param_names = char(M_.param_names, 'lam_y');
M_.param_names_tex = char(M_.param_names_tex, 'lam\_y');
M_.param_names_long = char(M_.param_names_long, 'lam_y');
M_.param_names = char(M_.param_names, 'bigtheta');
M_.param_names_tex = char(M_.param_names_tex, 'bigtheta');
M_.param_names_long = char(M_.param_names_long, 'bigtheta');
M_.param_names = char(M_.param_names, 'bby');
M_.param_names_tex = char(M_.param_names_tex, 'bby');
M_.param_names_long = char(M_.param_names_long, 'bby');
M_.param_names = char(M_.param_names, 'nu');
M_.param_names_tex = char(M_.param_names_tex, 'nu');
M_.param_names_long = char(M_.param_names_long, 'nu');
M_.param_names = char(M_.param_names, 'rhoz');
M_.param_names_tex = char(M_.param_names_tex, 'rhoz');
M_.param_names_long = char(M_.param_names_long, 'rhoz');
M_.param_names = char(M_.param_names, 'rhomp');
M_.param_names_tex = char(M_.param_names_tex, 'rhomp');
M_.param_names_long = char(M_.param_names_long, 'rhomp');
M_.param_names = char(M_.param_names, 'rhobk');
M_.param_names_tex = char(M_.param_names_tex, 'rhobk');
M_.param_names_long = char(M_.param_names_long, 'rhobk');
M_.param_names = char(M_.param_names, 'sigmaz');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaz');
M_.param_names_long = char(M_.param_names_long, 'sigmaz');
M_.param_names = char(M_.param_names, 'sigmamp');
M_.param_names_tex = char(M_.param_names_tex, 'sigmamp');
M_.param_names_long = char(M_.param_names_long, 'sigmamp');
M_.param_names = char(M_.param_names, 'sigmabk');
M_.param_names_tex = char(M_.param_names_tex, 'sigmabk');
M_.param_names_long = char(M_.param_names_long, 'sigmabk');
M_.param_names = char(M_.param_names, 'alpha_ss');
M_.param_names_tex = char(M_.param_names_tex, 'alpha\_ss');
M_.param_names_long = char(M_.param_names_long, 'alpha_ss');
M_.param_names = char(M_.param_names, 'mu_ss');
M_.param_names_tex = char(M_.param_names_tex, 'mu\_ss');
M_.param_names_long = char(M_.param_names_long, 'mu_ss');
M_.param_names = char(M_.param_names, 'omega');
M_.param_names_tex = char(M_.param_names_tex, 'omega');
M_.param_names_long = char(M_.param_names_long, 'omega');
M_.param_names = char(M_.param_names, 'varsigma');
M_.param_names_tex = char(M_.param_names_tex, 'varsigma');
M_.param_names_long = char(M_.param_names_long, 'varsigma');
M_.param_names = char(M_.param_names, 'epsb');
M_.param_names_tex = char(M_.param_names_tex, 'epsb');
M_.param_names_long = char(M_.param_names_long, 'epsb');
M_.param_names = char(M_.param_names, 'Btopbar');
M_.param_names_tex = char(M_.param_names_tex, 'Btopbar');
M_.param_names_long = char(M_.param_names_long, 'Btopbar');
M_.param_names = char(M_.param_names, 'Blowbar');
M_.param_names_tex = char(M_.param_names_tex, 'Blowbar');
M_.param_names_long = char(M_.param_names_long, 'Blowbar');
M_.param_names = char(M_.param_names, 'Chi');
M_.param_names_tex = char(M_.param_names_tex, 'Chi');
M_.param_names_long = char(M_.param_names_long, 'Chi');
M_.param_names = char(M_.param_names, 'CBBSp');
M_.param_names_tex = char(M_.param_names_tex, 'CBBSp');
M_.param_names_long = char(M_.param_names_long, 'CBBSp');
M_.param_names = char(M_.param_names, 'sigmaqe');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaqe');
M_.param_names_long = char(M_.param_names_long, 'sigmaqe');
M_.param_names = char(M_.param_names, 'rhoqe');
M_.param_names_tex = char(M_.param_names_tex, 'rhoqe');
M_.param_names_long = char(M_.param_names_long, 'rhoqe');
M_.exo_det_nbr = 0;
M_.exo_nbr = 4;
M_.endo_nbr = 130;
M_.param_nbr = 49;
M_.orig_endo_nbr = 130;
M_.aux_vars = [];
M_.Sigma_e = zeros(4, 4);
M_.Correlation_matrix = eye(4, 4);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('noqeFMzlb_static');
erase_compiled_function('noqeFMzlb_dynamic');
M_.lead_lag_incidence = [
 1 14 144;
 0 15 0;
 0 16 0;
 0 17 0;
 2 18 0;
 0 19 0;
 0 20 0;
 0 21 0;
 0 22 0;
 0 23 145;
 0 24 146;
 3 25 147;
 4 26 0;
 0 27 0;
 0 28 0;
 0 29 148;
 0 30 0;
 0 31 149;
 0 32 0;
 5 33 0;
 0 34 0;
 0 35 0;
 0 36 0;
 0 37 0;
 0 38 0;
 0 39 0;
 0 40 0;
 0 41 0;
 0 42 150;
 0 43 0;
 6 44 0;
 7 45 0;
 8 46 0;
 0 47 0;
 0 48 151;
 0 49 152;
 0 50 0;
 0 51 153;
 0 52 154;
 9 53 0;
 10 54 0;
 11 55 0;
 0 56 0;
 0 57 0;
 0 58 0;
 0 59 0;
 12 60 0;
 0 61 0;
 13 62 0;
 0 63 155;
 0 64 156;
 0 65 157;
 0 66 158;
 0 67 159;
 0 68 160;
 0 69 161;
 0 70 162;
 0 71 163;
 0 72 164;
 0 73 165;
 0 74 166;
 0 75 167;
 0 76 168;
 0 77 169;
 0 78 170;
 0 79 171;
 0 80 172;
 0 81 173;
 0 82 174;
 0 83 175;
 0 84 176;
 0 85 177;
 0 86 178;
 0 87 179;
 0 88 180;
 0 89 181;
 0 90 182;
 0 91 183;
 0 92 184;
 0 93 185;
 0 94 186;
 0 95 187;
 0 96 188;
 0 97 189;
 0 98 190;
 0 99 191;
 0 100 192;
 0 101 193;
 0 102 0;
 0 103 194;
 0 104 195;
 0 105 196;
 0 106 197;
 0 107 198;
 0 108 199;
 0 109 200;
 0 110 201;
 0 111 202;
 0 112 203;
 0 113 204;
 0 114 205;
 0 115 206;
 0 116 207;
 0 117 208;
 0 118 209;
 0 119 210;
 0 120 211;
 0 121 212;
 0 122 213;
 0 123 214;
 0 124 215;
 0 125 216;
 0 126 217;
 0 127 218;
 0 128 219;
 0 129 220;
 0 130 221;
 0 131 222;
 0 132 223;
 0 133 224;
 0 134 225;
 0 135 226;
 0 136 227;
 0 137 228;
 0 138 229;
 0 139 230;
 0 140 231;
 0 141 232;
 0 142 0;
 0 143 0;]';
M_.nstatic = 30;
M_.nfwrd   = 87;
M_.npred   = 11;
M_.nboth   = 2;
M_.nsfwrd   = 89;
M_.nspred   = 13;
M_.ndynamic   = 100;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:4];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(130, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(4, 1);
M_.params = NaN(49, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 364;
M_.NNZDerivatives(2) = 335;
M_.NNZDerivatives(3) = -1;
M_.params( 48 ) = 1.00000000000000;
sigmaqe = M_.params( 48 );
M_.params( 49 ) = 0.9000;
rhoqe = M_.params( 49 );
M_.params( 1 ) = 0.6500;
habit = M_.params( 1 );
M_.params( 2 ) = 0.9950;
bet = M_.params( 2 );
M_.params( 3 ) = 0.0200;
delta = M_.params( 3 );
M_.params( 4 ) = 0.00183221900000;
reta = M_.params( 4 );
M_.params( 5 ) = 0.4550;
psi_l = M_.params( 5 );
M_.params( 6 ) = 0.4290;
elas_l = M_.params( 6 );
M_.params( 9 ) = 1.0049629;
pi_ss = M_.params( 9 );
M_.params( 27 ) = 0.9000;
lam_r = M_.params( 27 );
M_.params( 28 ) = 1.8000;
lam_pi = M_.params( 28 );
M_.params( 29 ) = 0.0000;
lam_y = M_.params( 29 );
Rd_temp = pi_ss/bet;
M_.params( 10 ) = 21.0;
xi_w = M_.params( 10 );
M_.params( 11 ) = 0.6400;
phi_w = M_.params( 11 );
M_.params( 16 ) = 6.0000;
xi_p = M_.params( 16 );
M_.params( 17 ) = 0.6000;
phi_p = M_.params( 17 );
M_.params( 18 ) = M_.params(16)/(M_.params(16)-1);
mark_p = M_.params( 18 );
M_.params( 12 ) = M_.params(10)/(M_.params(10)-1);
mark_w = M_.params( 12 );
M_.params( 13 ) = 0.9000;
eta_h = M_.params( 13 );
M_.params( 14 ) = 0.0700;
eta_e = M_.params( 14 );
M_.params( 15 ) = 1-M_.params(13)-M_.params(14);
eta_b = M_.params( 15 );
M_.params( 20 ) = 0.63999900000000;
theta_h = M_.params( 20 );
M_.params( 19 ) = 0.36000000000000;
theta_k = M_.params( 19 );
M_.params( 21 ) = (1.0-M_.params(20)-M_.params(19))/2.0;
theta_e = M_.params( 21 );
M_.params( 22 ) = 1.0-M_.params(20)-M_.params(19)-M_.params(21);
theta_b = M_.params( 22 );
M_.params( 24 ) = 0.35000000000000;
delalpha = M_.params( 24 );
M_.params( 25 ) = 0.70000000000000;
tau_b = M_.params( 25 );
M_.params( 26 ) = 0.70000000000000;
tau_e = M_.params( 26 );
M_.params( 23 ) = 1.00000000000000;
bigR = M_.params( 23 );
M_.params( 31 ) = 0.0000;
bby = M_.params( 31 );
M_.params( 32 ) = 1.00000000000000;
nu = M_.params( 32 );
M_.params( 33 ) = 0.95000000000000;
rhoz = M_.params( 33 );
M_.params( 36 ) = 0.00350000000000;
sigmaz = M_.params( 36 );
sigma_a = 0.50000000000000;
M_.params( 34 ) = 0.00000000000000;
rhomp = M_.params( 34 );
M_.params( 37 ) = 0.00160000000000;
sigmamp = M_.params( 37 );
M_.params( 35 ) = 0.90000000000000;
rhobk = M_.params( 35 );
M_.params( 38 ) = 0.90000000000000;
sigmabk = M_.params( 38 );
M_.params( 39 ) = 0.99300000000000;
alpha_ss = M_.params( 39 );
M_.params( 46 ) = 15.0000000000000;
Chi = M_.params( 46 );
M_.params( 43 ) = 10.0000000000000;
epsb = M_.params( 43 );
M_.params( 44 ) = 3.8046172500000;
Btopbar = M_.params( 44 );
M_.params( 45 ) = 0.0000;
Blowbar = M_.params( 45 );
M_.params( 41 ) = (-324.60000000000002);
omega = M_.params( 41 );
M_.params( 42 ) = 0.00000000000000;
varsigma = M_.params( 42 );
M_.params( 40 ) = 0.02500000000000;
mu_ss = M_.params( 40 );
M_.params( 47 ) = (-0.81000000000000);
CBBSp = M_.params( 47 );
alphatemp = alpha_ss;
mutemp = mu_ss;
smallbtemp = Btopbar*(1+Chi*mutemp)^(-epsb);
kk1 = delta*theta_k/(alphatemp*bigR*(1/bet-1+delta)); 
betatemp  = 1/bet;
infl_temp = pi_ss;
qtemp = kk1*(alphatemp*tau_b*mutemp*betatemp/delalpha-(1+mutemp+alphatemp*mutemp/(Rd_temp*delalpha)));
qtemp = qtemp/(alphatemp*smallbtemp*kk1/(Rd_temp*delalpha)-alphatemp*bigR*kk1/Rd_temp-theta_e-theta_b-bby-alphatemp*tau_e*smallbtemp*kk1*betatemp/delalpha);
Gtemp = 1+mutemp-(qtemp*alphatemp/Rd_temp)*(bigR-mutemp/(qtemp*delalpha)-smallbtemp/delalpha);
IY = kk1/qtemp;
KY = alphatemp*bigR*IY/delta; 
KeY = tau_e*alphatemp*smallbtemp*IY/delalpha;
KbY = tau_b*alphatemp*mutemp*IY/(qtemp*delalpha);
KhY = KY - KeY - KbY;
CeY = (1-tau_e)*alphatemp*smallbtemp*IY*qtemp/delalpha;
CbY = nu*(1-tau_b)*alphatemp*mutemp*IY/delalpha;
ChY = 1 + bby - CeY - CbY  - (1+mutemp)*IY;
essai = (1-bet*habit)*theta_h/((1-habit)*ChY*psi_l*mark_w);
smallh = essai^(1/(1+elas_l));
Htemp = smallh*eta_h^(xi_w/(xi_w-1));
kk2 = (KY^(1/(1-theta_k))) * (eta_e^(theta_e/(1-theta_k)))* (eta_b^(theta_b/(1-theta_k)))* mark_p^(1/(theta_k-1));
Ktemp = Htemp^(theta_h/(1-theta_k))* kk2;
Ytemp = (1/mark_p)* (Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^theta_b);
M_.params( 30 ) = (M_.params(18)-1)*Ytemp;
bigtheta = M_.params( 30 );
Itemp = IY*Ytemp;
rktemp = (1/mark_p)*theta_k*(Ktemp^(theta_k-1))*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^theta_b);
M_.params( 7 ) = rktemp;
chi1 = M_.params( 7 );
M_.params( 8 ) = sigma_a*M_.params(7);
chi2 = M_.params( 8 );
w_htemp = (1/mark_p)*theta_h*(Ktemp^theta_k)*(Htemp^(theta_h-1))*(eta_e^theta_e)*(eta_b^theta_b);
w_etemp = (1/mark_p)*theta_e*(Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^(theta_e-1))*(eta_b^theta_b);
w_btemp = (1/mark_p)*theta_b*(Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^(theta_b-1));
Ketemp = KeY*Ytemp;
Kbtemp = KbY*Ytemp;
Khtemp = KhY*Ytemp;
Cetemp = CeY*Ytemp;
Cbtemp = CbY*Ytemp;
chtemp = ChY*Ytemp/eta_h;
bigAtemp = eta_b*w_btemp+(rktemp+qtemp*(1-delta))*Kbtemp+bby*Ytemp;
bigNtemp = eta_e*w_etemp+(rktemp+qtemp*(1-delta))*Ketemp;
totCtemp = Cetemp+Cbtemp+eta_h*chtemp;
lamtemp = (1-bet*habit)/((1-habit)*(chtemp));
smalldtemp = alphatemp*qtemp*( bigR - smallbtemp/delalpha - mutemp/(qtemp*delalpha))*Itemp/(Rd_temp*eta_b);
ptemp = (Rd_temp-1)*lamtemp/ ( (Rd_temp-1)*lamtemp*eta_b*smalldtemp + eta_h*reta);
mctemp = reta*ptemp/((Rd_temp-1)*lamtemp);
Ratemp = alphatemp*qtemp*(mutemp/(qtemp*delalpha))*Itemp/bigAtemp;
TLtemp = Itemp-bigNtemp;
wtildetemp = w_htemp*infl_temp*eta_h^(-1/(1-xi_w));
numwtemp =  ( ( (w_htemp*infl_temp)^(xi_w*(1+elas_l)) * Htemp^(1+elas_l) )/(1-bet*phi_w) )^(1/(1+elas_l*xi_w));
denwtemp = ( (w_htemp^(xi_w)*infl_temp^(xi_w-1)*Htemp*lamtemp)/(1-bet*phi_w) )^(1/(1+elas_l*xi_w));
wtildetemp = (mark_w*psi_l)^(1/(1+xi_w*elas_l))*numwtemp/denwtemp;
stemp = 1/mark_p;
mgrowthtemp = infl_temp;
expinfltemp = infl_temp;
gYtemp = 1;
utemp = 1;
gammagtemp = (Itemp-bigNtemp )/bigAtemp;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(36)^2;
M_.Sigma_e(2, 2) = M_.params(37)^2;
M_.Sigma_e(3, 3) = M_.params(38)^2;
M_.Sigma_e(4, 4) = M_.params(48)^2;
M_.sigma_e_is_diagonal = 1;
save('noqeFMzlb_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('noqeFMzlb_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('noqeFMzlb_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('noqeFMzlb_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('noqeFMzlb_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
