max_heights_no_mass =[0.00830851576453556;0.0544376703915253;0.105801221915264;0.156733943057817;0.206733835251986;0.255943470152299;0.304382245441986;0.352432139756579;0.397889100655156;0.438489407609814;0.473879642345353;0.504353282402110;0.530744355983196;0.552846250607336;0.570366818228436;0.585662653177588;0.595687407275427;0.603205452245368;0.610697278226758;0.611923901703388;0.611679474251870;0.610679709869594;0.608143814451447;0.602010599028620;0.591495353684537;0.583647713815251;0.574901560376373;0.565355336672831;0.555098665211714;0.541901986718290;0.525478376481088;0.512005665356095;0.499178929775312;0.486152131100169;0.472945225825622;0.459569023175622;0.446027712616749;0.432319110186760;0.418388701671084;0.404407303552875;0.390206203811644;0.375850863359870;0.361355728261443;0.346742556600245;0.332037852279351;0.317272954361738;0.302482642545459;0.288557728051908;0.274998973631035;0.259704019334866;0.244677128969593;0.230284779221980;0.215677164080427;0.201526513683053;0.187684991481334;0.174239031894008;0.162189710770129;0.149497204204158;0.136827667702842;0.124518049060257;0.112547201863401;0.100885288022787;0.0577799500544462;0.0586175707853511;0.0594478185705252;0.0602696075681364;0.0610819239995006;0.0618838208983014;0.0626744136238872;0.0634528758822355;0.0642184361540361;0.0649703744672206;0.0657080194656212;0.0664307087944827;0.0671379343105236;0.0678291184909845;0.0685037971848256;0.0691614316222704;0.0698016302048052]
max_heights_with_mass = [0.243750763940343;0.283207223879392;0.319557203108224;0.352622883776024;0.382024557184375;0.407979473938731;0.430548025058278;0.449933234871425;0.466299708722265;0.479805776993207;0.490462927195320;0.499173037819728;0.504862428821477;0.509017679713779;0.510971338171835;0.511283284325734;0.510606346573238;0.507295727258118;0.503117073921786;0.498846520154363;0.493553959554492;0.484789058651260;0.477538472885227;0.469943890675421;0.461863666776726;0.453350540995443;0.440873509695920;0.431199200707334;0.421415155668732;0.411437110062095;0.401286600104786;0.390980886177811;0.380533351653185;0.369954488761614;0.359252563317410;0.348434000396211;0.337506178038454;0.326473517165574;0.315342242149203;0.304119180269963;0.292810720738617;0.281445569367173;0.271062118576116;0.259354553468294;0.247593568508363;0.236192568718280;0.224980416053398;0.212799113573009;0.200873385321275;0.188868960293419;0.176804680526242;0.164842039174920;0.152999224667579;0.141286842203240;0.129766740232652;0.118421176350614;0.0339607249731624;0.0343899031760770;0.0348167191353781;0.0352406513326175;0.0356612129044604;0.0360779485089823;0.0364904319560285;0.0368982641704778;0.0373010713838281;0.0376985035109967;0.0380902326806628;0.0384759519002348;0.0388553738374983;0.0392282297050356;0.0395942682351520]
gear_ratios_no_mass = 1:0.5:40;
gear_ratios_mass = 5:0.5:40;

heights1 = max_heights_no_mass(9:end);
heights2 = max_heights_with_mass;

mass1 = 0;
mass2 = 3;

% try writing to csv

% format: each row: (mass, then heights)
massHeight1 = [mass1 heights1.'];
massHeight2 = [mass2 heights2.'];

data = [massHeight1; massHeight2];

writematrix(data, 'M.xls')



