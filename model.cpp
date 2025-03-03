// this file is distributed under 
// MIT license
#include <math_h/tabledata.h>
#include "model.h"
const MathTemplates::LinearInterpolation<> BC420_lambda = MathTemplates::Points<>{
    // Photon wavelength distribution
    // lambda(nm),  probability (a.u.)
    {359.686,0.686036},
    {361.518,1.66851},
    {363.089,2.79232},
    {364.66,3.91613},
    {365.707,5.18177},
    {367.016,6.44692},
    {367.801,7.85389},
    {368.586,8.97918},
    {370.157,11.6523},
    {370.681,13.0598},
    {371.204,14.4672},
    {371.728,15.7338},
    {372.251,17.1413},
    {373.822,23.3355},
    {374.084,24.7435},
    {374.346,26.1515},
    {375.393,33.0509},
    {375.916,34.4584},
    {376.178,35.8663},
    {376.963,42.2029},
    {377.225,43.6108},
    {378.796,54.1713},
    {379.058,55.5792},
    {380.105,63.0421},
    {380.366,64.45},
    {380.628,65.858},
    {381.675,73.0391},
    {382.199,74.4466},
    {382.461,75.8545},
    {383.246,81.346},
    {383.508,82.754},
    {384.031,84.1614},
    {384.555,85.4281},
    {385.079,86.8355},
    {385.602,88.243},
    {386.649,90.6354},
    {387.696,91.901},
    {388.743,93.3075},
    {390.576,94.29},
    {392.932,94.4264},
    {394.503,93.4376},
    {395.812,92.1675},
    {396.597,90.8984},
    {398.168,88.2194},
    {398.691,86.81},
    {399.215,85.4005},
    {400,84.1315},
    {400.524,82.722},
    {401.309,81.4529},
    {402.618,78.3519},
    {403.403,77.0828},
    {404.188,75.6729},
    {404.712,74.2635},
    {405.236,72.854},
    {405.497,71.4451},
    {407.33,66.6529},
    {407.853,65.2435},
    {408.639,63.9744},
    {408.901,62.5654},
    {409.424,61.156},
    {410.471,58.1963},
    {410.995,56.7869},
    {411.518,55.3774},
    {412.304,54.1083},
    {413.089,52.8393},
    {413.874,51.4293},
    {415.707,48.0456},
    {416.492,46.6357},
    {417.277,45.3666},
    {417.801,43.9572},
    {418.586,42.6881},
    {419.372,41.2782},
    {420.419,38.741},
    {421.204,37.4719},
    {421.99,36.062},
    {422.775,34.6521},
    {423.822,33.3825},
    {426.963,28.7287},
    {428.01,27.4591},
    {429.058,26.1896},
    {430.366,24.9195},
    {431.414,23.6499},
    {432.723,22.3799},
    {433.77,21.1103},
    {435.602,19.2759},
    {436.911,18.0058},
    {438.22,16.8766},
    {439.791,15.7469},
    {441.361,14.6172},
    {442.932,13.4875},
    {443.717,12.9226},
    {445.288,11.7929},
    {447.12,10.8035},
    {448.953,9.53248},
    {451.047,8.68348},
    {452.88,7.83497},
    {454.974,6.98596},
    {457.33,6.27731},
    {459.686,5.70951},
    {462.042,5.00086},
    {463.089,4.7172},
    {465.445,4.1494},
    {467.801,3.72244},
    {470.157,3.29548},
    {472.513,2.86852},
    {475.131,2.58191},
    {477.487,2.15495},
    {480.105,2.00919},
    {481.937,1.86491},
    {484.293,1.57879},
    {486.911,1.43303},
    {489.267,1.14692},
    {491.885,0.86031},
    {494.241,0.574196},
    {496.859,0.428434}
};
const MathTemplates::LinearInterpolation<> polyester_absorp = MathTemplates::Points<>{
    // absorption coefficient depending on photon wavelength
    // lambda(nm),  coef(mm^(-1))
    {348.842, 0.0083},
    {351.158, 0.0077},
    {354.247, 0.0069},
    {356.564, 0.0063},
    {358.88,  0.0057},
    {366.216, 0.0047},
    {374.324, 0.0032},
    {381.274, 0.0019},
    {388.996, 0.0012},
    {395.946, 0.0009},
    {413.707, 0.0005},
    {505.985, 0.0003},
    {597.876, 0.0002},
    {625.29,  0.0002}
};
const MathTemplates::LinearInterpolation<> Si_Photo_QE = MathTemplates::Points<>{
    // Efficiency of silicon photosensor depending on photon wavelength
    // lambda(nm), probability(n.d.)
    {140,	0.0},
    {160,	0.0},
    {180,	0.0},
    {200,	0.33},
    {340,	0.43},
    {380,	0.45},
    {400,	0.49},
    {420,	0.54},
    {440,	0.58},
    {460,	0.64},
    {480,	0.68},
    {500,	0.72},
    {520,	0.73},
    {540,	0.76},
    {560,	0.79},
    {580,	0.82},
    {600,	0.82},
    {620,	0.83},
    {640,	0.84},
    {660,	0.84},
    {680,	0.83},
    {700,	0.82},
    {720,	0.80},
    {740,	0.78},
    {760,	0.75},
    {780,	0.73},
    {800,	0.69},
    {820,	0.66},
    {840,	0.61},
    {860,	0.54},
    {880,	0.48},
    {900,	0.43},
    {920,	0.36},
    {940,	0.30}
};

const MathTemplates::LinearInterpolation<> tube_QE = MathTemplates::Points<>{
    // Efficiency of photomultiplier depending on photon wavelength
    // lambda(nm), probability(n.d.)
    {270,	0.1040038519},
    {280,	0.1457088571},
    {290,	0.1817626207},
    {300,	0.2040502667},
    {310,	0.220952},
    {320,	0.229768125},
    {330,	0.2348109091},
    {340,	0.2364971765},
    {350,	0.236716},
    {360,	0.2373497778},
    {370,	0.237392973},
    {380,	0.2358773684},
    {390,	0.234229641},
    {400,	0.2307547},
    {410,	0.226617561},
    {420,	0.2218153333},
    {430,	0.215408186},
    {440,	0.208103},
    {450,	0.1983724444},
    {460,	0.1871645217},
    {470,	0.1748611064},
    {480,	0.1644349074},
    {490,	0.1537881179},
    {500,	0.1422390222},
    {510,	0.1280928105},
    {520,	0.1118305128},
    {530,	0.0945025577},
    {540,	0.0788395062},
    {550,	0.065920404},
    {560,	0.0554703175},
    {570,	0.0463368421},
    {580,	0.0389103448},
    {590,	0.0318259661},
    {600,	0.0249963333},
    {610,	0.0188297049},
    {620,	0.013492},
    {630,	0.009176},
    {640,	0.0058764375},
    {650,	0.0034968},
    {660,	0.001984},
    {670,	0.0010660299},
    {680,	0.0005434118},
    {690,	0.0002695652},
    {700,	0.0001310857 }
};
