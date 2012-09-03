 /* **********************************************************
 * Header file defining what Markov processes are
 * Each process needs with associated with it the following:
 * 	- Data
 * 	- Tree
 *	- A calculation function
 * Each process can be used to calculate partial likelihoods up to
 * given nodes in the tree or get full likelihoods (usual calcs + mixture models)
 * ********************************************************* */

#ifndef PROCESS_HEADER
#define PROCESS_HEADER

#include "Leaphy.h"
#include "data.h"
#include "tree.h"

static int ProcID = 0;
static int QMatID = 0;

class CBaseProcess;

// Some definitions for THMM models
#define DEBUG_HMP_MODEL 1
#define HMP_MODEL_2_SIMPLE 0	// Forces the THMM model to mimic a standard model, useful for debugging.

// Some other definitions
enum ECodonEqm { cEQU,F1X4,F3X4,F64 };				// Defines the different types of codon frequencies

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino acid model information
const double dWAGVal[190] = {0.551571,0.509848,0.635346,0.738998,0.147304,5.429420,1.027040,0.528191,0.265256,0.0302949,0.908598,3.035500,1.543640,0.616783,0.0988179,1.582850,0.439157,0.947198,6.174160,0.021352,5.469470,1.416720,0.584665,1.125560,0.865584,0.306674,0.330052,0.567717,0.316954,2.137150,3.956290,0.930676,0.248972,4.294110,0.570025,0.249410,0.193335,0.186979,0.554236,0.039437,0.170135,0.113917,0.127395,0.0304501,0.138190,0.397915,0.497671,0.131528,0.0848047,0.384287,0.869489,0.154263,0.0613037,0.499462,3.170970,0.906265,5.351420,3.012010,0.479855,0.0740339,3.894900,2.584430,0.373558,0.890432,0.323832,0.257555,0.893496,0.683162,0.198221,0.103754,0.390482,1.545260,0.315124,0.174100,0.404141,4.257460,4.854020,0.934276,0.210494,0.102711,0.0961621,0.0467304,0.398020,0.0999208,0.0811339,0.049931,0.679371,1.059470,2.115170,0.088836,1.190630,1.438550,0.679489,0.195081,0.423984,0.109404,0.933372,0.682355,0.243570,0.696198,0.0999288,0.415844,0.556896,0.171329,0.161444,3.370790,1.224190,3.974230,1.071760,1.407660,1.028870,0.704939,1.341820,0.740169,0.319440,0.344739,0.967130,0.493905,0.545931,1.613280,2.121110,0.554413,2.030060,0.374866,0.512984,0.857928,0.822765,0.225833,0.473307,1.458160,0.326622,1.386980,1.516120,0.171903,0.795384,4.378020,0.113133,1.163920,0.0719167,0.129767,0.717070,0.215737,0.156557,0.336983,0.262569,0.212483,0.665309,0.137505,0.515706,1.529640,0.139405,0.523742,0.110864,0.240735,0.381533,1.086000,0.325711,0.543833,0.227710,0.196303,0.103604,3.873440,0.420170,0.398618,0.133264,0.428437,6.454280,0.216046,0.786993,0.291148,2.485390,2.006010,0.251849,0.196246,0.152335,1.002140,0.301281,0.588731,0.187247,0.118358,7.821300,1.800340,0.305434,2.058450,0.649892,0.314887,0.232739,1.388230,0.365369,0.314730};
const double dWAGFreq[20] = {0.0866279,0.043972,0.0390894,0.0570451,0.0193078,0.0367281,0.0580589,0.0832518,0.0244313,0.048466,0.086209,0.0620286,0.0195027,0.0384319,0.0457631,0.0695179,0.0610127,0.0143859,0.0352742,0.0708956};
const double dJTTVal[190] = {58,54,45,81,16,528,56,113,34,10,57,310,86,49,9,105,29,58,767,5,323,179,137,81,130,59,26,119,27,328,391,112,69,597,26,23,36,22,47,11,17,9,12,6,16,30,38,12,7,23,72,9,6,56,229,35,646,263,26,7,292,181,27,45,21,14,54,44,30,15,31,43,18,14,33,479,388,65,15,5,10,4,78,4,5,5,40,89,248,4,43,194,74,15,15,14,164,18,24,115,10,102,21,16,17,378,101,503,59,223,53,30,201,73,40,59,47,29,92,285,475,64,232,38,42,51,32,33,46,245,25,103,226,12,118,477,9,126,8,4,115,18,10,55,8,9,52,10,24,53,6,35,12,11,20,70,46,209,24,7,8,573,32,24,8,18,536,10,63,21,71,298,17,16,31,62,20,45,47,11,961,180,14,323,62,23,38,112,25,16};
const double dJTTFreq[20] = {7.674791534958306e-2,5.169091033818206e-2,4.26450085290017e-2,5.1543910308782054e-2,1.980300396060079e-2,4.075200815040162e-2,6.182991236598247e-2,7.315191463038292e-2,2.2944004588800915e-2,5.3760910752182145e-2,9.190391838078367e-2,5.867591173518234e-2,2.382600476520095e-2,4.01260080252016e-2,5.090091018018203e-2,6.876491375298274e-2,5.856491171298234e-2,1.4261002852200569e-2,3.210200642040128e-2,6.600491320098263e-2}; //renormalised
const double dDAYVal[190] = {27,98,32,120,0,905,36,23,0,0,89,246,103,134,0,198,1,148,1153,0,716,240,9,139,125,11,28,81,23,240,535,86,28,606,43,10,65,64,77,24,44,18,61,0,7,41,15,34,0,0,73,11,7,44,257,26,464,318,71,0,153,83,27,26,46,18,72,90,1,0,0,114,30,17,0,336,527,243,18,14,14,0,0,0,0,15,48,196,157,0,92,250,103,42,13,19,153,51,34,94,12,32,33,17,11,409,154,495,95,161,56,79,234,35,24,17,96,62,46,245,371,26,229,66,16,53,34,30,22,192,33,136,104,13,78,550,0,201,23,0,0,0,0,0,27,0,46,0,0,76,0,75,0,24,8,95,0,96,0,22,0,127,37,28,13,0,698,0,34,42,61,208,24,15,18,49,35,37,54,44,889,175,10,258,12,48,30,157,0,28};
const double dDAYFreq[20] = {0.087127,0.040904,0.040432,0.046872,0.033474,0.038255,0.049530,0.088612,0.033618,0.036886,0.085357,0.080482,0.014753,0.039772,0.050680,0.069577,0.058542,0.010494,0.029916,0.064718};
const double dmtREVVal[190] = {23.18,26.95,13.24,17.67,1.90,794.38,59.93,103.33,58.94,1.90,1.90,220.99,173.56,55.28,75.24,9.77,1.90,63.05,583.55,1.90,313.56,120.71,23.03,53.30,56.77,30.71,6.75,28.28,13.90,165.23,496.13,113.99,141.49,582.40,49.12,1.90,96.49,1.90,27.10,4.34,62.73,8.34,3.31,5.98,12.26,25.46,15.58,15.16,1.90,25.65,39.70,1.90,2.41,11.49,329.09,8.36,141.40,608.70,2.31,1.90,465.58,313.86,22.73,127.67,19.57,14.88,141.88,1.90,65.41,1.90,6.18,47.37,1.90,1.90,11.97,517.98,537.53,91.37,6.37,4.69,15.20,4.98,70.80,19.11,2.67,1.90,48.16,84.67,216.06,6.44,90.82,54.31,23.64,73.31,13.43,31.26,137.29,12.83,1.90,60.97,20.63,40.10,50.10,18.84,17.31,387.86,6.04,494.39,69.02,277.05,54.11,54.71,125.93,77.46,47.70,73.61,105.79,111.16,64.29,169.90,480.72,2.08,238.46,28.01,179.97,94.93,14.82,11.17,44.78,368.43,126.40,136.33,528.17,33.85,128.22,597.21,1.90,21.95,10.68,19.86,33.60,1.90,1.90,10.92,7.08,1.90,32.44,24.00,21.71,7.84,4.21,38.58,9.99,6.48,1.90,191.36,21.21,254.77,38.82,13.12,3.21,670.14,25.01,44.15,51.17,39.96,465.58,16.21,64.92,38.73,26.25,195.06,7.64,1.90,1.90,1.90,19.00,21.14,2.53,1.90,1222.94,91.67,1.90,387.54,6.35,8.23,1.90,204.54,5.37,1.90};
const double dmtREVFreq[20] = {0.072,0.019,0.039,0.019,0.006,0.025,0.024,0.056,0.028,0.088,0.169,0.023,0.054,0.061,0.054,0.072,0.086,0.029,0.033,0.043};
const double dLGVal[190] = {0.425093,0.276818,0.751878,0.395144,0.123954,5.076149,2.489084,0.534551,0.528768,0.062556,0.969894,2.807908,1.695752,0.523386,0.084808,1.038545,0.363970,0.541712,5.243870,0.003499,4.128591,2.066040,0.390192,1.437645,0.844926,0.569265,0.267959,0.348847,0.358858,2.426601,4.509238,0.927114,0.640543,4.813505,0.423881,0.311484,0.149830,0.126991,0.191503,0.010690,0.320627,0.072854,0.044265,0.008705,0.108882,0.395337,0.301848,0.068427,0.015076,0.594007,0.582457,0.069673,0.044261,0.366317,4.145067,0.536518,6.326067,2.145078,0.282959,0.013266,3.234294,1.807177,0.296636,0.697264,0.159069,0.137500,1.124035,0.484133,0.371004,0.025548,0.893680,1.672569,0.173735,0.139538,0.442472,4.273607,6.312358,0.656604,0.253701,0.052722,0.089525,0.017416,1.105251,0.035855,0.018811,0.089586,0.682139,1.112727,2.592692,0.023918,1.798853,1.177651,0.332533,0.161787,0.394456,0.075382,0.624294,0.419409,0.196961,0.508851,0.078281,0.249060,0.390322,0.099849,0.094464,4.727182,0.858151,4.008358,1.240275,2.784478,1.223828,0.611973,1.739990,0.990012,0.064105,0.182287,0.748683,0.346960,0.361819,1.338132,2.139501,0.578987,2.000679,0.425860,1.143480,1.080136,0.604545,0.129836,0.584262,1.033739,0.302936,1.136863,2.020366,0.165001,0.571468,6.472279,0.180717,0.593607,0.045376,0.029890,0.670128,0.236199,0.077852,0.268491,0.597054,0.111660,0.619632,0.049906,0.696175,2.457121,0.095131,0.248862,0.140825,0.218959,0.314440,0.612025,0.135107,1.165532,0.257336,0.120037,0.054679,5.306834,0.232523,0.299648,0.131932,0.481306,7.803902,0.089613,0.400547,0.245841,3.151815,2.547870,0.170887,0.083688,0.037967,1.959291,0.210332,0.245034,0.076701,0.119013,10.649107,1.702745,0.185202,1.898718,0.654683,0.296501,0.098369,2.188158,0.189510,0.249313};
const double dLGFreq[20] = {7.906592093407908e-2,5.594094405905594e-2,4.197695802304198e-2,5.3051946948053055e-2,1.2936987063012939e-2,4.0766959233040766e-2,7.158592841407159e-2,5.733694266305734e-2,2.2354977645022357e-2,6.215693784306216e-2,9.908090091909909e-2,6.45999354000646e-2,2.2950977049022953e-2,4.2301957698042306e-2,4.403995596004404e-2,6.11969388030612e-2,5.328694671305329e-2,1.2065987934012068e-2,3.4154965845034156e-2,6.914693085306915e-2}; // renormalised
const double dmtMAMVal[190] = {32,2,4,11,0,864,0,186,0,0,0,246,8,49,0,0,0,0,569,0,274,78,18,47,79,0,0,22,8,232,458,11,305,550,22,0,75,0,19,0,41,0,0,0,0,21,6,0,0,27,20,0,0,26,232,0,50,408,0,0,242,215,0,0,6,4,76,0,21,0,0,22,0,0,0,378,609,59,0,0,6,5,7,0,0,0,0,57,246,0,11,53,9,33,2,0,51,0,0,53,5,43,18,0,17,342,3,446,16,347,30,21,112,20,0,74,65,47,90,202,681,0,110,0,114,0,4,0,1,360,34,50,691,8,78,614,5,16,6,0,65,0,0,0,0,0,12,0,13,0,7,17,0,0,0,156,0,530,54,0,1,1525,16,25,67,0,682,8,107,0,14,398,0,0,10,0,33,20,5,0,2220,100,0,832,6,0,0,237,0,0};
const double dmtMAMFreq[20] = {0.0692,0.0184,0.0400,0.0186,0.0065,0.0238,0.0236,0.0557,0.0277,0.0905,0.1675,0.0221,0.0561,0.0611,0.0536,0.0725,0.0870,0.0293,0.0340,0.0428};
const double drtREVVal[190] = {34,51,35,10,30,384,439,92,128,1,32,221,236,78,70,81,10,79,542,1,372,135,41,94,61,48,18,70,30,90,320,91,124,387,34,68,1,24,35,1,104,33,1,1,34,45,18,15,5,110,54,21,3,51,385,38,593,123,20,16,309,141,30,76,34,23,235,57,1,1,156,158,1,37,116,375,581,134,1,7,49,1,70,1,1,7,141,64,179,14,247,97,24,33,55,1,68,52,17,44,10,22,43,1,11,460,102,294,136,75,225,95,152,183,4,24,77,1,20,134,258,64,148,55,117,146,82,7,49,72,25,110,131,69,62,671,5,13,16,1,55,10,17,23,48,39,47,6,111,182,9,14,1,55,47,28,1,131,45,1,21,307,26,64,1,74,1017,14,31,34,176,197,29,21,6,295,36,35,3,1,1048,112,19,236,92,25,39,196,26,59};
const double drtREVFreq[20] = {0.0646,0.0453,0.0376,0.0422,0.0114,0.0606,0.0607,0.0639,0.0273,0.0679,0.1018,0.0751,0.0150,0.0287,0.0681,0.0488,0.0622,0.0251,0.0318,0.0619};
const double dmtArtFreq[20] = {0.054116,0.018227,0.039903,0.020160,0.009709,0.018781,0.024289,0.068183,0.024518,0.092638,0.148658,0.021718,0.061453,0.088668,0.041826,0.091030,0.049194,0.029786,0.039443,0.057700};
const double dmtArtVal[190] = {0.2,0.2,0.2,1,4,500,254,36,98,11,0.2,154,262,0.2,0.2,0.2,0.2,183,862,0.2,262,200,0.2,121,12,81,3,44,0.2,41,180,0.2,12,314,15,0.2,26,2,21,7,63,11,7,3,0.2,4,2,13,1,79,16,2,1,6,515,0.2,209,467,2,0.2,349,106,0.2,0.2,3,4,121,5,79,0.2,312,67,0.2,56,0.2,515,885,106,13,5,20,0.2,184,0.2,0.2,1,14,118,263,11,322,49,0.2,17,0.2,0.2,39,8,0.2,1,0.2,12,17,5,15,673,3,398,44,664,52,31,226,11,7,8,144,112,36,87,244,0.2,166,0.2,183,44,43,0.2,19,204,48,70,289,14,47,660,0.2,0.2,8,0.2,22,7,11,2,0.2,0.2,21,16,71,54,0.2,2,0.2,1,4,251,0.2,72,87,8,9,191,12,20,117,71,792,18,30,46,38,340,0.2,23,0.2,350,0.2,14,3,0.2,1855,85,26,281,52,32,61,544,0.2,2};
const double dcpRevFreq[20] = {0.0755,0.0621,0.0410,0.0371,0.0091,0.0382,0.0495,0.0838,0.0246,0.0806,0.1011,0.0504,0.0220,0.0506,0.0431,0.0622,0.0543,0.0181,0.0307,0.0660};
const double dcpRevVal[190] = {105,227,357,175,43,4435,669,823,538,10,157,1745,768,400,10,499,152,1055,3691,10,3122,665,243,653,431,303,133,379,66,715,1405,331,441,1269,162,19,145,136,168,10,280,92,148,40,29,197,203,113,10,396,286,82,20,66,1745,236,4482,2430,412,48,3313,2629,263,305,345,218,185,125,61,47,159,202,113,21,10,1772,1351,193,68,53,97,22,726,10,145,25,127,454,1268,72,327,490,87,173,170,285,323,185,28,152,117,219,302,100,43,2440,385,2085,590,2331,396,568,691,303,216,516,868,93,487,1202,1340,314,1393,266,576,241,369,92,32,1040,156,918,645,148,260,2151,14,230,40,18,435,53,63,82,69,42,159,10,86,468,49,73,29,56,323,754,281,1466,391,142,10,1971,89,189,247,215,2370,97,522,71,346,968,92,83,75,592,54,200,91,25,4797,865,249,475,317,122,167,760,10,119};
const double dVTVal[190] = {1.2412691067876198,1.2184237953498958,1.5720770753326880,1.3759368509441177,0.7550654439001206,7.8584219153689405,2.4731223087544874,1.4414262567428417,0.9784679122774127,0.2272488448121475,2.2155167805137470,5.5120819705248678,3.0143201670924822,1.6562495638176040,0.4587469126746136,2.3379911207495061,1.3542404860613146,2.0093434778398112,9.6883451875685065,0.4519167943192672,6.8124601839937675,3.3386555146457697,1.3121700301622004,2.4117632898861809,1.9142079025990228,1.1034605684472507,0.8776110594765502,1.3860121390169038,0.9615841926910841,4.9238668283945266,6.1974384977884114,2.1459640610133781,1.5196756759380692,7.9943228564946525,1.6360079688522375,0.8561248973045037,0.8908203061925510,0.4323005487925516,0.9179291175331520,0.2161660372725585,0.9126668032539315,0.4882733432879921,0.4035497929633328,0.2888075033037488,0.5787937115407940,1.0778497408764076,0.8386701149158265,0.4098311270816011,0.3574207468998517,1.4081315998413697,1.3318097154194044,0.5610717242294755,0.3578662395745526,1.0765007949562073,6.0019110258426362,1.4932055816372476,10.017330817366002,4.4034547578962568,1.4521790561663968,0.3371091785647479,6.0519085243118811,4.3290086529582830,0.8945563662345198,1.8085136096039203,0.6244297525127139,0.5642322882556321,1.9006455961717605,1.2488638689609959,0.9378803706165143,0.4075239926000898,1.2213054800811556,1.9106190827629084,0.7471936218068498,0.5954812791740037,1.3808291710019667,6.7597899772045418,8.0327792947421148,1.7129670976916258,0.6883439026872615,0.4224945197276290,0.5044944273324311,0.1675129724559251,1.6953951980808002,0.3573432522499545,0.2317194387691585,0.3693722640980460,1.3629765501081097,2.2864286949316077,4.3611548063555778,0.3910559903834828,2.3201373546296349,2.7355620089953550,1.3091837782420783,0.7103720531974738,1.0714605979577547,0.4326227078645523,2.3019177728300728,1.5132807416252063,0.7744933618134962,1.8370555852070649,0.4811402387911145,1.0084320519837335,1.3918935593582853,0.4953193808676289,0.3746821107962129,6.4208961859142883,1.9202994262316166,6.1234512396801764,2.2161944596741829,3.6366815408744255,2.3193703643237220,1.8273535587773553,3.0637776193717610,1.9699895187387506,0.6047491507504744,0.8953754669269811,1.9776630140912268,1.0657482318076852,1.1079144700606407,3.5465914843628927,5.2892514169776437,1.3363401740560601,3.8852506105922231,1.5066839872944762,1.7557065205837685,2.1576510103471440,1.5839981708584689,0.7147489676267383,1.6136654573285647,2.6344778384442731,1.0192004372506540,2.5513781312660280,3.3628488360462363,0.6882725908872254,1.9485376673137556,8.8479984061248178,0.5488578478106930,1.5170142153962840,0.1808525752605976,0.2496584188151770,1.6275179891253113,0.8959082681546182,0.4198391148111098,0.9349753595598769,0.6301954684360302,0.5604648274060783,1.5183114434679339,0.5851920879490173,1.4680478689711018,3.3448437239772266,0.4326058001438786,0.6791126595939816,0.4514203099376473,0.5411769916657778,0.8912614404565405,1.0894926581511342,0.7447620891784513,2.1579775140421025,0.9183596801412757,0.5818111331782764,0.3374467649724478,7.7587442309146040,0.8626796044156272,1.2452243224541324,0.7835447533710449,1.0899165770956820,10.384852333133459,0.4819109019647465,0.9547229305958682,0.8564314184691215,4.5377235790405388,4.6501894691803214,0.7807017855806767,0.4586061981719967,0.4594535241660911,2.2627456996290891,0.6366932501396869,0.8940572875547330,0.6193321034173915,0.5333220944030346,14.872933461519061,3.5458093276667237,0.7801080335991272,4.0584577156753401,1.7039730522675411,0.5985498912985666,0.9305232113028208,3.4242218450865543,0.5658969249032649,1.0000000000000000};
const double dVTFreq[20] = {0.0770764620135024,0.0500819370772208,0.0462377395993731,0.0537929860758246,0.0144533387583345,0.0408923608974345,0.0633579339160905,0.0655672355884439,0.0218802687005936,0.0591969699027449,0.0976461276528445,0.0592079410822730,0.0220695876653368,0.0413508521834260,0.0476871596856874,0.0707295165111524,0.0567759161524817,0.0127019797647213,0.0323746050281867,0.0669190817443274};
const double dBLOSUM62Val[190]={0.73579038969751,0.48539105546575,1.29744670513370,0.54316182089867,0.50096440855513,3.18010004821610,1.45999531047000,0.22782657420895,0.39735894989702,0.24083661480204,1.19970570460200,3.02083361006360,1.83921614699200,1.19094570339600,0.32980150463028,1.17094904279990,1.36057419042030,1.24048850863960,3.76162520836850,0.14074889181440,5.52891917792820,1.95588357495950,0.41876330851753,1.35587234448450,0.79847324896839,0.41820319228376,0.60984630538281,0.42357999217628,0.71624144499779,1.45614116633600,2.41450143420810,0.77814266402188,0.35405810983129,2.43534113114010,1.62689105698170,0.53985912495418,0.60589900368677,0.23203644514174,0.28301732627800,0.41855573246161,0.77489402279418,0.23620245120365,0.18684804693170,0.18929629237636,0.25271844788492,0.80001653051838,0.62271166969249,0.21188815961519,0.21813157759360,0.83184264014158,0.58073709318144,0.37262517508685,0.21772115923623,0.34807220979697,3.89096377330350,1.29520126678330,5.41111514148890,1.59313704345740,1.03244792495210,0.28507880090648,3.94527767451460,2.80242715167870,0.75204244030271,1.02250703588900,0.40619358664202,0.44557027426059,1.25375826666350,0.98369298745695,0.64844127878707,0.22262189795786,0.76768882347954,2.49489607711270,0.55541539747043,0.45943617357855,0.98431152535870,3.36479776310420,6.03055937957160,1.07306118433190,0.49296467974759,0.37164469320875,0.35486124922252,0.28173069420651,0.44133747118660,0.14435695975031,0.29140908416530,0.36816646445253,0.71453370392764,1.51735932595390,2.06483970323750,0.26692475051102,1.77385516883050,1.17327590092390,0.44813366171831,0.49488704370192,0.73062827299842,0.35600849876863,0.85857057567418,0.92656393484598,0.50408659952683,0.52700733915060,0.38835540920564,0.37455568747097,1.04738345072150,0.45412362510273,0.23359790962888,4.32509268705660,1.12278310420960,2.90410165645600,1.58275414206530,1.19718841509420,1.93487092459650,1.76989323893730,1.50932625322360,1.11702976291050,0.35754441245967,0.35296918452729,1.75216591781950,0.91872341574605,0.54002764482413,1.16912957771570,1.72917801948500,0.91466595456337,1.89817363453320,0.93418750943056,1.11983135851600,1.27748029459560,1.07109723600730,0.64143601140497,0.58540709022472,1.17909119726010,0.91525985769421,1.30387520079870,1.48854805372180,0.48820611879305,1.00545168314880,5.15155629227040,0.46583936772479,0.42638231012175,0.19148204624678,0.14534504627853,0.52766441887169,0.75865380864172,0.40763564893830,0.50835892463812,0.30124860078016,0.34198578754023,0.69147463459998,0.33224304063396,0.88810109815193,2.07432489349650,0.25221483002727,0.38792562209837,0.51312812689059,0.71820669758623,0.72051744121611,0.53822251903674,0.26142220896504,0.47023773369610,0.95898974285014,0.59671930034577,0.30805573703500,4.21895396938900,0.67461709322842,0.81124585632307,0.71799348690032,0.95168216224591,6.74726043080080,0.36940531935451,0.79675152076106,0.80101024319939,4.05441900655800,2.18777452200450,0.43838834377202,0.31285879799342,0.25812928941763,1.11635247860620,0.53078579012486,0.52425384633796,0.25334079019018,0.20155597175031,8.31183940545820,2.23140568891310,0.49813847530407,2.57585075531530,0.83811961017754,0.49690841067567,0.56192545744165,2.25307405117630,0.26650873142646,1.00000000000000};
const double dBLOSUM62Freq[20]={0.074,0.052,0.045,0.054,0.025,0.034,0.054,0.074,0.026,0.068,0.099,0.058,0.025,0.047,0.039,0.057,0.051,0.013,0.032,0.073};
const double dHIVbVal[190]={0.30750700,0.00500000,0.29554300,1.45504000,0.00500000,17.66120000,0.12375800,0.35172100,0.08606420,0.00500000,0.05511280,3.42150000,0.67205200,0.00500000,0.00500000,1.48135000,0.07492180,0.07926330,10.58720000,0.00500000,2.56020000,2.13536000,3.65345000,0.32340100,2.83806000,0.89787100,0.06191370,3.92775000,0.08476130,9.04044000,7.64585000,1.91690000,0.24007300,7.05545000,0.11974000,0.00500000,0.00500000,0.67728900,0.68056500,0.01767920,0.00500000,0.00500000,0.00609079,0.00500000,0.10311100,0.21525600,0.70142700,0.00500000,0.00876048,0.12977700,1.49456000,0.00500000,0.00500000,1.74171000,5.95879000,0.00500000,20.45000000,7.90443000,0.00500000,0.00500000,6.54737000,4.61482000,0.52170500,0.00500000,0.32231900,0.08149950,0.01866430,2.51394000,0.00500000,0.00500000,0.00500000,0.30367600,0.17578900,0.00500000,0.00500000,11.20650000,5.31961000,1.28246000,0.01412690,0.00500000,0.00500000,0.00500000,9.29815000,0.00500000,0.00500000,0.29156100,0.14555800,3.39836000,8.52484000,0.03426580,0.18802500,2.12217000,1.28355000,0.00739578,0.03426580,0.00500000,4.47211000,0.01202260,0.00500000,2.45318000,0.04105930,2.07757000,0.03138620,0.00500000,0.00500000,2.46633000,3.47910000,13.14470000,0.52823000,4.69314000,0.11631100,0.00500000,4.38041000,0.38274700,1.21803000,0.92765600,0.50411100,0.00500000,0.95647200,5.37762000,15.91830000,2.86868000,6.88667000,0.27472400,0.73996900,0.24358900,0.28977400,0.36961500,0.71159400,8.61217000,0.04376730,4.67142000,4.94026000,0.01412690,2.01417000,8.93107000,0.00500000,0.99133800,0.00500000,0.00500000,2.63277000,0.02665600,0.00500000,1.21674000,0.06951790,0.00500000,0.74884300,0.00500000,0.08907800,0.82934300,0.04445060,0.02487280,0.00500000,0.00500000,0.00991826,1.76417000,0.67465300,7.57932000,0.11303300,0.07926330,0.00500000,18.69430000,0.14816800,0.11198600,0.00500000,0.00500000,15.34000000,0.03043810,0.64802400,0.10565200,1.28022000,7.61428000,0.08124540,0.02665600,1.04793000,0.42002700,0.02091530,1.02847000,0.95315500,0.00500000,17.73890000,1.41036000,0.26582900,6.85320000,0.72327400,0.00500000,0.07492180,0.70922600,0.00500000,0.04105930};
const double dHIVbFreq[20]={0.060490222,0.066039665,0.044127815,0.042109048,0.020075899,0.053606488,0.071567447,0.072308239,0.022293943,0.069730629,0.098851122,0.056968211,0.019768318,0.028809447,0.046025282,0.050604330,0.053636813,0.033011601,0.028350243,0.061625237};
const double dHIVwVal[190]={0.0744808,0.6175090,0.1602400,4.4352100,0.0674539,29.4087000,0.1676530,2.8636400,0.0604932,0.0050000,0.0050000,10.6746000,0.3420680,0.0050000,0.0050000,5.5632500,0.0251632,0.2015260,12.1233000,0.0050000,3.2065600,1.8685000,13.4379000,0.0604932,10.3969000,0.0489798,0.0604932,14.7801000,0.0050000,6.8440500,8.5987600,2.3177900,0.0050000,18.5465000,0.0050000,0.0050000,0.0050000,1.3406900,0.9870280,0.1451240,0.0050000,0.0342252,0.0390512,0.0050000,0.0050000,0.1602400,0.5867570,0.0050000,0.0050000,0.0050000,2.8904800,0.1298390,0.0489798,1.7638200,9.1024600,0.5927840,39.8897000,10.6655000,0.8943130,0.0050000,13.0705000,23.9626000,0.2794250,0.2240600,0.8174810,0.0050000,0.0050000,3.2865200,0.2015260,0.0050000,0.0050000,0.0050000,0.0050000,0.0489798,0.0050000,17.3064000,11.3839000,4.0956400,0.5979230,0.0050000,0.0050000,0.0050000,0.3629590,0.0050000,0.0050000,0.0050000,0.0050000,1.4828800,7.4878100,0.0050000,0.0050000,1.0098100,0.4047230,0.3448480,0.0050000,0.0050000,3.0450200,0.0050000,0.0050000,13.9444000,0.0050000,9.8309500,0.1119280,0.0050000,0.0342252,8.5942000,8.3502400,14.5699000,0.4278810,1.1219500,0.1602400,0.0050000,6.2796600,0.7251570,0.7400910,6.1439600,0.0050000,0.3925750,4.2793900,14.2490000,24.1422000,0.9282030,4.5420600,0.6303950,0.0050000,0.2030910,0.4587430,0.0489798,0.9595600,9.3634500,0.0050000,4.0480200,7.4131300,0.1145120,4.3370100,6.3407900,0.0050000,5.9656400,0.0050000,0.0050000,5.4989400,0.0443298,0.0050000,2.8258000,0.0050000,0.0050000,1.3703100,0.0050000,0.0050000,0.0050000,0.0050000,1.1015600,0.0050000,0.0050000,0.0050000,5.0647500,2.2815400,8.3483500,0.0050000,0.0050000,0.0050000,47.4889000,0.1145120,0.0050000,0.0050000,0.5791980,4.1272800,0.0050000,0.9331420,0.4906080,0.0050000,24.8094000,0.2794250,0.0744808,2.9178600,0.0050000,0.0050000,2.1995200,2.7962200,0.8274790,24.8231000,2.9534400,0.1280650,14.7683000,2.2800000,0.0050000,0.8626370,0.0050000,0.0050000,1.3548200};
const double dHIVwFreq[20]={0.0377494,0.0573210,0.0891129,0.0342034,0.0240105,0.0437824,0.0618606,0.0838496,0.0156076,0.0983641,0.0577867,0.0641682,0.0158419,0.0422741,0.0458601,0.0550846,0.0813774,0.0195970,0.0205847,0.0515639};


////////////////////////////////////////////////////////
// Data compression routine
vector <int> GetCompressedData(CTree *T,CData *D);	// Returns vector of size (m_iNoIntNode * m_iSize) describing the compression (see implementation for more details)
void ColSorter(int NT, int NF, CTree *T, CData *D, vector < vector <int> > &ColInfo,bool First = true);
vector <int> DoCompress(vector <vector <int> > &CV, CTree *T, CData *D);

///////////////////////////////////////////////////////////////////////////
// Classes expanding on CPar for Q matrix specific parameters
///////////////////////////////////////////////////////////////////////////

// Standard Q matrix parameters
class CQPar : public CPar	{
public:
	CQPar(string Name, int Char, double Value, bool Optimise=true, double Lowbound = 0.0, double Upbound=MAX_PAR_VALUE,ParOp Oper=MULTIPLY);	// Core parameter construction routine
	// Variables mapping to the Q matrix
	int m_iChar;
	vector <int> m_viQMap;									// Elements of the Q matrix it applies to
	vector <int> m_viApply2Q;								// The list of Q matrix IDs that the parameter applies to (empty = all)
	// Functions
	void AddQij(int i, int j, bool Rev = true);				// Adds coordinates of Q to
	void Par2Q(double *QMat, int QMatID);					// Applies the parameter to the Q matrix of number QNum
	// Updated output function
	ostream &Output(ostream &os);
};

// The alpha parameter of the gamma distribution
#define INITIAL_GAMMA MODEL_INITIAL_GAMMA
#define MAX_ALFA 100
class CGammaPar : public CPar	{
public:
	CGammaPar(string Name, int Char, CPar* Rate,double Value = INITIAL_GAMMA, bool Optimise=true, double Lowbound = 0.01, double Upbound=MAX_ALFA,ParOp Oper=REPLACE);	// Constructor
	~CGammaPar();							// Destructor
	// public functions
	void GlobalApply();						// How the parameter is applied
	void AddRateToGamma(CBaseProcess *Proc);// Control this rate using the gamma distribution
	int NoCat() { return m_iNoCat; }		// Returns the number of rate categories
	ostream &Output(ostream &os);			// Output function
private:
	// Private variables
	int m_iNoCat;						// Number of categories in the gamma distribution
	vector <CQPar *> m_arpRates;		// Rate parameters
	CPar *m_pProcRate;					// Pointer to the rate of the process
	CQPar *m_pProb;						// Probability of the base process
};

// This is the Q matrix class
/////////////////////////////////////////////////////////////////////////
// This contains the space for a single Q matrix and the functions
// for turning it into a P(t) matrix
class CQMat {
public:
	CQMat(EDataType Type,string Name = "Unnamed Q");
	CQMat(int Char, EDataType Type = OTHER,string Name = "Unnamed Q");
	// Access functions
	double *Q(int i=0,int j=0)	{ return &m_ardQMat[(i*m_iChar)+j]; };		// Access to the Q matrix
	int Char()					{ return m_iChar; };						// Access to the number of characters
	vector <double> Eqm();													// Access to the processes equilibrium
	EDataType Type()			{ return m_DataType; };						// Access to the datatype
	int ID()	{ return m_iQMatID; }										// Returns the identifier of the matrix
	void OutQRoot(ostream &os=cout) { int i; FOR(i,m_iChar) { os << m_ardRoot[i] << " "; } };	// Outputs the eigen values
	void OutQ(ostream &os=cout)		{ int i; FOR(i,m_iChar2) { if(i%m_iChar==0) { os << endl; } os << m_ardQMat[i] << " "; } } // Outputs the Q matrix
	void OutU(ostream &os=cout)		{ int i; FOR(i,m_iChar2) { if(i%m_iChar==0) { os << endl; } os << m_ardU[i] << " "; } } // Outputs the U matrix
	void OutV(ostream &os=cout)		{ int i; FOR(i,m_iChar2) { if(i%m_iChar==0) { os << endl; } os << m_ardV[i] << " "; } } // Outputs the V matrix
	void OutRoot(ostream &os=cout)	{ int i; FOR(i,m_iChar)  { os << m_ardRoot[i] << " "; } }
	bool IsLocked() { return FlipBool(m_bAllowModelUpdate); }							// Returns whether the process is locked
	void Lock()		{ m_bAllowModelUpdate = false; }						// Locks the model so the Q matrix and the eigen vectors can't be changed
	void Unlock()	{ m_bAllowModelUpdate = true; }							// Unlock the model so Q and eigens can change
	// Linear algebra functions
	void Decompose(double *Eqm,bool Scale = true, bool REV = true,double Rate = 1.0);			// Decompose the matrix to its eigen vectors and values (Eqm contains equilibrium distribution if known)
	void Decompose(vector <double> Eqm,bool Scale = true,bool REV = true,double Rate = 1.0);	// Decompose the matrix to its eigen vectors and values (Eqm contains equilibrium distribution if known)
	bool MakePT(double T,double *PT);										// Returns the P(t) matrix
	virtual double OverallSubRate(bool ForceWholeProcess = false);			// Returns the overall substitution rate of the process
	virtual double OverallTransRate() { return 0.0; }						// Returns the overall transition rate between hidden states
	virtual double TransRateFrom(int From) { return 0.0; };					// Returns the overall transition rate from a hidden state
	virtual double TransRateTo(int To) { return 0.0; };						// Returns the overall transition rate to a hidden state
	virtual void ScaleQ(double Rate);												// Scale the Q matrix to the overall rate Rate

	// Parameter application functions
	void InitQ(double Val);								// Initialise Q matrix to a particular value
	void ApplyPar2Q(CQPar *Par);						// Apply the parameter to the Q matrix
	void DoQDiag();															// Make all rows sum to 0
	friend ostream &operator<<(ostream &os,CQMat &Mat) { return Mat.Output(os); };

protected:
	// Variables
	bool m_bScaleReady;					// Whether the process is ready to be scaled (requires the eqm distribution to be verified)
	double *m_ardEqm;					// The equilibrium of the process
	double m_dScale;					// Current overall rate of process
	bool m_bAlwaysI;					// When the rate is set to zero, P(t) will always produce an identity matrix
	// Variables
	int m_iQMatID;						// Unique identifier for the Q matrix
	string m_sName;						// Name for the process
	EDataType m_DataType;				// The datatype stored
	int m_iChar;						// The number of characters in the data
	int m_iChar2;						// The number of characters in the data squared
	bool m_bAllowModelUpdate;			// Whether the model can be updated
	double *m_ardQMat;					// The Q matrix
	double *m_ardU,*m_ardV,*m_ardRoot;	// The eigen information
	double *m_ardRootEqm;				// The square root of the eqm distribution (for fast eigen decomposition);
	// Memory management functions
	int GetQMatID()	{ return QMatID++; }									// Get the identifier for a Q matrix
	void MakeSpace(int Char,EDataType Type,string Name);
	void MakeSpace(EDataType Type,string Name)	{ MakeSpace(NumStates(Type),Type,Name); }
	void CleanSpace();
	virtual ostream &Output(ostream &os);									// Output function
};

///////////////////////////////////////////////////////////////////////////
// This is the Q matrix class for covarion models
// --
// The only difference between this and the class CQMat is how scaling
//  is performed
///////////////////////////////////////////////////////////////////////////
class CTHMMQMat : public CQMat	{
public:
	// Constructor
	CTHMMQMat(int DataChar, int m_iChar, EDataType Type,string Name = "Unnamed Q");
	// Destructor
	~CTHMMQMat() { /* empty */ };
	// Interaction functions
	double OverallSubRate(bool FWP);		// Returns the overall rate of the *SUBSTITUTION* process -- ignores transitions between hidden states
	double OverallTransRate();				// Returns the overall rate of transitions in the *HIDDEN* process -- ignores substitions between observable states
	double TransRateFrom(int From);			// Returns the overall transition rate from a state
	double TransRateTo(int To);				// Returns the overall transition rate to a state
	void ScaleQ(double Rate);				// Scales only the substitution process in Q
private:
	int m_iDataChar;							// Size of the data (used for calculating
};

///////////////////////////////////////////////////////////////////////////
// Eqm defining classes - Used for applying a defined Eqm to a CQMat
///////////////////////////////////////////////////////////////////////////
// equ			= equiprobable (e.g. Jukes and Cantor)
// obs			= single frequencies produced from data
// complex		= seperate frequency vectors for each hidden process
// complex_GC	= single frequency vector with different GC content for each process

class CBaseEqm	{
public:
	// Core functions
	CBaseEqm(int Char, vector <CQPar *> *ProcPar);			// The constructor function
	~CBaseEqm();											// The destructor function
	// Interaction functions
	void ApplyEqm2QMat(double *Q, int MatID);				// Applies the equilibrium distribution to these matrices (based on m_viQMatID)
	bool IsID(int ID);										// Whether the equilibrium distribution applies to a particular matrix
	void AddMatID(int ID);									// Add matrix ID that equilibrium applies to
	void SetOpt(bool Optimise);								// Whether to optimise the Eqm distribution
	int Char() { return m_iChar; }							// Returns the number of characters in teh eqm distribution
	virtual vector <double> Eqm() {							// Returns the equilibrium distribution for the whole process
		vector <double> v; Error("Empty eqm distribution"); return v; }
	// HMM based eqm functions
	virtual vector <double> SubEqm(int HiddenState) {		// Returns the equilibrium distribution for a substitution state
		vector <double> v; Error("Empty CBaseEqm::SubEqm() function"); return v; }
	virtual vector <double> TransEqm() {					// Returns the equilibrium for the transition matrix
		vector <double> v; Error("Empty CBaseEqm::TransEqm() function"); return v; }
	vector <double *> OptimiserValues();					// Returns a vector containing the optimiser values for the eqm distribution
	void SetDoBasicNonReversibleCovarion(bool Val) { DoBasicNonReversibleCovarion = Val; }
	virtual void Shuffle()	{};										// Blank shuffle routine for when there are multiple eqms
	virtual void ResetEqm(vector <double> Eqm, bool RandomFactor) { Error("Empty CBaseEqm::ResetEqm function...\n"); };// Resets the eqm parameters to the values in the vector
protected:
	// Variables
	int m_iChar;											// Number of character states in the equilibrium
	vector <int> m_viQMatID;								// The CQMat objects that the eqm process is applied to
	vector <CPar *> m_vpEqmPar;								// Pointers to the parameters that define the equilibrium distribution
	vector <CQPar *> *m_pProcPar;							// Pointer to process parameters
	double *m_ardQModifier;									// The modifier matrix (element-wise)
	// Elements redundant for some kinds of models
	bool DoBasicNonReversibleCovarion;					// Produce a model that relies on \pi^{k}_i == \pi^{l}_i i.e. standard covarion model
	// Checks as to whether a recalculation of m_ardQModifier is required
	vector <double> m_vdOldParVals;							// Old parameters values (used for whether updating is required)
	bool CheckQ();											// Checks whether the equilibrium distribution has changed
	// The function that actually applies the distribution
	virtual void AdapterFunc() { Error("Empty adaptor function"); }	// Creates m_ardQModifier matrix for the equilibrium distribution

};

/* ********************** Class describing a simple eqm distribution *********************** */
class CSimpleEqm : public CBaseEqm {
public:
	// Core functions
	CSimpleEqm(int Char, vector <CQPar *> *ProcPar, CData *Data);			// Constructor
	CSimpleEqm(int Char, vector <CQPar *> *ProcPar, vector <double> Eqm);	// Constructor
	~CSimpleEqm();															// Destructor
	// Interaction
	void ResetEqm(vector <double> Eqm, bool RandomFactor);					// Resets the equlibrium distribution
	vector <double> Eqm();									// Returns equilibrium distribution
private:
	void ConstructSimpleEqm(vector <double> eqm);			// Construct the equilibrium
	void AdapterFunc();										// Applies the simple eqm distribution to Q
};

/* ************************ Class describing codon equ distribution ************************* */
class CCodonEqm : public CBaseEqm {
public:
	// Core functions
	CCodonEqm(ECodonEqm CE, int GenCode, vector <CQPar *> *ProcPar,CData *D);			// Constructor
	~CCodonEqm();																// Destructor
	// Interaction
	ECodonEqm m_EqmType;														// Stores the eqm type
	int m_iGenCode;																// Store the genetic code
	vector <double> Eqm();														// Returns eqm distribution
private:
	ECodonEqm m_CE;
	void ConstructCodonEqm(vector <double> eqm);								// Construct the eqm
	void AdapterFunc();															// Applies the eqm distribution to Q
	CData *m_pData;																// Pointer to data
};

/* ************************ Class describing covarion eqm process *************************** */
class CCovEqm : public CBaseEqm {
public:
	// Core functions
	// Two constructors to create simple covarion equilibriums (1 observable distribution for all hidden states)
	CCovEqm(int Char, vector <CQPar *> *ProcPar, CData *Data, vector <CPar *> *StateProbs);			// Constructor
	CCovEqm(int Char, vector <CQPar *> *ProcPar, vector <double> Eqm, vector <CPar *> *StateProbs);	// Constructor
	// Constructor functions to create complex covarion equilibriums (>1 observable distribution for all hidden states)
	CCovEqm(int Char, vector <CQPar *> *ProcPar, vector <vector <double> > Eqm, vector <vector <int> > EqmMap, vector <CPar *> *StateProbs, bool AllowOpt = false);
	// Destructor function
	~CCovEqm();															// Destructor
	// Interaction
	vector <double> Eqm();								// Returns equilibrium distribution for the whole process
	vector <double> SubEqm(int HiddenState);			// Returns the equilibrium distribution for a specific hidden state (size == m_iDataChar);
	vector <double> TransEqm();							// Returns the equilibrium distribution for the C matrix (describes transitions between hidden states; size == m_iNoCovState)
	void Shuffle();										// Shuffle the states that the eqm applies to
	void ResetEqm(vector <double> Eqm, bool RandomBit);	// Resets the eqm
private:
	// General variables
	int m_iNoCovState;								// Number of covarion states
	int m_iDataChar;									// Nunber of characters in the data (m_iChar == m_iNoCovStates * m_iDataChar)
	vector <CPar *> m_vpStateProbs;						// Probability of each of the m_iDataChar states
	// Variables relating to observed equilibrium distributions (i.e. those of characters)
	int m_iNoCovObsEqm;									// Number of different equilibrium distributions for the hidden-states
	vector <vector <int> > m_vviEqm2State;				// Mapping of these eqm distributions onto hidden states
	// Functions
	void ConstructSimpleCovEqm(vector <double> eqm,vector <CPar *> *StateProbs);	// Construct the equilibrium
	void AdapterFunc();									// Applies the simple eqm distribution to Q
};

///////////////////////////////////////////////////////////////////////////
// Class defining the space for an individual site at a node
///////////////////////////////////////////////////////////////////////////
class CSite {
public:
	CSite(int *Char);						// Basic constructor
	CSite(const CSite &Site);						// Copy constructor
	~CSite();								// Basic destructor
	// Interaction functions
	inline void Overwrite(CSite &Site,int SiteNum);	// Function that makes a site a copy of another
	inline void CopyVals(CSite *Site, bool ForceReal = false);	// Copies values in Site->m_pSpacePointer and Site->m_pScalePointer to this
	inline bool IsReal()	{ return m_bReal; }		// Whether computations should be performed for a site
	inline double *Sp() { return m_pSpacePointer;}	// Returns the space for the site
	inline double *RealSp() { return m_ardSpace; }	// Returns the real space regardless (for derivative calculations)
	inline int *Sc()	 { return m_pScalePointer;}	// Returns the space for the scaling factor
	inline int *RealSc() { return &m_iScale; }		// Returns the real scaling factor regardless
	inline int OriSite()	{ return m_iOriSite; }	// Returns the site a sequence was based on
	inline int Copies()	{ return m_iCopyNum; }	// Returns the number of copies of a site
	inline void ResetSite();						// Resets the site to a real site
	inline void ZeroSite(bool ForceAll = true);	// Sets all values in a site to zero
	friend ostream &operator<<(ostream &os, const CSite &Site);	// Output function
private:
	// Variables
	int *m_Char;							// The size
	double *m_pSpacePointer;				// Pointer to the space used
	double *m_ardSpace;						// The space for the node calculations
	int m_iScale;							// Value holding the scaling factor
	int *m_pScalePointer;					// The Scaling factor
	int m_iOriSite;							// if(m_bReal == false) The site it has been copied from
	int m_iCopyNum;							// if(m_bReal == true) the number of copies of the site
	bool m_bReal;							// Whether the node is actually real
	CSite *m_pOrSite;						// Pointer to original site
};

// This is the base process class
/////////////////////////////////////////////////////////////////////////
// Provides basic functionality for all likelihood calculations
// It will include information held by the individual processes, including:
//
//	- Q matrix information
//  - Partial likelihood information
//
// The process will
//	- Scale the process to mean rate 1 or leave it unscaled (for special cases)
//
// The model will expect to be able to:
//  - call the process calculation functions and receive back a final vector of likelihoods

class CBaseProcess {
public:
	// Constructor
	CBaseProcess(CData *Data, CTree *Tree,string Name = "Unnamed process");
	// Destructor
	~CBaseProcess();

	////////////////////////////////////////////////////////////////////////////////////////
	// Main interaction functions for calculations and so on
	////////////////////////////////////////////////////////////////////////////////////////
	void PrepareFastCalc(vector <int> *C = NULL);					// Prepares the columns for fast calculations from vector (if blank it makes its own
	void DecompressSpace();																// Reverts to normal memory allocation

	bool PrepareLikelihood(bool DoQ, bool ForceRemake);		// Prepare for the main likelihood function
	bool CreatePTMats(int Bra = -1);				// Creates the PT matrices ready for for branch [Bra]; if(Bra == -1) Do all branches...
	bool Likelihood(bool ForceReal = false);		// Main likelihood function
	virtual double Penalty() { return 0.0; }		// Returns a penalty value for penalized likelihood (if appropriate)
	double LogL();													// Calculate log likelihood of process (used after Likelihood is called);
	CProb &Lsum(int Site);						// Get the sitewise likelihood

	// Data and tree memory manipulation routines
	void ApplyNewData(CData *D, bool RedoSpace = true);	// Apply new data to the process (RedoSpace will delete existing space and replace it)
	void ApplyNewTree(CTree *T);						// Apply new tree to the process (RedoSpace will delete existing space and replace it)
	void CleanCalcSpace();											// Clean space used by calculatios
	void MakeCalcSpace(bool AllowRemake = false);					// Make the calculation space (or remake it if necessary) from character size

	// Functions for getting branch derivatives
	virtual void PrepareBraDer();						// Prepare the process for get branch derivatives function
	bool GetBraDer(CProb *ModelL);				// Gets the branch derivatives. ModelL = array of partial likelihoods for whole model


	// Parameter access functions
	////////////////////////////////////////////////////////////////////////////////////////
	CQPar * Kappa(EDataType Type);					// Adds kappa to a DNA (4x4) or Codon (64x64) process
	int Size() { return m_iSize; }						// Returns the size
	int PatOcc(int i) { assert(m_pData != NULL); assert(InRange(i,0,m_pData->m_iSize)); return m_pData->m_ariPatOcc[i]; }
	inline string Name()	{ return m_sName; }
	inline int NoPar()		{ return (int)m_vpPar.size(); }			// Returns the number of parameters
	inline int Char() { return m_iChar; };				// Returns the number of states in the model
	inline int HiddenChar() { return m_iHiddenChar; }	// Returns the number of hidden states in the process
	inline int NoSeq() { return m_pData->m_iNoSeq; }	// Returns the number of sequences in the process
	inline int DataChar() { return NumStates(m_DataType); }	// Returns the data character
	bool PseudoProcess() { return m_bPseudoProcess; }	// Returns whether process is a pseudoprocess
	CQPar *RatePar() { return m_pRate; }			// Returns a pointer to the rate parameter
	CQPar *ProbPar() { return m_pProcProb; }		// Returns a pointer to the probability parameters
	CPar *pPar(int i){ assert(i<NoPar()); return m_vpPar[i]; }	// Return a parameter
	CPar *AddQPar(CQPar *P) { assert(P!=NULL); m_vpPar.push_back(P); return P; }	// Adds a parameter
	CTree *Tree() { if(IsSubTree()) { return m_pSubTree; } return m_pTree; }	// Returns a pointer to the tree
	vector <CBaseEqm *> EqmPar() { return m_vpEqm; }			// Returns the Eqm parameters
	vector<double> Eqm(int QMatNum) { return m_vpQMat[QMatNum]->Eqm(); };	// Returns the numerical value of equilibrium
	virtual vector <double> RootEqm();				// Returns the root equilibrium for likelihood calculations
	void OutEqm(int QMat = -1, ostream &os = cout);	// Returns the eqm distribution
	double Rate(double NewRate = -1.0,bool MakeRateOpt = false);// Returns the rate and if NewRate >= 0 sets m_pRate = NewRate
	bool MaxRate() { return m_bMaxRate;}						// Whether the process is meant to be max rate
	bool MakeMaxRate() { m_bMaxRate = true; }					// Enforces max rate
	bool MakeNormalRate() { m_bMaxRate = false;	}  				// Goes back to a 'normal' rate
	CQPar *AddRatePar2Opt();									// Adds rate parameter to the optimised parameters
	void RemovePar(string Name);						// Removes a parameter from the process

	CTree *MainTree(CTree *T = NULL);							// Pointer to underlying tree
	CData *MainData() { return m_pData; }						// Pointer to underlying data
	inline bool IsSubTree() { if(m_pSubTree == NULL) { return false; } return true; }
	inline bool IsCompressed() { return m_bCompressedSpace; }

	// Functions dealing with the processes probability
	double Prob() { return m_pProcProb->Val() / *m_piProbFactor; }					// Returns the probability of the process
	void SetProbFactor(int i) { *m_piProbFactor = i; }

	// Functions for shuffling parameters
	virtual void ShuffleParameters() { cout << "\nDoing an empty shuffle..."; exit(-1); };				// Useful for mixture models where parameters combinations cause problems...

	// Some useful functions for debugging
	vector <double> GetQMatSums();						// Returns the sum of the off diagonal entries in the Q matrix (useful for detecting silly changes)

	// Functions for getting copies of the process (used in model.cxx)
	// These *MUST* be overloaded when other virtual functions are used
	virtual CBaseProcess *RateProcessCopy();								// Returns a pseudoprocess with a seperate rate category
	virtual CBaseProcess *GammaRateProcessCopy();							// Returns a pseudoprocess with a seperate rate category and tied probability (used for gamma distribution)

	// interaction functions for subtree assignment and calculations
	////////////////////////////////////////////////////////////////////////////////////////
	void ApplyCPMapping(vector <int> LeafMap, vector <int> NodeFr, bool DoPartial = true);	// Functions that applies CP mapping
	void CleanCPMapping();												// Function that cleans everything up and removes the partial likelihood mappings
	// Functions that allocate for a specific tree
	void ApplySubTree(CTree *Tree);							// Applies a subtree for likelihood computation
	void CleanSubTree();									// Clears subtree from likelihood computation
	bool OldGetBraDer(CTree *Tree, bool Init = true);			// Calculate branch derivatives (returns whether it was successful)
	// Function that delivers the likelihood of a site
	CProb &L(int Site)	{ assert(IsProb(Prob())); return m_ardL[Site].Multiply(Prob(),false); }; 			// The final likelihoods
	CProb &ModelL(int Site) { return m_arModelL[Site]; }
	bool LOkay() { if(MATCH_PAML == 1) { return true; } return FlipBool(m_bFailedL); };			// Returns whether the likelihood has computed
	bool IsGamma() { return m_bIsGamma;	}					// returns whether the process is gamma distributed
	void MakeGamma() { m_bIsGamma = true; }					// Flags that the process has a gamma distributed rate associated with it
	// Miscellaneous space access and output functions
	void OutQ(int QNum = -1,ostream &os = cout);								// Output the Q matrices for the process
	ostream &SiteOut(ostream &os, int Node, int PosBeg = -1, int PosEnd = -1);	// Output function to output details of a site
	void OutPT(ostream &os, int Branch);										// Output the PT for a branch
	void ZeroSpace();															// Function that removes all values currently held in the space
	vector <double> AggregateQ(int Position);				// Returns the aggregated Q matrix for the Position-th character
//	vector <double> AggregatePT(int Position);				// Returns the aggregated P(t) matrix for the Position-th character

	ostream &SiteL(int Site, ostream &os = cout);								// Output likelihood of a site
	void SetOutputDetail(bool V) { m_bDetailedOutput = V; }						// Set level of detail in general output
	void OutPartL(int Site, int Node, ostream &out = cout);						// Formatted output of partial likelihood for a node and site
	vector <double> GetPartL(int Site, int Node);								// Returns the partial likelihood for a node and a site
	int GetPartLScale(int Site, int Node);										// Returns the partial likelihood scale factor for a node and a site
	vector <double> GetPT(int Branch);											// Returns the PT matrix
	// Space updating routines
	virtual bool Make_PT(int Branch, bool RedoRate = false);				// Makes the P(t) matrix for the likelihood function
	void ScaleQ() { int i; FOR(i,(int) m_vpQMat.size()) { m_vpQMat[i]->ScaleQ(m_pRate->Val()); } };	// Scale Q matrices ready for calcs
	void LeafNode_Update(int NTo, int NFr, int Br, CTree *T, int First, bool DoCompleteUpdate = false);
	void BranNode_Update(int NTo, int NFr, int Br, CTree *T, int First, bool DoNTo = true, bool DoNFr = true, bool DoCompleteUpdate = false);
	// Space functions that calculate partial likelihoods from a centre point from SNAP paper
	void PreparePartialL(int Node, int NodeFrom,int Node2Trans=-1);	// Prepare a partial likelihood: if(Node2Trans != -1) { transfer space info to Node2Trans space; }
	// Function for FastBranchCalc
	void GetBranchPartL(CProb **P, int NTo, int NFr, int Br);	// Get the sitewise likelihoods
	// Reset the calculation space...
	void ResetCalcSpace();
protected:
	// Critical variables
	int m_iChar;						// Number of characters in the process
	int m_iHiddenChar;					// Number of hidden states in the process (only matters for THMMs)
	int m_iDataChar;					// Number of characters in the observable data (m in the model)
	int m_iSize;						// Length of the sequences
	int m_iSiCh;						// Is m_iChar * m_iLength
	int m_iChar2;						// Is m_iChar * m_iChar
	string m_sName;						// Name of process
	vector <CPar *> m_vpCovProbs;				// The probability parameters of each hidden state (size == m)
	EDataType m_DataType;				// The type of data the process is meant to apply to
	int m_iSpaceSize;					// Size allocated in space
	int m_iSpaceNoSeq;					// NoSeq allocated in space
	// Some flags of various things
	bool m_bPseudoProcess;				// Whether the process is real or not (i.e. needs no recomputing of QMat and m_vdEqm and only scales)
	bool m_bFailedL;					// Likelihood computations failed or is 0.0
	bool m_bIsGamma;					// Whether the process has a gamma distributed rate associated with it.
	bool m_bBraDerReady;				// Flag to check whether branch derivative calculations were ready
	bool m_bDetailedOutput;				// Flag as to whether to do detailed output
	bool m_bModelPerBranch;				// Flag highlighting whether parameters need estimating per branch
	bool m_bDoStandardDecompose;		// Whether eigen vectors/values for the model are produced in the usual manner

	////////////////////////////////////////////////////////////////////////////////////////
	//			Parameter related definitions - Including some standard model stuff
	//			Mostly used for preparing and defining models
	//
	// ********* NOTE: Parameters are held in two arrays: m_vpPar and m_vpEqmPar ***********
	////////////////////////////////////////////////////////////////////////////////////////

	// Some Q related values including parameters
	vector <CQMat *> m_vpQMat;			// Vector of the Q matrices used by the process (1 for most, but many when process varies across branches)
	vector <double> m_vdEqm;			// Pointer to the equilibrium distribution at the root (as used for likelihood computations)
	int m_iRootQ;						// The Q matrix of the process at the root
	string m_sABET;						// The string defining the processes alphabet
	vector <int> m_viQ2Bra;				// Vector describing which Q matrices describe which branches (if empty then use m_vpQMat[0] for all branches
	int m_iABET_length;					// The length of each character in the alphabet
	vector <CQPar *> m_vpPar;			// The parameters in the process
	CQPar *m_pRate;						// The rate parameter -- Usually set to one, but varies when rate varies between processes
	bool m_bMaxRate;					// Whether the process is meant to be max rate;
	CQPar *m_pProcProb;					// The probability of the process occurring
	int *m_piProbFactor;				// The number of processes a whole probability is divided between (important for gamma rates)
	double m_dBaseVal;					// The value that a Q is initialised to (doesn't matter is its scaled...)
	int m_iProcID;

	// Functions for creating Q matrices

	virtual bool PrepareQMats(vector <int> Qs2do, bool DoScale = true);			// Prepares the Q matrices in a process
	virtual bool PrepareQMats(bool DoScale = true)								// Overloaded function for doing all the Q matrices
				{ vector <int> v; return PrepareQMats(v,DoScale); }
	virtual void DoPostHocQMatUpdate() {};								// Function available if necessary

	// Variables and functions relating to the eqm distribution of Q matrices
	vector <CBaseEqm *> m_vpEqm;										// The equilibrium functions
	vector <double> SimpleEqm(int QMatID,bool AllowSumNotOne = false);	// Function returning a vector of the eqm distribution for a Q matrix

	// Functions relating to parameters and Q matrix
	void CleanPar();					// Removes all the parameters in m_arpPar and their mappings
	void CleanQ();						// Removes all Q matrices and their mappings
	void MakeABET(EDataType Type);		// Define it as a standard type
	int GetProcID() { return ProcID++; }// Gets the unique process ID;

	// Functions for adding Q matrices to the model
	virtual CQMat *Add_QMat(string Name,EDataType Type);					// Function for adding simple Q matrices based on standard data types
	virtual CQMat *Add_QMat(string Name, int Char);							// Function for adding Q matrices of size Char
	virtual CQMat *Add_CodRedQMat(string Name, int Char);			// Function for adding Q matrices for codon data with the genetic code enforced
	// Functions for mapping Q matrices back to the likelihood function
	int QMat4Bra(int Branch);										// Returns the number of the QMat object for a specific branch

	// The variables carrying the information about the process
	// Uses structure CSite for vector space
	bool m_bCheckSpace;					// Checks that the space has been initialised
	double *m_ardPT;					// The PT matrices
	double *m_ardQP;					// The Q * P(t) matrices

	bool m_bCompressedSpace;			// Whether the space has been compressed or not
	bool m_bDoingPartial;				// Set to true when making a partial likelihood (stops anything being transferred to m_iStartCalc);
	vector <CSite> m_vSpace;			// Space used for calculations
	vector <CSite> m_vBackSp;			// Space used for backwards calculations (for branch derivatives)

	CProb *m_ardL;						// The final likelihoods
	CProb *m_arModelL;					// Pointer to the final model likelihoods (set in PrepareBraDer and calculated in CBaseModel::Sitewise
	// Space access helpers
	inline int InitNodePos(int Node)	{ return (Node * m_iSize); }
	inline int PartLNode()				{ return m_pTree->NoNode(); }

	// Functions that define access to calculation space
	void CleanSpace();												// Cleans all the currently used space
	void MakeBasicSpace(int Char);									// Make the basic properties of the process
	void MakeBasicSpace(EDataType Type);							// Make basic properties under a certain standard data type

	// Space access functions
	inline double *PT(int Br) { return &m_ardPT[m_iChar2 * Br]; };					// The PT matrix for a branch
	inline double *QP(int Br) { return &m_ardQP[m_iChar2 * Br]; };					// The Q * PT matrix for a branch
	inline double *Fd(int Node,int Pos) { return QkFd(InitNodePos(Node) + Pos); };	// Get partial likelihood vector for specific node
	inline double *QkFd(int NodePos) { return m_vSpace[NodePos].Sp(); };			// Quick access function to Forward space

	inline bool FdReal(int Node, int Pos) {return QkFdReal(InitNodePos(Node)+Pos);};// Returns whether forward space needs calculating
	inline bool QkFdReal(int NodePos) { return m_vSpace[NodePos].IsReal(); };		// Quick checker to see if space needs calculating
	inline bool BkReal(int Node, int Pos) {return QkBkReal(InitNodePos(Node)+Pos);};// Returns whether back space needs calculating
	inline bool QkBkReal(int NodePos) { return m_vBackSp[NodePos].IsReal(); };		// Quick checker to see if back space needs calculating
	inline int OrFd(int NodePos) { return m_vSpace[NodePos].OriSite(); }			// Return the original label for a node in Fdspace
	inline int OrBk(int NodePos) { return m_vBackSp[NodePos].OriSite(); }			// Return the original lable for a node in the BkSpace
	inline int FdCopies(int Node, int Pos) { return QkFdCopies(InitNodePos(Node) + Pos); }	// Returns the number of copies of a site
	inline int QkFdCopies(int NodePos) { return m_vSpace[NodePos].Copies(); }				// Returns the number of copies of a site
	inline double *Bk(int Node,int Pos){ return QkBk(InitNodePos(Node) + Pos); };	// Get the Backwards likelihood vector for specific node
	inline double *QkBk(int NodePos) { return m_vBackSp[NodePos].Sp(); };			// Quick access function to backward space
	inline int *FdScale(int NodeNum, int Pos) { return QkFdSc(InitNodePos(NodeNum)+Pos); };	// Access to the forwards space held in m_ardSpace
	inline int *BkScale(int NodeNum, int Pos) { return QkBkSc(InitNodePos(NodeNum)+Pos); };	// Access to the backwards space held in m_ardBackSpace
	inline int *QkFdSc(int NodePos) { return m_vSpace[NodePos].Sc(); }; 			// Quick forward scale access
	inline int *QkBkSc(int NodePos)	{ return m_vBackSp[NodePos].Sc(); };			// Quick backward scale access
	// Force space functions
	inline double *ForceRealFd(int Node,int Pos) { return QkForceRealFd(InitNodePos(Node) + Pos); }	// Force use of real backward space (required for derivative calculations)
	inline double *QkForceRealFd(int NodePos) { return m_vSpace[NodePos].RealSp(); }
	inline int *ForceRealFdSc(int Node, int Pos) { return QkForceRealFdSc(InitNodePos(Node) + Pos); }
	inline int *QkForceRealFdSc(int NodePos) { return m_vSpace[NodePos].RealSc(); }
	inline double *ForceRealBk(int Node,int Pos) { return QkForceRealBk(InitNodePos(Node) + Pos); }	// Force use of real backward space (required for derivative calculations)
	inline double *QkForceRealBk(int NodePos) { return m_vBackSp[NodePos].RealSp(); }
	inline int *ForceRealBkSc(int Node, int Pos) { return QkForceRealBkSc(InitNodePos(Node) + Pos); }
	inline int *QkForceRealBkSc(int NodePos) { return m_vBackSp[NodePos].RealSc(); }

	void CopyNode(int NodeFr, int NTo);										// Copy all space in NodeFr to NTo
	void CleanScale(int NodeNum, bool ForceAll = true);						// Removes all scales in a node
	void TransScale(int NodeTo,int NodeFr,bool First, bool Partial = true, bool ForceReal = false);	// Transfer (sum) values from NodeFr to NodeTo
	void DoScale(int Node, bool ForceRealScale = false);					// Performs scaling on the fly; ForceRealScale = true will treat all nodes as real

	inline int *LScale(int Site) { return FdScale(PartLNode(),Site); }		// Final node and m_ardPerSitelnL scaling factor
	inline double *PartL(int Site) { return Fd(PartLNode(),Site); };		// The Partial likelihoods


	/////////////////////////////////////////////////////////////////////////////////////////
	// Variables and functions for calculating partial likelihoods upto subnodes, &c
	bool m_bSubTreeActive;									// Flag stating that the likelihood function is expecting to work with a subtree
	vector <int> m_viLeafMap;								// The mapping of the full data onto the subtrees nodes (0<i<m_iNoSeq) == Real sequences; -1 = Partial likelihood from Space i;
	CTree *m_pSubTree;										// The tree being used with the partial likelihood

	/////////////////////////////////////////////////////////////////////////////////////////
	// Likelihood computation functions
	// --------------------------------
	// These are more complicated than neccessary, but at least they work!

	// Normal calculations
	void PartialL(CTree *pTree, int iNoTo, int iNoFr, int Branch, bool LeafFirst = true, bool BlockWriteback = false);
	void BranchNodePartialL(int SpFrom, int SpTo, double *PT, bool First);
	void LeafNodePartialL(CTree *pTree, int LeafNode, int Branch, int SpFlag, double *PT,bool First);		// Usual leaf node calculation
	void MakeZeroRateLikelihood();				// Makes likelihood for Rate = 0 models
	void MakeMaxRateLikelihood();				// Makes likelihood for Rate = inf model (garbage collector models)

	// Branch derivative calculations
	void Branch_dT(int NTo, int NFr, int Branch, CTree *pTree, int First, int *BrError);
	void LeafNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError);
	void BranNode_dT(int NTo, int NFr, int Br, CTree *pTree, int First, int *BrError);
	double PartialGrad(int site,double Total,int SiteScale);

	// Some pointers to other data based other stuff
	CTree *m_pTree;				// Pointer to the tree the model should be using
	CData *m_pData;				// Pointer to the data

	// Some virtual functions that will vary considerably within the calculations
	virtual void Data2PartL(int Char, double *PT, double *RetSpace, vector <double> *eqm);
	virtual void Vec_by_Data(int Char, double *Vec);		// Multiplies a vector by the characters (i.e. set all values not of type Char to 0); Returns the sum of these values
	virtual double Sum_Vec(int Char, double *Vec, vector <double> eqm);	// Returns the sum of elements * eqm matching Char
	void *Make_PT(double T,double *PT,double *U,double *V,double *Root);	// Makes the P(t) matrix													// Makes the Q matrix
	// Output functions
	friend ostream &operator<<(ostream &os, CBaseProcess &Proc) { return Proc.Output(os); }
	virtual ostream &Output(ostream &os);

	///////////////////////////////////////////////////////////////////////////////////
	// Some useful functions for building the equilibrium distribution of models
	// AddSimpleEqm adds a simple eqm distribution with only a single eqm distribution.
	void AddSimpleEqm(bool Opt = true);										// Add a simple eqm distribution from sequences
	void AddSimpleEqm(int Char, vector <double> Freq, bool Opt = true);		// Add a simple eqm distribution from Freq
	void AddSimpleCovEqm(vector <CPar *> *StateProbs, bool ForceEqu = false,bool Opt = false);						// Add a simple eqm for Covarion models
	void AddSimpleCovEqm(int Char, vector <double> Freq, vector <CPar *> *StateProbs, bool ForceEqu = false, bool Opt = false);	// Add a simple eqm for Covarion models
	// AddComplexEqm adds multiple eqm distributions
	// 1. Add seperate eqm for each hidden state. The initial eqm distributions are guessed from looking at the data
	void AddComplexCovEqm(int Char, vector <CPar *> *StateProbs, bool AllowOpt = true);
	// 2. Add seperate eqm for each hidden state. The initial eqm distributions are taken from Freq
	void AddComplexCovEqm(int Char, vector <vector <double> > Freq, vector <CPar *> *StateProbs, bool AllowOpt);
	// 3. Equilibrium distribution specified in Freq (the actual frequencies) and FreqMap (the list of states that each vector in Freq maps to)
	void AddComplexCovEqm(int Char, vector <vector <double> > Freq, vector <vector <int> > FreqMap,vector <CPar *> *StateProbs, bool AllowOpt); // Adds an equilibrium distribution with different probabilities for each CS state
	// AddCodonEqm adds a codon equilibrium distribution
	void AddCodonEqm(int GenCode, int Char, ECodonEqm CE, bool Opt);
};


///////////////////////////////////////////////////////////////////////////////////////
// Some basic process types that the processes can be assigned

enum DNAProc	{pRY,pJC,pFEL,pK2P,pHKY,pREV,pCOV,pDNA_THMM,pBRACOV,DNA_OTHER};
enum AAProc		{pEQU,pJTT,pWAG,pDAY,pMTREV,AA_OTHER};
enum CodonProc  {pM0};

class CDNAProcess : public CBaseProcess	{
public:
	CDNAProcess(CData *Data, CTree *Tree, DNAProc Model, string name);		// Constructor to produce a defaulted model
	~CDNAProcess();
	// Specific routines for building different types of model
	vector <CQPar *> MakeDNAREV(EDataType Type);		// Make a full reversible model
};

class CAAProcess : public CBaseProcess {
public:
	CAAProcess(CData *Data, CTree *Tree, AAProc Model, bool AddF = true);								// Constructor to produce a defaulted model
	CAAProcess(CData *D, CTree *T, string Name, bool AddF, double *S_ij, double *pi_j);			// Constructor to produce a general empirical model
	~CAAProcess();
	// General empirical model creation routine
	void CreateEMPmodel(double *S_ij,double *Freq,bool DoF);
	// Model specific routines
	void MakeEQU(bool AddF);
	void MakeDAY(bool AddF);
	void MakeWAG(bool AddF);
	void MakeJTT(bool AddF);
	void MakeMTREV(bool AddF);
};



class CCodonProcess : public CBaseProcess {
public:
	// Constructor function
	CCodonProcess(CData *Data, CTree *Tree, CodonProc Model,ECodonEqm CE,int GenCode = 0);
	// Destructor function
	~CCodonProcess();
	// Parameter specific bits and pieces
	CQPar * AddOmega(int GenCode);
	void AddMultiChangeZeros(int GenCode);
	// Functions for getting overall rates of certain types of events
	double ObsSynRate();					// Gets the observed rate of synonymous substitutions
	double ObsNonsynRate();					// Gets the observed rate of non-synonymous substitutions
	// Other interaction functions
	int GenCode() { return m_iGenCode; }
	// Output
	ostream &Output(ostream &os);
private:
	int	m_iGenCode;							// Stores the genetic code
};


#endif
