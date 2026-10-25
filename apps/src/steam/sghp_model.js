/* SGHP cycle model in JavaScript: a port of HeatPumpLab ConventionalHeatPump (IHX, compressor
   efficiency map, UA shell heat loss, compressor speed). Pure functions over a property object
   with PT/PH/PS/PQ/TQ(…) -> {T (C), P (MPa), h (kJ/kg), s (kJ/kg/K), D (kg/m3), q, phase}.
   Units: MPa, C, kJ/kg, kg/s, kW — the same as the Python model. */
(function (root) {
  "use strict";

  // ---------------------------------------------------------------- compressor (R600a fits)
  // Semi-hermetic reciprocating compressors, R600a fits (compressor_data.py), smallest first
  const COEFFS = {
    "Recip-S": {
      vol: [1.12984011409572, -0.01000517868509311, 9.959408735935163e-05, -5.429511239235132e-07,
        0.003567133668441818, 9.479581042021442e-05, -9.565274823906224e-07, 7.170717115542697e-09,
        -0.0001232257031541566, 1.028269483251338e-06, -2.282557219145039e-08, 7.981615275880793e-11,
        2.550844912615862e-06, -6.231819262761143e-08, 8.149157481312754e-10, -3.246732661005363e-12],
      is: [0.4116836162879849, 0.005039822811869186, -2.631708318422278e-05, -2.046577662134925e-08,
        -0.003314441060494978, 0.0001339324043644695, -8.108428787902513e-07, 2.761876674705595e-09,
        -0.0008637056504366176, 3.207454510381633e-05, -4.196854713876601e-07, 1.713676254812856e-09,
        4.34014145040436e-06, -3.48778095942875e-07, 5.507361160638747e-09, -2.42783380539514e-11],
    },
    "Recip-M": {
      vol: [1.096611343131835, -0.00726227358913289, 5.309047868640223e-05, -3.096676393991497e-07,
        0.006733364512289938, -3.458999501344552e-05, 1.169742469456288e-06, -5.714608237170313e-09,
        -0.0002840846510202909, 6.967940930769206e-06, -9.7576674515825e-08, 4.673105300527226e-10,
        5.100400187957733e-06, -1.539976618200062e-07, 1.904689442801216e-09, -8.23309880717932e-12],
      is: [0.2820090966150441, 0.0120891735264767, -0.0001194130317419781, 3.921449567398475e-07,
        0.001020597758718241, -0.0001552582461411895, 3.819292261440514e-06, -2.125798153420446e-08,
        -0.001021385474733074, 4.213732121164092e-05, -5.832829781115806e-07, 2.561191347258478e-09,
        8.812924593464052e-06, -5.206597937866331e-07, 7.768970794352923e-09, -3.446962012422449e-11],
    },
    "Recip-L": {
      vol: [0.9456698645061833, -0.0009241333917138055, -1.981682414911695e-05, -2.974086664626571e-08,
        0.006439330476049707, -0.0001232201170140139, 2.390487484843521e-06, -8.830117868237834e-09,
        -0.0001643829586845571, 5.825566412853439e-06, -9.913209075311744e-08, 4.150372919033212e-10,
        2.535938980404551e-06, -9.615602095713242e-08, 1.436561992150825e-09, -6.097111927442945e-12],
      is: [0.365732694101855, 0.008380402635408853, -7.40105429843498e-05, 1.995263626032545e-07,
        -0.0006020506345257311, -2.396568897230092e-05, 1.303701754995347e-06, -6.770570118114042e-09,
        -0.0009300664051148016, 3.286808772249003e-05, -3.98272290261303e-07, 1.544452281550714e-09,
        1.954486640420136e-06, -2.661964600946924e-07, 4.30163369026706e-09, -1.872496284512809e-11],
    },
  };
  // swept volume at 50 Hz (1450 rpm), m3/h
  const DISPLACEMENT_50HZ = { "Recip-S": 13.42, "Recip-M": 56.25, "Recip-L": 151.6 };
  // ---------------------------------------------------------------- water and steam (steam from bled hot water)
  // IAPWS-95 saturation (CoolProp 8.0.0 "Water"), 1 K grid 1-220 C: [T C, P MPa, h_f, h_g kJ/kg, v_f m3/kg];
  // the same table as water_steam.py
  const WATER_SAT = [[1,0.0006570856,4.1767,2502.7312,0.001000149],[2,0.0007059856,8.3918,2504.5650,0.001000107],[3,0.0007580805,12.6039,2506.3980,0.001000083],[4,0.0008135484,16.8132,2508.2303,0.001000075],[5,0.0008725751,21.0200,2510.0619,0.001000083],[6,0.0009353551,25.2243,2511.8928,0.001000106],[7,0.001002091,29.4265,2513.7229,0.001000145],[8,0.001072995,33.6266,2515.5522,0.001000197],[9,0.001148288,37.8248,2517.3807,0.001000265],[10,0.001228199,42.0213,2519.2083,0.001000345],[11,0.001312969,46.2161,2521.0351,0.00100044],[12,0.001402848,50.4094,2522.8610,0.001000547],[13,0.001498096,54.6013,2524.6860,0.001000667],[14,0.001598984,58.7920,2526.5101,0.0010008],[15,0.001705793,62.9815,2528.3332,0.001000945],[16,0.001818816,67.1698,2530.1553,0.001001101],[17,0.001938358,71.3572,2531.9764,0.00100127],[18,0.002064735,75.5437,2533.7965,0.001001449],[19,0.002198275,79.7293,2535.6155,0.00100164],[20,0.002339318,83.9141,2537.4334,0.001001842],[21,0.002488219,88.0983,2539.2502,0.001002054],[22,0.002645344,92.2818,2541.0658,0.001002277],[23,0.002811072,96.4648,2542.8802,0.00100251],[24,0.002985798,100.6472,2544.6934,0.001002753],[25,0.003169929,104.8292,2546.5054,0.001003006],[26,0.003363888,109.0108,2548.3161,0.001003268],[27,0.003568112,113.1920,2550.1254,0.001003541],[28,0.003783053,117.3729,2551.9334,0.001003822],[29,0.004009178,121.5536,2553.7400,0.001004113],[30,0.004246971,125.7340,2555.5452,0.001004413],[31,0.004496931,129.9142,2557.3489,0.001004722],[32,0.004759574,134.0943,2559.1511,0.00100504],[33,0.005035433,138.2742,2560.9518,0.001005367],[34,0.005325058,142.4541,2562.7509,0.001005702],[35,0.005629016,146.6339,2564.5484,0.001006045],[36,0.005947893,150.8136,2566.3443,0.001006397],[37,0.006282292,154.9934,2568.1384,0.001006757],[38,0.006632835,159.1732,2569.9309,0.001007126],[39,0.007000164,163.3531,2571.7215,0.001007502],[40,0.007384938,167.5330,2573.5103,0.001007887],[41,0.007787838,171.7131,2575.2973,0.001008279],[42,0.008209563,175.8933,2577.0823,0.001008679],[43,0.008650834,180.0736,2578.8655,0.001009087],[44,0.009112392,184.2541,2580.6466,0.001009502],[45,0.009594999,188.4348,2582.4257,0.001009925],[46,0.01009944,192.6157,2584.2027,0.001010355],[47,0.01062651,196.7968,2585.9775,0.001010793],[48,0.01117705,200.9782,2587.7502,0.001011238],[49,0.01175191,205.1598,2589.5207,0.00101169],[50,0.01235195,209.3418,2591.2889,0.00101215],[51,0.01297807,213.5240,2593.0548,0.001012616],[52,0.01363119,217.7065,2594.8183,0.00101309],[53,0.01431225,221.8894,2596.5794,0.001013571],[54,0.01502223,226.0727,2598.3380,0.001014058],[55,0.0157621,230.2563,2600.0941,0.001014553],[56,0.01653289,234.4402,2601.8476,0.001015055],[57,0.01733564,238.6246,2603.5985,0.001015563],[58,0.01817142,242.8094,2605.3467,0.001016078],[59,0.01904131,246.9947,2607.0922,0.0010166],[60,0.01994643,251.1804,2608.8349,0.001017128],[61,0.02088794,255.3665,2610.5747,0.001017663],[62,0.021867,259.5531,2612.3117,0.001018205],[63,0.0228848,263.7402,2614.0457,0.001018753],[64,0.02394259,267.9279,2615.7767,0.001019308],[65,0.0250416,272.1160,2617.5047,0.00101987],[66,0.02618312,276.3047,2619.2295,0.001020438],[67,0.02736846,280.4939,2620.9511,0.001021012],[68,0.02859896,284.6838,2622.6695,0.001021593],[69,0.02987599,288.8742,2624.3846,0.00102218],[70,0.03120093,293.0652,2626.0964,0.001022773],[71,0.03257522,297.2568,2627.8047,0.001023373],[72,0.03400031,301.4491,2629.5096,0.001023979],[73,0.03547769,305.6420,2631.2109,0.001024592],[74,0.03700887,309.8356,2632.9087,0.001025211],[75,0.0385954,314.0299,2634.6027,0.001025836],[76,0.04023884,318.2249,2636.2931,0.001026467],[77,0.04194082,322.4207,2637.9796,0.001027104],[78,0.04370297,326.6171,2639.6623,0.001027748],[79,0.04552695,330.8143,2641.3411,0.001028398],[80,0.04741447,335.0124,2643.0159,0.001029054],[81,0.04936727,339.2112,2644.6866,0.001029716],[82,0.0513871,343.4108,2646.3533,0.001030385],[83,0.05347577,347.6112,2648.0157,0.001031059],[84,0.05563511,351.8125,2649.6739,0.00103174],[85,0.05786697,356.0147,2651.3278,0.001032427],[86,0.06017326,360.2178,2652.9773,0.00103312],[87,0.06255591,364.4218,2654.6223,0.001033819],[88,0.06501687,368.6267,2656.2628,0.001034524],[89,0.06755815,372.8325,2657.8988,0.001035235],[90,0.07018177,377.0394,2659.5300,0.001035952],[91,0.07288979,381.2472,2661.1565,0.001036676],[92,0.07568433,385.4561,2662.7782,0.001037405],[93,0.0785675,389.6660,2664.3951,0.001038141],[94,0.08154148,393.8769,2666.0069,0.001038883],[95,0.08460847,398.0890,2667.6138,0.001039631],[96,0.0877707,402.3021,2669.2155,0.001040384],[97,0.09103044,406.5164,2670.8120,0.001041144],[98,0.09439,410.7318,2672.4033,0.001041911],[99,0.09785173,414.9484,2673.9893,0.001042683],[100,0.101418,419.1662,2675.5699,0.001043461],[101,0.1050912,423.3852,2677.1450,0.001044246],[102,0.1088738,427.6054,2678.7145,0.001045036],[103,0.1127683,431.8269,2680.2784,0.001045833],[104,0.1167772,436.0498,2681.8366,0.001046636],[105,0.1209031,440.2739,2683.3890,0.001047445],[106,0.1251485,444.4994,2684.9355,0.00104826],[107,0.1295161,448.7262,2686.4761,0.001049081],[108,0.1340085,452.9544,2688.0106,0.001049909],[109,0.1386285,457.1841,2689.5391,0.001050742],[110,0.1433787,461.4152,2691.0613,0.001051582],[111,0.148262,465.6478,2692.5773,0.001052428],[112,0.1532811,469.8818,2694.0870,0.001053281],[113,0.158439,474.1174,2695.5902,0.001054139],[114,0.1637384,478.3546,2697.0869,0.001055004],[115,0.1691824,482.5933,2698.5770,0.001055875],[116,0.1747738,486.8337,2700.0604,0.001056753],[117,0.1805156,491.0756,2701.5371,0.001057636],[118,0.186411,495.3193,2703.0069,0.001058526],[119,0.1924629,499.5646,2704.4699,0.001059423],[120,0.1986744,503.8117,2705.9258,0.001060326],[121,0.2050488,508.0605,2707.3746,0.001061235],[122,0.2115891,512.3111,2708.8162,0.00106215],[123,0.2182986,516.5635,2710.2506,0.001063072],[124,0.2251805,520.8178,2711.6777,0.001064001],[125,0.2322382,525.0739,2713.0973,0.001064936],[126,0.2394749,529.3319,2714.5094,0.001065877],[127,0.246894,533.5919,2715.9139,0.001066825],[128,0.2544989,537.8538,2717.3108,0.00106778],[129,0.2622931,542.1177,2718.6999,0.001068741],[130,0.27028,546.3836,2720.0811,0.001069709],[131,0.2784631,550.6516,2721.4544,0.001070683],[132,0.286846,554.9217,2722.8197,0.001071664],[133,0.2954323,559.1939,2724.1769,0.001072652],[134,0.3042255,563.4683,2725.5259,0.001073646],[135,0.3132294,567.7449,2726.8666,0.001074648],[136,0.3224476,572.0236,2728.1989,0.001075656],[137,0.3318839,576.3047,2729.5228,0.001076671],[138,0.3415421,580.5880,2730.8382,0.001077693],[139,0.3514258,584.8737,2732.1449,0.001078721],[140,0.3615391,589.1617,2733.4430,0.001079757],[141,0.3718857,593.4521,2734.7322,0.001080799],[142,0.3824696,597.7449,2736.0126,0.001081849],[143,0.3932947,602.0402,2737.2840,0.001082906],[144,0.4043651,606.3381,2738.5463,0.00108397],[145,0.4156846,610.6384,2739.7996,0.00108504],[146,0.4272575,614.9413,2741.0435,0.001086118],[147,0.4390878,619.2468,2742.2782,0.001087204],[148,0.4511796,623.5550,2743.5035,0.001088296],[149,0.4635371,627.8659,2744.7193,0.001089396],[150,0.4761645,632.1794,2745.9255,0.001090503],[151,0.4890661,636.4958,2747.1220,0.001091618],[152,0.5022461,640.8149,2748.3088,0.00109274],[153,0.5157088,645.1369,2749.4857,0.001093869],[154,0.5294587,649.4617,2750.6527,0.001095006],[155,0.5435,653.7895,2751.8097,0.001096151],[156,0.5578372,658.1202,2752.9566,0.001097303],[157,0.5724747,662.4539,2754.0933,0.001098463],[158,0.5874171,666.7907,2755.2197,0.001099631],[159,0.6026689,671.1305,2756.3358,0.001100806],[160,0.6182346,675.4735,2757.4413,0.00110199],[161,0.6341189,679.8196,2758.5364,0.001103181],[162,0.6503263,684.1689,2759.6208,0.00110438],[163,0.6668615,688.5214,2760.6944,0.001105587],[164,0.6837293,692.8773,2761.7572,0.001106802],[165,0.7009344,697.2364,2762.8091,0.001108026],[166,0.7184815,701.5990,2763.8500,0.001109257],[167,0.7363754,705.9649,2764.8799,0.001110497],[168,0.7546211,710.3343,2765.8985,0.001111745],[169,0.7732233,714.7072,2766.9058,0.001113002],[170,0.792187,719.0837,2767.9018,0.001114267],[171,0.8115171,723.4638,2768.8863,0.001115541],[172,0.8312187,727.8475,2769.8592,0.001116823],[173,0.8512967,732.2348,2770.8205,0.001118113],[174,0.8717561,736.6260,2771.7700,0.001119413],[175,0.8926021,741.0209,2772.7076,0.001120721],[176,0.9138398,745.4196,2773.6334,0.001122038],[177,0.9354742,749.8222,2774.5471,0.001123364],[178,0.9575107,754.2288,2775.4486,0.001124699],[179,0.9799544,758.6393,2776.3379,0.001126044],[180,1.002811,763.0538,2777.2149,0.001127397],[181,1.026084,767.4724,2778.0795,0.00112876],[182,1.049781,771.8952,2778.9315,0.001130132],[183,1.073907,776.3221,2779.7709,0.001131513],[184,1.098466,780.7533,2780.5976,0.001132904],[185,1.123464,785.1887,2781.4114,0.001134304],[186,1.148907,789.6285,2782.2123,0.001135715],[187,1.1748,794.0726,2783.0002,0.001137135],[188,1.201149,798.5212,2783.7749,0.001138564],[189,1.227959,802.9743,2784.5363,0.001140004],[190,1.255236,807.4320,2785.2844,0.001141454],[191,1.282985,811.8942,2786.0191,0.001142914],[192,1.311213,816.3611,2786.7402,0.001144384],[193,1.339924,820.8328,2787.4476,0.001145865],[194,1.369124,825.3092,2788.1412,0.001147356],[195,1.398819,829.7904,2788.8209,0.001148857],[196,1.429016,834.2766,2789.4866,0.001150369],[197,1.459719,838.7677,2790.1382,0.001151892],[198,1.490935,843.2638,2790.7755,0.001153426],[199,1.522669,847.7650,2791.3985,0.001154971],[200,1.554928,852.2713,2792.0070,0.001156526],[201,1.587717,856.7828,2792.6010,0.001158094],[202,1.621042,861.2996,2793.1802,0.001159672],[203,1.654909,865.8217,2793.7446,0.001161262],[204,1.689325,870.3492,2794.2940,0.001162863],[205,1.724295,874.8822,2794.8284,0.001164476],[206,1.759826,879.4207,2795.3476,0.001166101],[207,1.795923,883.9648,2795.8514,0.001167738],[208,1.832592,888.5145,2796.3398,0.001169387],[209,1.869841,893.0700,2796.8127,0.001171048],[210,1.907675,897.6312,2797.2698,0.001172721],[211,1.9461,902.1984,2797.7110,0.001174407],[212,1.985123,906.7714,2798.1363,0.001176106],[213,2.02475,911.3505,2798.5455,0.001177817],[214,2.064987,915.9357,2798.9385,0.001179541],[215,2.105841,920.5270,2799.3150,0.001181278],[216,2.147318,925.1246,2799.6751,0.001183029],[217,2.189424,929.7284,2800.0184,0.001184793],[218,2.232167,934.3387,2800.3450,0.00118657],[219,2.275552,938.9554,2800.6546,0.001188361],[220,2.319586,943.5787,2800.9471,0.001190166]];
  function waterSat(T) {
    const n = WATER_SAT.length;
    if (!(T >= WATER_SAT[0][0] && T <= WATER_SAT[n - 1][0])) throw new Error(`water saturation table covers 1-220 C, not ${T} C`);
    const i = Math.min(Math.max(Math.floor(T) - 1, 0), n - 2), a = WATER_SAT[i], b = WATER_SAT[i + 1], t = (T - a[0]) / (b[0] - a[0]);
    const f = k => a[k] + t * (b[k] - a[k]);
    return { P: f(1), h_f: f(2), h_g: f(3), v_f: f(4) };
  }
  // saturated steam at T_steam raised from feedwater at T_feed (pumped to the steam pressure) with duty Q_kW
  function steamFromFeedwater(Q_kW, T_steam, T_feed) {
    if (!(T_feed < T_steam)) throw new Error(`feedwater at ${T_feed} C is not below the steam temperature ${T_steam} C`);
    const s = waterSat(T_steam), w = waterSat(T_feed);
    const h_feed = w.h_f + w.v_f * (s.P - w.P) * 1000, dh = s.h_g - h_feed;
    return { m_kgs: Q_kW / dh, dh_kJkg: dh, h_feed_kJkg: h_feed, h_steam_kJkg: s.h_g, P_steam_MPa: s.P, preheat_kJkg: s.h_f - h_feed };
  }
  const ENVELOPE = [[-10, 30], [27, 42], [48, 70], [48, 120], [6, 120], [-10, 97]];

  function bicubic(a, to, tc) {
    let total = 0;
    for (let i = 0; i < 4; i++) {
      const row = a[4 * i] + tc * (a[4 * i + 1] + tc * (a[4 * i + 2] + tc * a[4 * i + 3]));
      total += row * Math.pow(to, i);
    }
    return total;
  }

  function inEnvelope(to, tc, env = ENVELOPE) {
    let inside = false;
    for (let i = 0, j = env.length - 1; i < env.length; j = i++) {
      const [x1, y1] = env[i], [x2, y2] = env[j];
      const cross = (x2 - x1) * (tc - y1) - (y2 - y1) * (to - x1);
      if (Math.abs(cross) < 1e-9 && to >= Math.min(x1, x2) - 1e-9 && to <= Math.max(x1, x2) + 1e-9 &&
          tc >= Math.min(y1, y2) - 1e-9 && tc <= Math.max(y1, y2) + 1e-9) return true;
      if ((y1 > tc) !== (y2 > tc) && to < x1 + (tc - y1) * (x2 - x1) / (y2 - y1)) inside = !inside;
    }
    return inside;
  }

  // ---------------------------------------------------------------- density basis (density_basis.py)
  // R600a saturated vapour density (kg/m3) against T (C), CoolProp 8.0.0 HEOS IsoButane; monotone
  // cubic (PCHIP) in ln rho. The maps are read at the R600a (to*, tc*) with the working fluid's
  // saturated-vapour densities at the evaporating and condensing pressures.
  const R600A_SAT_VAPOUR = satVapourDensity([-50, -49, -48, -47, -46, -45, -44, -43, -42, -41, -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 125.5, 126, 126.5, 127, 127.5, 128, 128.5, 129, 129.5, 130, 130.5, 131, 131.5, 132, 132.25, 132.5, 132.75, 133, 133.25, 133.5, 133.75, 134, 134.25, 134.5],
    [0.53214474, 0.56059848, 0.5902529, 0.62114383, 0.65330766, 0.68678135, 0.72160239, 0.75780885, 0.79543935, 0.83453306, 0.87512968, 0.9172695, 0.96099334, 1.0063426, 1.0533591, 1.1020854, 1.1525646, 1.2048401, 1.2589562, 1.3149575, 1.3728892, 1.4327973, 1.4947279, 1.5587282, 1.6248456, 1.6931281, 1.7636244, 1.8363838, 1.9114561, 1.9888916, 2.0687415, 2.1510572, 2.2358911, 2.323296, 2.4133254, 2.5060334, 2.6014747, 2.6997047, 2.8007796, 2.9047561, 3.0116916, 3.1216443, 3.2346731, 3.3508375, 3.4701978, 3.5928153, 3.7187517, 3.8480697, 3.9808328, 4.1171053, 4.2569524, 4.4004401, 4.5476352, 4.6986057, 4.8534203, 5.0121488, 5.1748618, 5.3416312, 5.5125297, 5.6876312, 5.8670106, 6.050744, 6.2389089, 6.4315835, 6.6288477, 6.8307824, 7.0374701, 7.2489944, 7.4654404, 7.6868946, 7.9134452, 8.1451817, 8.3821953, 8.6245789, 8.8724271, 9.1258362, 9.3849044, 9.6497318, 9.9204203, 10.197074, 10.479799, 10.768704, 11.0639, 11.365498, 11.673616, 11.98837, 12.309882, 12.638274, 12.973673, 13.316208, 13.666011, 14.023218, 14.387967, 14.7604, 15.140664, 15.528908, 15.925286, 16.329955, 16.743078, 17.164821, 17.595355, 18.034858, 18.483509, 18.941496, 19.409011, 19.886253, 20.373426, 20.870741, 21.378416, 21.896675, 22.425752, 22.965885, 23.517324, 24.080325, 24.655155, 25.242088, 25.841411, 26.453418, 27.078418, 27.716728, 28.36868, 29.034617, 29.714896, 30.409889, 31.119982, 31.84558, 32.587101, 33.344984, 34.119686, 34.911684, 35.721477, 36.549588, 37.396563, 38.262976, 39.149427, 40.056548, 40.985003, 41.935488, 42.90874, 43.905533, 44.926686, 45.973065, 47.045585, 48.145217, 49.27299, 50.430001, 51.617414, 52.836473, 54.088505, 55.374931, 56.697274, 58.057169, 59.456379, 60.896802, 62.380493, 63.909679, 65.486782, 67.114439, 68.795537, 70.533241, 72.331035, 74.192771, 76.122721, 78.125648, 80.206886, 82.372443, 84.629123, 86.984684, 89.448035, 92.029492, 94.741105, 97.597099, 100.61446, 103.81377, 107.22031, 110.86578, 112.79036, 114.79079, 116.87418, 119.0488, 121.3244, 123.71257, 126.22725, 128.88551, 131.70853, 134.72319, 137.96442, 141.47894, 145.33162, 149.61695, 151.96461, 154.48162, 157.20256, 160.17507, 163.46775, 167.18532, 171.50123, 176.73987, 183.64597, 194.90209], "R600a");

  function satVapourDensity(T, rho, fluid) {
    const n = T.length, y = rho.map(Math.log), h = [], d = [], m = new Array(n).fill(0);
    for (let i = 0; i < n - 1; i++) { h.push(T[i + 1] - T[i]); d.push((y[i + 1] - y[i]) / h[i]); }
    for (let i = 1; i < n - 1; i++) {
      if (d[i - 1] * d[i] > 0) { const w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1]; m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]); }
    }
    const end = (h0, h1, d0, d1) => {
      const v = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
      if (v * d0 <= 0) return 0;
      if (d0 * d1 <= 0 && Math.abs(v) > Math.abs(3 * d0)) return 3 * d0;
      return v;
    };
    m[0] = end(h[0], h[1], d[0], d[1]); m[n - 1] = end(h[n - 2], h[n - 3], d[n - 2], d[n - 3]);
    const upper = (arr, v) => { let lo = 0, hi = arr.length; while (lo < hi) { const mid = (lo + hi) >> 1; if (arr[mid] <= v) lo = mid + 1; else hi = mid; } return lo; };
    const lnRho = t => {
      const i = Math.min(Math.max(upper(T, t) - 1, 0), n - 2), hh = T[i + 1] - T[i], s = (t - T[i]) / hh;
      return (2 * s ** 3 - 3 * s ** 2 + 1) * y[i] + (s ** 3 - 2 * s ** 2 + s) * hh * m[i]
        + (-2 * s ** 3 + 3 * s ** 2) * y[i + 1] + (s ** 3 - s ** 2) * hh * m[i + 1];
    };
    return {
      fluid, T, rho, T_range: [T[0], T[n - 1]], rho_range: [rho[0], rho[n - 1]],
      rho_dew: t => Math.exp(lnRho(t)),
      T_dew(r) {                                   // -> [T, clamped]
        if (r <= rho[0]) return [T[0], r < rho[0]];
        if (r >= rho[n - 1]) return [T[n - 1], r > rho[n - 1]];
        const i = upper(rho, r) - 1;
        if (rho[i] === r) return [T[i], false];
        const target = Math.log(r);
        let a = T[i], b = T[i + 1];
        for (let k = 0; k < 100 && b - a >= 1e-11; k++) { const mid = 0.5 * (a + b); if (lnRho(mid) <= target) a = mid; else b = mid; }
        return [0.5 * (a + b), false];
      },
    };
  }

  // R600a density at its dew pressure for T and T + 15 K: the rated suction state (same T grid)
  const R600A_RATED_SUPERHEAT_K = 15;
  const R600A_SUCTION_15K = satVapourDensity(R600A_SAT_VAPOUR.T,
    [0.49750999, 0.52421361, 0.55204793, 0.58104671, 0.61124424, 0.64267531, 0.67537522, 0.70937977, 0.74472529, 0.78144859, 0.81958697, 0.85917827, 0.90026077, 0.9428733, 0.98705513, 1.0328461, 1.0802864, 1.1294168, 1.1802787, 1.2329136, 1.2873639, 1.3436723, 1.401882, 1.4620366, 1.5241803, 1.5883578, 1.6546143, 1.7229953, 1.7935471, 1.8663162, 1.9413498, 2.0186955, 2.0984015, 2.1805164, 2.2650893, 2.35217, 2.4418086, 2.5340559, 2.6289632, 2.7265821, 2.8269652, 2.9301652, 3.0362356, 3.1452305, 3.2572045, 3.3722128, 3.4903111, 3.6115559, 3.7360042, 3.8637136, 3.9947424, 4.1291496, 4.2669948, 4.4083383, 4.553241, 4.7017646, 4.8539717, 5.0099253, 5.1696894, 5.3333287, 5.5009086, 5.6724955, 5.8481566, 6.0279598, 6.2119741, 6.4002692, 6.5929159, 6.7899858, 6.9915515, 7.1976867, 7.4084661, 7.6239654, 7.8442615, 8.0694322, 8.2995567, 8.5347152, 8.7749893, 9.0204617, 9.2712165, 9.5273391, 9.7889162, 10.056036, 10.328788, 10.607264, 10.891556, 11.181759, 11.477968, 11.780281, 12.088798, 12.40362, 12.724851, 13.052594, 13.386958, 13.728051, 14.075984, 14.430872, 14.79283, 15.161976, 15.538431, 15.922318, 16.313761, 16.712891, 17.119837, 17.534734, 17.957719, 18.388932, 18.828515, 19.276615, 19.733382, 20.198969, 20.673533, 21.157234, 21.650237, 22.152709, 22.664824, 23.186759, 23.718693, 24.260814, 24.813311, 25.37638, 25.95022, 26.535039, 27.131046, 27.738458, 28.357499, 28.988396, 29.631385, 30.286707, 30.954611, 31.635351, 32.329191, 33.036402, 33.75726, 34.492053, 35.241076, 36.004633, 36.783037, 37.57661, 38.385686, 39.210608, 40.05173, 40.909416, 41.784045, 42.676005, 43.585698, 44.51354, 45.459958, 46.425396, 47.410312, 48.415179, 49.440486, 50.486739, 51.554461, 52.644194, 53.756497, 54.891951, 56.051153, 57.234726, 58.443311, 59.677574, 60.938203, 62.225911, 63.541437, 64.885546, 66.25903, 67.662711, 69.097441, 70.564102, 72.06361, 73.596919, 75.165016, 76.768932, 78.409742, 80.088569, 81.806592, 83.565055, 84.459862, 85.365282, 86.281496, 87.208691, 88.147064, 89.096824, 90.05819, 91.031395, 92.016693, 93.014358, 94.024694, 95.048041, 96.08479, 97.135402, 97.666078, 98.20044, 98.738587, 99.280631, 99.826706, 100.37698, 100.93168, 101.49111, 102.05575, 102.62656], "R600a at 15 K suction superheat");

  // tc* on the condensing saturated vapour; to* on the suction density (suctionRef) or the evaporating saturated vapour
  function equivalentPoint(ref, rhoEvap, rhoCond, suctionRef = null, rhoSuction = null) {
    const [eRef, rE, eName] = suctionRef ? [suctionRef, rhoSuction, "suction"] : [ref, rhoEvap, "evaporating vapour"];
    const [to, loC] = eRef.T_dew(rE), [tc, hiC] = ref.T_dew(rhoCond), warnings = [];
    for (const [name, R, clamped, r, t] of [[eName, eRef, loC, rE, to], ["condensing vapour", ref, hiC, rhoCond, tc]])
      if (clamped) warnings.push(`${name} density ${r.toPrecision(4)} kg/m3 is beyond the ${R.fluid} range: map read at ${t} C`);
    return [to, tc, warnings];
  }

  // map point for a working-fluid point: [to, tc, warnings]. map_basis "density": suction density (R600a at
  // 15 K superheat) and condensing saturated vapour; "saturated": saturated vapour both sides; "temperature": as given
  function mapPoint(comp, to, tc, rhoEvap, rhoCond, rhoSuction = null) {
    if (comp.map_basis === "density") return equivalentPoint(comp.density_reference, rhoEvap, rhoCond, comp.suction_reference, rhoSuction);
    if (comp.map_basis === "saturated") return equivalentPoint(comp.density_reference, rhoEvap, rhoCond);
    return [to, tc, []];
  }

  // ---- maps in a working fluid's (to, tc) (density_basis.fluid_curve / envelope_in_fluid, EfficiencyMap.fluid_grid)
  function fluidCurve(rhoOfT, T_lo, T_hi, step = 1, fluid = "") {
    const T = [], rho = [], n = Math.floor((T_hi - T_lo) / step + 1e-9);
    for (let k = 0; k <= n; k++) {
      const t = T_lo + k * step; let r;
      try { r = rhoOfT(t); } catch (e) { break; }
      if (!(r > 0) || !isFinite(r) || (rho.length && r <= rho[rho.length - 1])) break;
      T.push(t); rho.push(r);
    }
    if (T.length < 3) throw new Error(`could not tabulate the saturated vapour density of ${fluid || "the fluid"}`);
    return satVapourDensity(T, rho, fluid);
  }
  function toReference(curve, ref, T) {
    if (!(T >= curve.T_range[0] && T <= curve.T_range[1])) return [null, false];
    return ref.T_dew(curve.rho_dew(T));
  }
  function fromReference(curve, ref, Tref) { const [T, clamped] = curve.T_dew(ref.rho_dew(Tref)); return clamped ? null : T; }
  function envelopeInFluid(env, ref, curve, step = 1, evapRef = null, evapCurve = null) {
    const [er, ec] = evapRef ? [evapRef, evapCurve] : [ref, curve];
    const pts = [];
    for (let i = 0; i < env.length; i++) {
      const [x1, y1] = env[i], [x2, y2] = env[(i + 1) % env.length];
      const n = Math.max(1, Math.ceil(Math.max(Math.abs(x2 - x1), Math.abs(y2 - y1)) / step));
      for (let k = 0; k < n; k++) {
        const fa = fromReference(ec, er, x1 + (x2 - x1) * k / n), fb = fromReference(curve, ref, y1 + (y2 - y1) * k / n);
        if (fa !== null && fb !== null) pts.push([fa, fb, i]);      // i: the edge (EDGE_KINDS)
      }
    }
    return pts;
  }
  function fluidGrid(comp, toAxis, tcAxis, curve = null, suctionCurve = null) {
    const dens = comp.map_basis === "density" || comp.map_basis === "saturated";
    const to_map = comp.map_basis === "density" ? toAxis.map(t => toReference(suctionCurve, comp.suction_reference, t)[0])
      : dens ? toAxis.map(t => toReference(curve, comp.density_reference, t)[0]) : toAxis.slice();
    const tc_map = dens ? tcAxis.map(t => toReference(curve, comp.density_reference, t)[0]) : tcAxis.slice();
    // outside the window the cells show the extension (extendedEfficiencies), with the ratio solve() uses
    const ratio = (i, j) => comp.map_basis === "density" ? extensionRatio(comp, to_map[i], tc_map[j], curve.rho_dew(tcAxis[j]), suctionCurve.rho_dew(toAxis[i]))
      : extensionRatio(comp, to_map[i], tc_map[j]);
    const ext = tc_map.map((b, j) => to_map.map((a, i) => (a === null || b === null ? null : extendedEfficiencies(comp, a, b, ratio(i, j)))));
    return {
      to: toAxis.slice(), tc: tcAxis.slice(), to_map, tc_map,
      eta_is: ext.map(row => row.map(e => (e ? e.eta_is : null))),
      eta_vol: ext.map(row => row.map(e => (e ? e.eta_vol : null))),
      in_env: tc_map.map(b => to_map.map(a => (a === null || b === null ? null : inEnvelope(a, b)))),
    };
  }

  // ---------------------------------------------------------------- beyond the application window
  /* The window's edges are of two kinds (compressor_data.py, same numbers):
     ratio edges, where efficiency falls away (high: density ratio ~19, volumetric efficiency ~0.5; low: pressure
     ratio ~1.5-1.7), and mechanical edges (seals, materials, motor), which say nothing about efficiency.
     Outside the window the maps are read at the nearest window point (continuous across the edge), then:
       above the high-ratio line (density ratio r over its limit r_hi): clearance re-expansion,
         eta_vol = 1 - (1 - eta_vol,q)(r - 1)/(r_hi - 1), and eta_is = eta_is,q (eta_vol/eta_vol,q)^beta
         (both reach zero together; beta 0.4 gives about 25 % lower eta_is at 1.5 r_hi);
       below the low-ratio line: valve losses take over as the lift vanishes, with x = r - 1,
         eta_is = eta_is,q [x/(x + delta)] / [x_lo/(x_lo + delta)], delta 0.3; eta_vol held (at most 1).
     The ratio lines follow the ratio edges (piecewise linear in to, flat beyond their ends), so the falloff starts
     exactly at the window. The mechanical limits are checked on the real cycle (limitFlags). */
  const EDGE_KINDS = ["mech", "ratio_low", "mech", "mech", "ratio_high", "mech"];     // ENVELOPE edge i -> i + 1
  const EDGE_NAMES = ["minimum condensing", "low-ratio", "maximum evaporating", "maximum condensing",
    "high-ratio", "minimum evaporating"];
  const EXTENSION = { beta: 0.4, delta_low: 0.3, eta_vol_floor: 0.02, eta_is_floor: 0.05 };
  // mechanical limits, from R600a over the window at 15 K suction superheat with no shell loss (compressor_data.py):
  // high-side pressure (dew at 120 C), suction pressure (dew at -10 C), discharge temperature and electrical work
  // per swept volume (the shaft-torque load), each the window's maximum rounded up
  const LIMITS = {
    "Recip-S": { P_dis_max_bar: 28.4, P_suc_min_bar: 1.08, T_dis_max: 140, w_swept_max_kJm3: 1085 },
    "Recip-M": { P_dis_max_bar: 28.4, P_suc_min_bar: 1.08, T_dis_max: 135, w_swept_max_kJm3: 1015 },
    "Recip-L": { P_dis_max_bar: 28.4, P_suc_min_bar: 1.08, T_dis_max: 137, w_swept_max_kJm3: 1015 },
  };
  // density ratio at a map point (R600a): condensing saturated vapour over the 15 K-superheat suction
  // density ratio for the extension (EfficiencyMap.extension_ratio): the fluid's own on the density basis, else
  // the map point's R600a ratio, so the falloff starts at the window on every basis
  function extensionRatio(comp, to, tc, rhoCond = null, rhoSuction = null) {
    return comp.map_basis === "density" && rhoCond && rhoSuction ? rhoCond / rhoSuction : mapDensityRatio(to, tc);
  }
  function mapDensityRatio(to, tc) {
    const cl = (t, R) => Math.min(Math.max(t, R.T_range[0]), R.T_range[1]);
    return R600A_SAT_VAPOUR.rho_dew(cl(tc, R600A_SAT_VAPOUR)) / R600A_SUCTION_15K.rho_dew(cl(to, R600A_SUCTION_15K));
  }
  // [to, density ratio] along the edges of one kind, by to
  function ratioLine(kind, env = ENVELOPE, kinds = EDGE_KINDS, n = 16) {
    const pts = [];
    env.forEach((a, i) => {
      if (kinds[i] !== kind) return;
      const b = env[(i + 1) % env.length];
      for (let k = 0; k <= n; k++) { const f = k / n, to = a[0] + f * (b[0] - a[0]), tc = a[1] + f * (b[1] - a[1]); pts.push([to, mapDensityRatio(to, tc)]); }
    });
    return pts.sort((p, q) => p[0] - q[0]);
  }
  function lineAt(line, x) {
    if (x <= line[0][0]) return line[0][1];
    if (x >= line[line.length - 1][0]) return line[line.length - 1][1];
    let i = 0; while (line[i + 1][0] < x) i++;
    const [x0, y0] = line[i], [x1, y1] = line[i + 1];
    return x1 === x0 ? y1 : y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  }
  const RATIO_LINES = { high: ratioLine("ratio_high"), low: ratioLine("ratio_low") };
  // nearest point of the window's boundary to (to, tc), and the edge it lies on
  function nearestOnEnvelope(to, tc, env = ENVELOPE) {
    let best = null;
    env.forEach((a, i) => {
      const b = env[(i + 1) % env.length], dx = b[0] - a[0], dy = b[1] - a[1];
      const t = Math.min(Math.max(((to - a[0]) * dx + (tc - a[1]) * dy) / (dx * dx + dy * dy), 0), 1);
      const x = a[0] + t * dx, y = a[1] + t * dy, d = Math.hypot(to - x, tc - y);
      if (!best || d < best.d - 1e-12) best = { to: x, tc: y, d, edge: i, t };
    });
    return best;
  }
  /* efficiencies at map point (to, tc), with the working fluid's density ratio r (condensing saturated vapour over
     suction): the map inside the window; outside, the extension above */
  function extendedEfficiencies(comp, to, tc, r) {
    const env = comp.envelope || ENVELOPE, kinds = comp.edge_kinds || EDGE_KINDS;
    if (inEnvelope(to, tc, env)) return { inside: true, eta_is: comp.eta_is(to, tc), eta_vol: comp.eta_vol(to, tc), read: [to, tc], ratio: r, ratio_beyond: 0 };
    const q = nearestOnEnvelope(to, tc, env), X = EXTENSION;
    let eis = comp.eta_is(q.to, q.tc), ev = comp.eta_vol(q.to, q.tc);
    const out = { inside: false, read: [q.to, q.tc], edge: q.edge, edge_kind: kinds[q.edge], edge_name: EDGE_NAMES[q.edge] || `edge ${q.edge + 1}`,
      eta_is_edge: eis, eta_vol_edge: ev, ratio: r, ratio_beyond: 0 };
    const rHi = lineAt(RATIO_LINES.high, to), rLo = lineAt(RATIO_LINES.low, to);
    if (r > rHi) {
      const evq = ev;
      ev = 1 - (1 - evq) * (r - 1) / (rHi - 1);
      ev = Math.max(ev, X.eta_vol_floor);
      eis = eis * Math.pow(ev / evq, X.beta);
      Object.assign(out, { ratio_limit: rHi, ratio_beyond: r / rHi - 1 });
    } else if (r < rLo) {
      const x = Math.max(r - 1, 0), xl = rLo - 1, d = X.delta_low;
      eis = eis * (x / (x + d)) / (xl / (xl + d));
      ev = Math.min(ev, 1);
      Object.assign(out, { ratio_limit: rLo, ratio_beyond: r / rLo - 1 });
    }
    out.eta_is = Math.max(eis, X.eta_is_floor); out.eta_vol = ev;
    return out;
  }
  // mechanical limits on the real cycle -> [{key, label, value, limit, unit, over}] for those exceeded
  function limitFlags(limits, v) {
    if (!limits) return [];
    const chk = [["P_dis", "High-side pressure", v.P_dis_bar, limits.P_dis_max_bar, "bar", 1],
      ["T_dis", "Discharge temperature", v.T_dis, limits.T_dis_max, "°C", 1],
      ["w_swept", "Work per swept volume (shaft load)", v.w_swept_kJm3, limits.w_swept_max_kJm3, "kJ/m³", 1],
      ["P_suc", "Suction pressure", v.P_suc_bar, limits.P_suc_min_bar, "bar", -1]];
    return chk.filter(([, , val, lim, , sgn]) => lim != null && isFinite(val) && sgn * (val - lim) > 0)
      .map(([key, label, value, limit, unit, sgn]) => {
        const f = { key, label, value, limit, unit, kind: sgn > 0 ? "max" : "min", over_pct: 100 * Math.abs(value / limit - 1) };
        f.detail = limitDetail(f, v);
        return f;
      });
  }
  // a sentence for each exceeded limit: what it is, by how much, and why it matters
  function limitDetail(f, v) {
    const d = f.unit === "bar" ? 2 : 0, n = x => x.toFixed(d), pct = `${f.over_pct.toFixed(1)} %`;
    if (f.key === "P_dis") return `High-side pressure ${n(f.value)} bar is ${pct} beyond the operating window (${n(f.limit)} bar): check with the supplier that the casing, gaskets and valves are rated for it.`;
    if (f.key === "T_dis") return `Discharge temperature ${n(f.value)} °C is ${(f.value - f.limit).toFixed(0)} K beyond the operating window (${n(f.limit)} °C): check with the supplier that the oil and discharge valves are suitable.`;
    if (f.key === "w_swept") {
      const what = v.rev_m3 ? `Shaft torque ≈ ${(f.value * (v.shaft_fraction ?? 1) * v.rev_m3 / (2 * Math.PI) * 1000).toFixed(0)} N·m (${n(f.value)} kJ/m³ swept)` : `Shaft load ${n(f.value)} kJ/m³ swept`;
      return `${what} is ${pct} beyond the operating window (${n(f.limit)} kJ/m³): check with the supplier that the motor, crank and bearings are suitable.`;
    }
    if (f.key === "P_suc") return `Suction pressure ${n(f.value)} bar is ${pct} below the operating window (${n(f.limit)} bar): check with the supplier for air ingress and motor cooling.`;
    return `${f.label} ${n(f.value)} ${f.unit} is beyond the ${n(f.limit)} ${f.unit} limit.`;
  }

  function compressorMap(model = "Recip-M", opts = {}) {
    const c = COEFFS[model];
    if (!c) throw new Error(`unknown compressor ${model}`);
    return Object.assign({
      name: model, displacement_m3h: DISPLACEMENT_50HZ[model], rated_Hz: 50, rated_rpm: 1450, frequency_Hz: 50,
      basis: "electrical", motor_eff: 0.93, eta_is_bounds: [0.2, 1.0], eta_vol_bounds: [0.0, 1.1],
      envelope: ENVELOPE, edge_kinds: EDGE_KINDS, limits: LIMITS[model] || null,
      eta_is: (to, tc) => bicubic(c.is, to, tc), eta_vol: (to, tc) => bicubic(c.vol, to, tc),
      map_basis: "density", density_reference: R600A_SAT_VAPOUR, suction_reference: R600A_SUCTION_15K,
    }, opts);
  }

  function clip(v, [lo, hi], name, warnings) {
    if (v < lo || v > hi) { warnings.push(`${name} ${v.toFixed(3)} clipped to ${lo}-${hi}`); return Math.min(Math.max(v, lo), hi); }
    return v;
  }

  // ---------------------------------------------------------------- numerics
  function brent(f, a, b, tol = 1e-10, maxit = 100) {
    let fa = f(a), fb = f(b);
    if (fa * fb > 0) throw new Error("root not bracketed");
    if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }
    let c = a, fc = fa, d = b - a, mflag = true;
    for (let it = 0; it < maxit; it++) {
      if (fb === 0 || Math.abs(b - a) < tol) return b;
      let s;
      if (fa !== fc && fb !== fc) s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
      else s = b - fb * (b - a) / (fb - fa);
      const cond = (s - (3 * a + b) / 4) * (s - b) >= 0 || (mflag && Math.abs(s - b) >= Math.abs(b - c) / 2) ||
        (!mflag && Math.abs(s - b) >= Math.abs(c - d) / 2) || (mflag && Math.abs(b - c) < tol) || (!mflag && Math.abs(c - d) < tol);
      if (cond) { s = (a + b) / 2; mflag = true; } else mflag = false;
      const fs = f(s); d = c; c = b; fc = fb;
      if (fa * fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
      if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }
    }
    return b;
  }

  // ---------------------------------------------------------------- configuration
  const DEFAULTS = {
    T_source: 55, dTmin_source: 5, T_steam: 120, dTmin_steam: 5, T_water_in: 80, T_water_out: 95, dTmin_water: 5,
    T_evap: null, T_cond: null, dT_subcool: null, dT_superheat: 15, evap_basis: "dew", cond_basis: "dew",
    dTmin_ihx: 5,
    compressor_Hz: 50, T_ambient: 15, loss_dT_basis: "mean", loss_UA_kWK: null,
    loss_UA_cond_kWK: 0, loss_UA_subcool_kWK: 0, loss_UA_evap_kWK: 0, preheat_subcool_share: 0, dP_K_cond: 0, dP_K_subcool: 0, dP_K_evap: 0,     // exchanger heat losses (as the Python model)
    dT_discharge_superheat: null, subcooler_heat_useful: true, dP_crit_min: 0,
  };

  // ---------------------------------------------------------------- hot-water subcooler
  // (ConventionalHeatPump._water_subcooler) outlet as cold as dTmin_water allows; water linear in duty
  const WATER_CP = 4.19;
  function waterSubcooler(pr, c, P_H, scIn, n = 20) {
    const T_in = scIn.T, h_in = scIn.h, dT = c.dTmin_water;
    if (T_in - c.T_water_out < dT - 1e-9)
      return { message: `hot water outlet ${c.T_water_out} C is less than ${dT} K below the refrigerant entering the subcooler (${T_in.toFixed(1)} C)` };
    const approach = T_out => {          // n samples, then golden section around the smallest interior one
      const h_out = pr.PT(P_H, T_out).h;
      const gap = x => pr.PH(P_H, h_out + (h_in - h_out) * x).T - (c.T_water_in + (c.T_water_out - c.T_water_in) * x);
      const ends = Math.min(T_out - c.T_water_in, T_in - c.T_water_out);
      let g = Infinity, i = 0;
      for (let k = 1; k < n; k++) { const v = gap(k / n); if (v < g) { g = v; i = k; } }
      if (g >= ends) return ends;
      let a = (i - 1) / n, b = (i + 1) / n; const r = (Math.sqrt(5) - 1) / 2;
      for (let k = 0; k < 40; k++) { const x1 = b - r * (b - a), x2 = a + r * (b - a); if (gap(x1) < gap(x2)) b = x2; else a = x1; }
      return Math.min(g, gap(0.5 * (a + b)), ends);
    };
    let T_out = c.T_water_in + dT;
    if (T_out >= T_in) return { out: scIn, min_approach: T_in - c.T_water_in };
    if (approach(T_out) < dT - 1e-6) T_out = brent(T => approach(T) - dT, T_out, T_in - 1e-6, 1e-6);
    const out = pr.PT(P_H, T_out);
    if (!liquidLike(out)) return { message: `subcooler outlet not liquid (${out.phase})` };
    return { out, min_approach: approach(T_out) };
  }

  // ---------------------------------------------------------------- saturation
  function saturationPressure(pr, T, basis) {
    if (basis === "dew") return pr.TQ(T, 1).P;
    if (basis === "bubble") return pr.TQ(T, 0).P;
    if (basis !== "mean") throw new Error(`unknown saturation basis ${basis}`);
    const Pd = pr.TQ(T, 1).P, Pb = pr.TQ(T, 0).P;
    if (Math.abs(Pb / Pd - 1) < 1e-9) return Pd;
    const f = lnP => { const P = Math.exp(lnP); return 0.5 * (pr.PQ(P, 1).T + pr.PQ(P, 0).T) - T; };
    return Math.exp(brent(f, Math.log(Pd), Math.log(Pb), 1e-10));
  }

  // liquid, or the dense fluid above the critical pressure (supercritical high side)
  const liquidLike = st => st.phase === "liquid" || st.phase === "saturated_liquid" || st.phase === "supercritical";

  /* Supercritical high side. Above the critical temperature there is no bubble or dew point; the high-side
     "saturation" continues along the critical isochore (rho = rho_crit), which meets the saturation line at the
     critical point. critIsochore(P): the state on that isochore at pressure P (P >= P_crit). */
  function critIsochore(pr, P) {
    const rc = pr.rho_crit;
    if (!(rc > 0)) throw new Error("critical density not available for this refrigerant");
    const g = T => pr.TD(T, rc).P - P;
    let lo = pr.T_crit - 0.5, hi = pr.T_crit + 20;
    for (let k = 0; k < 12 && g(hi) < 0; k++) hi += 40;
    const T = brent(g, lo, hi, 1e-9);
    const st = Object.assign({}, pr.TD(T, rc)); st.P = P; st.phase = "supercritical"; st.q = -1;
    return st;
  }

  const wallT = (cfg, T1, T2) => cfg.loss_dT_basis === "discharge" ? T2 : 0.5 * (T1 + T2);

  // ---------------------------------------------------------------- one operating point
  /* Refrigerant pressure drops: dP [kPa] = K (V [m3/h])^1.8 in the condenser, subcooler and evaporator, with V
     the volume flow at the mean of the inlet and outlet specific volumes. The exits stay pinned by the approaches
     (condenser exit bubble point, evaporator exit dew point), so the compressor discharges at P_H + dP_cond and the
     valve feeds the evaporator at P_L + dP_evap; the drops depend on the flow, so the point is re-solved until
     they settle. */
  function solve(pr, cfgIn, comp, UA = null) {
    const c = Object.assign({}, DEFAULTS, cfgIn);
    const K = { cond: c.dP_K_cond || 0, sub: c.dP_K_subcool || 0, evap: c.dP_K_evap || 0 };
    let dP = { cond: 0, sub: 0, evap: 0 };
    let res = solveOnce(pr, c, comp, UA, dP);
    if (!(K.cond || K.sub || K.evap) || !res.success) return res;
    for (let it = 0; it < 15; it++) {
      const r = res.report;
      const next = { cond: K.cond * Math.pow(r.Vmean_cond_m3h, 1.8) / 1000, sub: K.sub * Math.pow(r.Vmean_subcool_m3h, 1.8) / 1000,
        evap: K.evap * Math.pow(r.Vmean_evap_m3h, 1.8) / 1000 };            // MPa
      const err = Math.abs(next.cond - dP.cond) + Math.abs(next.sub - dP.sub) + Math.abs(next.evap - dP.evap);
      dP = next;
      res = solveOnce(pr, c, comp, UA, dP);
      if (!res.success || err < 1e-7) break;
    }
    return res;
  }
  function solveOnce(pr, c, comp, UA, dP) {
    const warnings = [];
    const fail = (message, extra = {}) => Object.assign({ success: false, message, warnings }, extra);
    try {
      if ((c.dT_subcool != null && c.dT_subcool < 0) || c.dT_superheat < 0) return fail("subcooling and superheat must be >= 0");
      if (c.dT_subcool == null && !(c.T_water_out > c.T_water_in)) return fail("hot water outlet must be above its inlet");
      // evaporator outlet (dew point) dTmin_source below the heat source, unless T_evap is given directly
      const T_hot = c.T_cond != null ? c.T_cond : c.T_steam + c.dTmin_steam;
      // at or above the critical temperature the high side is supercritical: the gas-cooler exit (where the bubble
      // point would be) at T_hot on the critical isochore sets the pressure
      const sc = isFinite(pr.T_crit) && T_hot >= pr.T_crit;
      if (sc && !(pr.rho_crit > 0))
        return fail(`condensing at ${T_hot.toFixed(1)} C is at or above the critical temperature of this refrigerant (${pr.T_crit.toFixed(1)} C), and its critical density is not available`);
      const hiSat = (P, q) => (sc ? critIsochore(pr, P) : pr.PQ(P, q));
      const P_L = c.T_evap != null ? saturationPressure(pr, c.T_evap, c.evap_basis) : saturationPressure(pr, c.T_source - c.dTmin_source, "dew");
      // condenser outlet (bubble point) dTmin_steam above the steam, unless T_cond is given directly
      const gcOut = sc ? pr.TD(T_hot, pr.rho_crit) : null;
      const P_H = sc ? gcOut.P : c.T_cond != null ? saturationPressure(pr, c.T_cond, c.cond_basis) : saturationPressure(pr, c.T_steam + c.dTmin_steam, "bubble");
      if (!(P_H > P_L)) return fail("condensing pressure not above evaporating pressure", { P_H, P_L });
      if (!sc && pr.P_crit && isFinite(pr.P_crit) && P_H > pr.P_crit - c.dP_crit_min)
        return fail(`condensing pressure ${(P_H * 10).toFixed(2)} bar is above the critical pressure ${(pr.P_crit * 10).toFixed(2)} bar`, { P_H, P_L });
      // pressure drops (exits pinned): compressor discharge and evaporator inlet pressures
      const P_dis = P_H + dP.cond, P_evin = P_L + dP.evap, P_valve = P_H - dP.sub;
      if (!sc && pr.P_crit && isFinite(pr.P_crit) && P_dis > pr.P_crit - c.dP_crit_min)
        return fail(`discharge pressure ${(P_dis * 10).toFixed(2)} bar (with the condenser pressure drop) is above the critical pressure ${(pr.P_crit * 10).toFixed(2)} bar`, { P_H, P_L });
      if (!(P_valve > P_evin)) return fail("pressure drops leave no pressure difference across the valve", { P_H, P_L });
      // condDew is the dew point at the condenser inlet (compressor discharge pressure); supercritical: the
      // critical-isochore states at the discharge pressure and at the gas-cooler exit (T_hot)
      const evapDew = pr.PQ(P_L, 1), evapBub = pr.PQ(P_L, 0), condDew = hiSat(P_dis, 1);
      const condBub = sc ? Object.assign({}, gcOut, { P: P_H, phase: "supercritical", q: -1 }) : pr.PQ(P_H, 0);
      let condOut = condBub;                 // 3 condenser exit (subcooled below, after the compressor)
      const share = Math.min(Math.max(c.preheat_subcool_share || 0, 0), 1);
      // steam from bled hot water: needs the hot-water subcooler (dT_subcool null), as the steam bookkeeping does
      const steamMode = c.T_cond == null && c.dT_subcool == null && c.T_water_out != null && c.T_water_out < c.T_steam;

      const suction = c.dT_superheat > 0 ? pr.PT(P_L, evapDew.T + c.dT_superheat) : evapDew;
      const evapOut = evapDew;              // saturated: the IHX gives all the suction superheat
      const q_ihx = suction.h - evapOut.h, ihx = c.dT_superheat > 0;

      // compressor
      const isen = pr.PS(P_dis, suction.s), h_is = isen.h;
      const [to, tc, mapWarn] = mapPoint(comp, evapDew.T, condDew.T, evapDew.D, condDew.D, suction.D);
      warnings.push(...mapWarn);
      // condensing saturated vapour (critical isochore when supercritical) over the suction
      const ext = extendedEfficiencies(comp, to, tc, extensionRatio(comp, to, tc, condDew.D, suction.D));
      if (!ext.inside) {
        const r = ext.ratio_beyond > 0 ? `; density ratio ${ext.ratio.toFixed(2)}, ${(100 * ext.ratio_beyond).toFixed(0)} % above the high-ratio limit: efficiency falls off`
          : ext.ratio_beyond < 0 ? `; density ratio ${ext.ratio.toFixed(2)}, ${(-100 * ext.ratio_beyond).toFixed(0)} % below the low-ratio limit: efficiency falls off` : "";
        warnings.push(`outside the ${comp.name} envelope (${comp.map_basis !== "temperature" ? "R600a-equivalent " : ""}to ${to.toFixed(1)} C, tc ${tc.toFixed(1)} C): read at the nearest window point, beyond the ${ext.edge_name} edge${r}`);
      }
      // beyond the window the extension's falloff (down to its own floor) replaces the lower guard rail
      const eta = clip(ext.eta_is, ext.inside ? comp.eta_is_bounds : [Math.min(comp.eta_is_bounds[0], EXTENSION.eta_is_floor), comp.eta_is_bounds[1]], "eta_is", warnings);
      const eta_vol = clip(ext.eta_vol, comp.eta_vol_bounds, "eta_vol", warnings);
      const dh_is = h_is - suction.h;
      const w_el = comp.basis === "electrical" ? dh_is / eta : dh_is / eta / comp.motor_eff;
      const w_sh = comp.basis === "electrical" ? w_el * comp.motor_eff : dh_is / eta;
      const Hz = c.compressor_Hz ?? comp.frequency_Hz;
      const swept = comp.displacement_m3h * Hz / comp.rated_Hz;
      const m = eta_vol * swept / 3600 * suction.D;
      const rpm = comp.rated_rpm * Hz / comp.rated_Hz;

      let discharge, q_loss, mode;
      if (c.dT_discharge_superheat != null) {
        mode = "discharge_superheat";
        discharge = pr.PT(P_dis, condDew.T + c.dT_discharge_superheat);
        q_loss = w_el - (discharge.h - suction.h);
        if (q_loss < 0) return fail(`discharge superheat ${c.dT_discharge_superheat} K needs ${(-q_loss).toFixed(1)} kJ/kg more than the compressor's electrical input`, { P_H, P_L });
      } else {
        mode = "UA";
        const ua = UA ?? c.loss_UA_kWK ?? 0;
        const h_ad = suction.h + w_el;
        if (ua === 0) { discharge = pr.PH(P_dis, h_ad); q_loss = 0; }
        else {
          const g = h2 => (h2 - suction.h) - w_el + ua * (wallT(c, suction.T, pr.PH(P_dis, h2).T) - c.T_ambient) / m;
          const g0 = g(h_ad);
          let a = h_ad, b = h_ad, step = Math.abs(g0) + 1, ok = false;
          for (let i = 0; i < 60; i++) {
            if (g0 > 0) { a = b - step; if (g(a) <= 0) { ok = true; break; } b = a; }
            else { b = a + step; if (g(b) >= 0) { ok = true; break; } a = b; }
            step *= 2;
          }
          if (!ok) return fail("compressor heat-loss balance not bracketed", { P_H, P_L });
          const h2 = brent(g, a, b, 1e-9);
          discharge = pr.PH(P_dis, h2);
          q_loss = w_el - (discharge.h - suction.h);
        }
        UA = ua;
      }
      if (discharge.phase !== "superheated_vapour" && discharge.phase !== "supercritical" && discharge.phase !== "saturated_vapour")
        return fail(`wet compression (discharge ${discharge.phase.replace("_", " ")})`, { P_H, P_L });

      // liquid line: subcooler (hot water), IHX
      // 3 condenser (steam generator) exit. With preheat_subcool_share > 0 the liquid is subcooled so that it
      // supplies that share of the feedwater preheat: per kg of refrigerant, q_x = share fp (q_c0 + q_x - l_c),
      // with fp = preheat / (steam - feed) enthalpy rise of the water and l_c the condenser loss per kg
      if (share > 0 && steamMode) {
        const s1 = steamFromFeedwater(1, c.T_steam, c.T_water_out), fp = s1.preheat_kJkg / s1.dh_kJkg;
        const q_c0 = discharge.h - condBub.h;
        for (let k = 0; k < 4; k++) {         // the condenser loss depends a little on the exit temperature
          const l_c = (c.loss_UA_cond_kWK || 0) * ((discharge.T + condOut.T) / 2 - c.T_ambient) / m;
          const q_x = share * fp * (q_c0 - l_c) / (1 - share * fp);
          condOut = pr.PH(P_H, condBub.h - q_x);
        }
        if (!liquidLike(condOut)) return fail(`condenser exit not liquid (${condOut.phase})`, { P_H, P_L });
      } else if (share > 0) warnings.push("preheat from subcooling needs the steam and hot-water specification; ignored");
      let scOut = condOut, water = null;
      if (c.dT_subcool != null) {
        if (c.dT_subcool > 0) {
          scOut = pr.PT(P_H, condOut.T - c.dT_subcool);
          if (!liquidLike(scOut)) return fail("subcooler outlet not liquid", { P_H, P_L });
        }
      } else {
        const ws = waterSubcooler(pr, c, P_H, condOut);
        if (ws.message) return fail(ws.message, { P_H, P_L });
        scOut = ws.out; water = ws;
      }
      const h1 = suction.h, h2 = discharge.h, h3 = condOut.h, h4 = scOut.h;
      // IHX after the subcooler; its liquid outlet is the valve inlet
      const ihxIn = scOut, hIhxOut = h4 - q_ihx, hValve = hIhxOut;
      const ihxOut = q_ihx > 0 ? pr.PH(P_H, hIhxOut) : ihxIn;
      let ihxApproach = null;
      if (ihx) {
        ihxApproach = Infinity;
        for (let i = 0; i <= 6; i++) {
          const d = pr.PH(P_H, hIhxOut + q_ihx * i / 6).T - pr.PH(P_L, evapOut.h + q_ihx * i / 6).T;
          ihxApproach = Math.min(ihxApproach, d);
        }
        if (ihxApproach < c.dTmin_ihx - 1e-6)
          return fail(`IHX approach ${ihxApproach.toFixed(1)} K below ${c.dTmin_ihx} K: the liquid (${ihxIn.T.toFixed(1)} C) is too cold to superheat the suction gas to ${suction.T.toFixed(1)} C`, { P_H, P_L });
      }
      const evapIn = pr.PH(P_evin, hValve);

      const q_desup = h2 - condDew.h, q_cond = h2 - h3;
      const q_sc = h3 - h4;
      const q_evap = evapOut.h - hValve;
      const q_heat = q_cond + (c.subcooler_heat_useful ? q_sc : 0);
      // exchanger heat losses to ambient (kW), from the mean refrigerant temperature in each exchanger:
      // off the heat delivered (condenser, subcooler); the source supplies the evaporator duty plus its loss
      const L_cond = (c.loss_UA_cond_kWK || 0) * ((discharge.T + condOut.T) / 2 - c.T_ambient);
      const L_sc = (c.loss_UA_subcool_kWK || 0) * ((condOut.T + scOut.T) / 2 - c.T_ambient);
      const L_evap = (c.loss_UA_evap_kWK || 0) * ((evapIn.T + evapOut.T) / 2 - c.T_ambient);
      const L_heat = L_cond + (c.subcooler_heat_useful ? L_sc : 0);
      const cop = (m * q_heat - L_heat) / (m * w_el);
      const T_hot2 = c.T_cond != null ? c.T_cond : condBub.T;     // Carnot between the condensing ...
      const T_low = c.T_evap != null ? c.T_evap : evapDew.T;      // ... and evaporating temperatures
      const r = {
        COP: cop, COP_cooling: q_evap / w_el, COP_carnot: T_hot2 > T_low ? (T_hot2 + 273.15) / (T_hot2 - T_low) : null,
        T_source: c.T_evap != null ? null : c.T_source, T_evap: T_low,
        T_steam: c.T_cond != null ? null : c.T_steam, T_cond: T_hot2, dT_subcool: condOut.T - scOut.T,
        T_water_in: water ? c.T_water_in : null, T_water_out: water ? c.T_water_out : null,
        subcooler_min_approach_K: water ? water.min_approach : null,
        VHC: (m * q_heat - L_heat) / m * suction.D / 1000, rho_suction_kgm3: suction.D, P_evap_bar: P_L * 10, P_cond_bar: P_H * 10, PR: P_dis / P_L,
        supercritical: sc, rho_crit_kgm3: sc ? pr.rho_crit : null,
        P_discharge_bar: P_dis * 10, P_evap_in_bar: P_evin * 10, P_valve_in_bar: P_valve * 10,
        dP_cond_kPa: dP.cond * 1000, dP_subcool_kPa: dP.sub * 1000, dP_evap_kPa: dP.evap * 1000,
        Vmean_cond_m3h: m * (1 / discharge.D + 1 / condOut.D) / 2 * 3600, Vmean_subcool_m3h: m * (1 / condOut.D + 1 / scOut.D) / 2 * 3600,
        Vmean_evap_m3h: m * (1 / evapIn.D + 1 / evapOut.D) / 2 * 3600,
        T_evap_dew: evapDew.T, T_evap_bubble: evapBub.T, evap_glide_K: evapDew.T - evapBub.T,
        T_cond_dew: condDew.T, T_cond_bubble: condBub.T, cond_glide_K: condDew.T - condBub.T,
        T_cond_out: condOut.T, dT_subcool_cond: condBub.T - condOut.T, preheat_subcool_share: steamMode ? share : 0,
        T_evap_in: evapIn.T, x_evap_in: evapIn.q, T_evap_out: evapOut.T, T_suction: suction.T,
        T_discharge: discharge.T, dT_discharge_superheat: discharge.T - condDew.T, T_subcool_out: scOut.T,
        T_valve_in: ihxOut.T, eta_is_overall: dh_is / w_el, eta_is_shaft: dh_is / w_sh, eta_vol,
        map_to: to, map_tc: tc,
        compressor_Hz: Hz, compressor_rpm: rpm, V_suction_m3h: m / suction.D * 3600, V_discharge_m3h: m / discharge.D * 3600,
        q_heat_kJkg: q_heat, q_cond_kJkg: q_cond, q_desuperheat_kJkg: q_desup, q_subcool_kJkg: q_sc, q_evap_kJkg: q_evap,
        q_ihx_kJkg: q_ihx, w_electric_kJkg: w_el, w_to_gas_kJkg: h2 - h1, q_compressor_loss_kJkg: q_loss,
        compressor_loss_fraction: q_loss / w_el, loss_mode: mode, UA_loss_kWK: UA, loss_dT_K: wallT(c, suction.T, discharge.T) - c.T_ambient,
        ihx, ihx_min_approach_K: ihxApproach,
        energy_balance_residual_kJkg: q_evap + (h2 - h1) - (q_cond + q_sc),
        m_dot_kgs: m, Q_heat_kW: m * q_heat - L_heat, Q_cond_kW: m * q_cond - L_cond, Q_subcool_kW: m * q_sc - L_sc, Q_evap_kW: m * q_evap,
        Q_cond_refrigerant_kW: m * q_cond, Q_subcool_refrigerant_kW: m * q_sc, Q_source_kW: m * q_evap + L_evap,
        Q_loss_cond_kW: L_cond, Q_loss_subcool_kW: L_sc, Q_loss_evap_kW: L_evap,
        UA_loss_cond_kWK: c.loss_UA_cond_kWK || 0, UA_loss_subcool_kWK: c.loss_UA_subcool_kWK || 0, UA_loss_evap_kWK: c.loss_UA_evap_kWK || 0,
        Q_ihx_kW: m * q_ihx, W_electric_kW: m * w_el, W_shaft_kW: m * w_sh, Q_compressor_loss_kW: m * q_loss,
        m_water_kgs: water ? (m * q_sc - L_sc) / (WATER_CP * (c.T_water_out - c.T_water_in)) : null,
        in_envelope: inEnvelope(to, tc),
        // beyond the window: where the maps were read, the density ratio and how far past a ratio limit it is
        window_edge: ext.inside ? null : ext.edge_name, window_edge_kind: ext.inside ? null : ext.edge_kind,
        to_read: ext.read[0], tc_read: ext.read[1], density_ratio: ext.ratio, density_ratio_beyond: ext.inside ? 0 : ext.ratio_beyond,
        w_swept_kJm3: m * w_el / (swept / 3600),
      };
      // mechanical limits (seals, materials, motor) on the real cycle: flagged, not enforced
      r.limit_flags = limitFlags(comp.limits, { P_dis_bar: P_dis * 10, T_dis: discharge.T, w_swept_kJm3: r.w_swept_kJm3, P_suc_bar: P_L * 10,
        rev_m3: comp.displacement_m3h / (comp.rated_rpm * 60), shaft_fraction: w_sh / w_el });
      // steam from hot water bled at the subcooler outlet: the steam generator (desuperheating and condensing)
      // heats it from T_water_out to saturated steam at T_steam
      r.m_steam_kgs = r.steam_feed_T = r.steam_bleed_fraction = r.m_hot_water_net_kgs = null;
      if (water && c.T_steam != null && c.T_water_out < c.T_steam) {
        const sf = steamFromFeedwater(r.Q_cond_kW, c.T_steam, c.T_water_out);
        Object.assign(r, { m_steam_kgs: sf.m_kgs, steam_feed_T: c.T_water_out, steam_dh_kJkg: sf.dh_kJkg, P_steam_bar: sf.P_steam_MPa * 10,
          steam_preheat_fraction: sf.preheat_kJkg / sf.dh_kJkg,
          steam_bleed_fraction: sf.m_kgs / r.m_water_kgs, m_hot_water_net_kgs: r.m_water_kgs - sf.m_kgs });
        // steam generator approach, counter-current: the feedwater is heated to boiling over the first
        // preheat fraction of the duty from the refrigerant's cold (exit) end, then boils at T_steam. With a
        // subcooled exit the liquid must not fall within dTmin_steam of the water anywhere.
        r.steam_gen_min_approach_K = condBub.T - c.T_steam;
        if (h3 < condBub.h - 1e-9) {
          const fp = r.steam_preheat_fraction, n = 40, fr = [];
          for (let i = 0; i <= n; i++) fr.push(i / n);
          fr.push(fp); fr.sort((a, b) => a - b);
          let mn = Infinity;
          for (const f of fr) {
            const Tr = pr.PH(P_H + dP.cond * f, h3 + (h2 - h3) * f).T;
            const Tw = f < fp ? c.T_water_out + (c.T_steam - c.T_water_out) * f / fp : c.T_steam;
            mn = Math.min(mn, Tr - Tw);
          }
          r.steam_gen_min_approach_K = mn;
          if (mn < c.dTmin_steam - 1e-6)
            return fail(`preheat from subcooling leaves the steam generator approach at ${mn.toFixed(1)} K, below ${c.dTmin_steam} K`, { P_H, P_L });
        }
      }
      if (r.m_water_kgs != null && r.m_water_kgs <= 0)
        return fail(`the subcooler delivers no hot water (${r.Q_subcool_kW.toFixed(1)} kW net of its heat loss)`, { P_H, P_L });
      if (r.m_steam_kgs != null && r.m_steam_kgs > r.m_water_kgs)
        return fail(`steam ${r.m_steam_kgs.toFixed(3)} kg/s needs more feedwater than the subcooler heats (${r.m_water_kgs.toFixed(3)} kg/s)`, { P_H, P_L });
      const states = { evap_out: evapOut, suction, discharge, cond_dew: condDew, cond_bub: condBub, cond_out: condOut, subcool_out: scOut,
        ihx_out: ihxOut, evap_in: evapIn, evap_dew: evapDew, isentropic: isen };
      return { success: true, message: "success", COP: cop, P_H, P_L, P_dis, P_evin, report: r, states, warnings };
    } catch (e) {
      const raw = e && e.message ? e.message : String(e);
      const Th = c.T_cond != null ? c.T_cond : c.T_steam + c.dTmin_steam;
      const near = isFinite(pr.T_crit) ? (Th > pr.T_crit - 10 ? ` — condensing at ${Th.toFixed(1)} C is within ${(pr.T_crit - Th).toFixed(1)} K of the critical point (${pr.T_crit.toFixed(1)} C)` : "")
        : (pr.isMix ? ` — the mixture may be near its critical point at ${Th.toFixed(1)} C condensing` : "");
      return fail(`property calculation failed${near}: ${raw}`);
    }
  }


  /* Diagram paths and heat-exchange profiles for a solved point (n points per segment). */
  function paths(pr, res, n = 10) {
    const st = res.states, P_H = res.P_H, P_L = res.P_L, P_dis = res.P_dis ?? P_H, P_evin = res.P_evin ?? P_L;
    const pt = s => [s.T, s.P * 10, s.h, s.s];
    const r = res.report, ihx = r.ihx && r.q_ihx_kJkg > 1e-9;
    const h1 = st.suction.h, h2 = st.discharge.h, h3 = st.cond_out.h, h4 = st.subcool_out.h, hi = st.ihx_out.h, heo = st.evap_out.h;
    const scl = [h3, h4], ihxLiq = [h4, hi];
    const hv = ihx ? hi : h4;
    // pressure along the condenser and evaporator, linear in enthalpy between inlet and (pinned) exit
    const lin = (Pa, ha, Pb, hb) => h => Math.abs(hb - ha) < 1e-12 ? Pb : Pa + (Pb - Pa) * (h - ha) / (hb - ha);
    const Pc = lin(P_dis, h2, P_H, h3), Pe = lin(P_evin, hv, P_L, heo);
    const iso = (P, a, b, k = n) => Array.from({ length: k + 1 }, (_, i) => { const h = a + (b - a) * i / k; return pt(pr.PH(typeof P === "function" ? P(h) : P, h)); });
    // the condenser through its exact dew and bubble points (sharp corners), with each part resolved;
    // a supercritical gas cooler has no corners
    const sc = !!r.supercritical;
    const satAlong = (P, q, h0) => { let h = h0; for (let k = 0; k < 4; k++) h = pr.PQ(P(h), q).h; return h; };
    const hBub0 = st.cond_bub ? st.cond_bub.h : pr.PQ(P_H, 0).h;
    const hDew = sc || P_dis === P_H ? st.cond_dew.h : satAlong(Pc, 1, st.cond_dew.h), hBub = sc || P_dis === P_H ? hBub0 : satAlong(Pc, 0, hBub0);
    const isoParts = (P, parts) => parts.filter(([a, b]) => Math.abs(b - a) > 1e-9)
      .flatMap(([a, b, k], j) => iso(P, a, b, k).slice(j ? 1 : 0));
    const comp = Array.from({ length: n + 1 }, (_, i) => { const f = i / n; const P = Math.exp(Math.log(P_L) + f * Math.log(P_dis / P_L)); return pt(pr.PH(P, h1 + f * (h2 - h1))); });
    const valve = Array.from({ length: n + 1 }, (_, i) => { const P = Math.exp(Math.log(P_H) + i / n * Math.log(P_evin / P_H)); return pt(pr.PH(P, hv)); });
    const out = [
      { c: "Compressor", pts: comp }, { c: "Condenser", pts: isoParts(Pc, sc ? [[h2, h3, 3 * n]] : [[h2, hDew, n], [hDew, hBub, 2 * n], [hBub, h3, n]]) }, { c: "Subcooler", pts: iso(P_H, scl[0], scl[1]) },
      { c: "Expansion valve", pts: valve }, { c: "Evaporator", pts: iso(Pe, hv, heo, 2 * n) },
    ];
    if (ihx) { out.push({ c: "Internal HX", pts: iso(P_H, ihxLiq[0], ihxLiq[1], n / 2) }); out.push({ c: "Internal HX", pts: iso(P_L, heo, h1, n / 2) }); }
    // profiles: T vs cumulative kW from each exchanger's cold end
    const m = r.m_dot_kgs, hot = [];
    // from the cold end: subcooler, liquid subcooled in the condenser, condensing (bubble to dew), desuperheating;
    // each segment ends exactly at the bubble and dew points
    const segsHot = sc ? [["Subcooling", scl[1], scl[0], n], ["Gas cooling", h3, h2, 3 * n]]
      : [["Subcooling", scl[1], scl[0], n], ["Subcooled liquid", h3, hBub, n], ["Condensing", hBub, hDew, 2 * n], ["Desuperheating", hDew, h2, n]];
    for (const [sec, a, b, k] of segsHot) {
      if (Math.abs(b - a) < 1e-9) continue;
      for (let i = 0; i <= k; i++) { const h = a + (b - a) * i / k; const q = sec === "Subcooling" ? m * (h - scl[1]) : m * ((scl[0] - scl[1]) + (h - h3)); hot.push([q, pr.PH(sec === "Subcooling" ? P_H : Pc(h), h).T, sec]); }
    }
    const evap = Array.from({ length: 2 * n + 1 }, (_, i) => { const h = hv + (heo - hv) * i / (2 * n); return [m * (h - hv), pr.PH(Pe(h), h).T]; });
    const prof = { hot, evap };
    if (ihx) {
      const k = n / 2, [li, lo] = ihxLiq;
      prof.ihx_liq = Array.from({ length: k + 1 }, (_, i) => [m * (li - lo) * i / k, pr.PH(P_H, lo + (li - lo) * i / k).T]);
      prof.ihx_vap = Array.from({ length: k + 1 }, (_, i) => [m * (h1 - heo) * i / k, pr.PH(P_L, heo + (h1 - heo) * i / k).T]);
    }
    return { paths: out, prof };
  }

  /* Bubble and dew lines by temperature, stopping where the flash fails (mixture critical region). */
  function dome(pr, T_lo = -10, n = 40) {
    const T_hi = isFinite(pr.T_crit) ? pr.T_crit - 0.5 : 230;
    const liq = [], vap = [];
    for (let i = 0; i <= n; i++) {
      const T = T_lo + (T_hi - T_lo) * (1 - Math.pow(1 - i / n, 1.6));
      try {
        const l = pr.TQ(T, 0), v = pr.TQ(T, 1);
        liq.push([l.T, l.P * 10, l.h, l.s]); vap.push([v.T, v.P * 10, v.h, v.s]);
      } catch (e) { break; }
    }
    return { liq, vap };
  }

  const api = { EDGE_KINDS, EDGE_NAMES, EXTENSION, LIMITS, mapDensityRatio, extensionRatio, ratioLine, RATIO_LINES, nearestOnEnvelope, extendedEfficiencies, limitFlags, critIsochore, WATER_SAT, waterSat, steamFromFeedwater, COEFFS, DISPLACEMENT_50HZ, ENVELOPE, R600A_SAT_VAPOUR, R600A_SUCTION_15K, R600A_RATED_SUPERHEAT_K, satVapourDensity, equivalentPoint, mapPoint, fluidCurve, toReference, fromReference, envelopeInFluid, fluidGrid, bicubic, inEnvelope, compressorMap, brent, DEFAULTS, saturationPressure, solve, paths, dome };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.SGHP = api;
})(typeof self !== "undefined" ? self : this);
