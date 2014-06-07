using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SwephNet
{
    /// <summary>
    /// Library
    /// </summary>
    public class SweLib
    {
        #region Embedded types

        /// <summary>
        /// Precession coefficients for remote past and future.
        /// </summary>
        public enum PrecessionCoefficients
        {
            Vondrak2011,
            Williams1994,
            Simon1994,
            Laskar1986,
            Bretagnon2003,
        }

        /// <summary>
        /// IAU precession 1976 or 2003 for recent centuries.
        /// </summary>
        public enum PrecessionIAU
        {
            None,
            IAU_1976,
            IAU_2000,
            /// <summary>
            /// precession model P03
            /// </summary>
            IAU_2006
        }

        /// <summary>
        /// Nutation models
        /// </summary>
        public enum NutationModel
        {
            None,
            IAU_1980,
            /// <summary>
            /// very time consuming !
            /// </summary>
            IAU_2000A,
            /// <summary>
            /// fast, but precision of milli-arcsec
            /// </summary>
            IAU_2000B
        }

        #endregion

        #region Constants

        /// <summary>
        /// Value to convert Degrees to Radian
        /// </summary>
        public const double DEGTORAD = 0.0174532925199433;

        /// <summary>
        /// Value to convert Radian to Degrees
        /// </summary>
        public const double RADTODEG = 57.2957795130823;

        /// <summary>
        /// PI * 2
        /// </summary>
        public const Double TWOPI = Math.PI * 2.0;

        /// <summary>
        /// J2000 +/- two centuries
        /// </summary>
        const double PrecessionIAU_1976_Centuries = 2.0;
        /// <summary>
        /// J2000 +/- two centuries
        /// </summary>
        const double PrecessionIAU_2000_Centuries = 2.0;
        /// <summary>
        /// J2000 +/- 75 centuries
        /// </summary>
        /// <remarks>
        /// we use P03 for whole ephemeris
        /// </remarks>
        const double PrecessionIAU_2006_Centuries = 75.0;

        #endregion

        ISwephData _Swed;

        /// <summary>
        /// New library
        /// </summary>
        public SweLib(ISwephData swed)
        {
            _Swed = swed;
        }

        #region Calculation tools

        /// <summary>
        /// Reduce x modulo 360 degrees
        /// </summary>
        public static double DegNorm(double x)
        {
            double y;
            y = (x % 360.0);
            if (Math.Abs(y) < 1e-13) y = 0;	// Alois fix 11-dec-1999
            if (y < 0.0) y += 360.0;
            return (y);
        }

        /// <summary>
        /// Reduce x modulo TWOPI degrees
        /// </summary>
        public static double RadNorm(double x)
        {
            double y;
            y = (x % TWOPI);
            if (Math.Abs(y) < 1e-13) y = 0;	/* Alois fix 11-dec-1999 */
            if (y < 0.0) y += TWOPI;
            return (y);
        }

        /// <summary>
        /// Convert a degrees value to radians
        /// </summary>
        public static double DegToRad(double x)
        {
            return x * DEGTORAD;
        }

        #endregion

        #region Precession and ecliptic obliquity

        const double AS2R = (DEGTORAD / 3600.0);
        const double D2PI = Math.PI * 2.0;
        const double EPS0 = (84381.406 * AS2R);
        const int NPOL_PEPS = 4;
        const int NPER_PEPS = 10;
        const int NPOL_PECL = 4;
        const int NPER_PECL = 8;
        const int NPOL_PEQU = 4;
        const int NPER_PEQU = 14;

        /// <summary>
        /// for pre_peps(): polynomials
        /// </summary>
        static double[,] pepol = new double[,]{
          {+8134.017132, +84028.206305},
          {+5043.0520035, +0.3624445},
          {-0.00710733, -0.00004039},
          {+0.000000271, -0.000000110}
        };

        /// <summary>
        /// for pre_peps(): periodics 
        /// </summary>
        static double[,] peper = new double[,]{
          {+409.90, +396.15, +537.22, +402.90, +417.15, +288.92, +4043.00, +306.00, +277.00, +203.00},
          {-6908.287473, -3198.706291, +1453.674527, -857.748557, +1173.231614, -156.981465, +371.836550, -216.619040, +193.691479, +11.891524},
          {+753.872780, -247.805823, +379.471484, -53.880558, -90.109153, -353.600190, -63.115353, -28.248187, +17.703387, +38.911307},
          {-2845.175469, +449.844989, -1255.915323, +886.736783, +418.887514, +997.912441, -240.979710, +76.541307, -36.788069, -170.964086},
          {-1704.720302, -862.308358, +447.832178, -889.571909, +190.402846, -56.564991, -296.222622, -75.859952, +67.473503, +3.014055}
        };

        /// <summary>
        /// for pre_pecl(): polynomials 
        /// </summary>
        static double[,] pqpol = new double[,]{
          {+5851.607687, -1600.886300},
          {-0.1189000, +1.1689818},
          {-0.00028913, -0.00000020},
          {+0.000000101, -0.000000437}
        };

        /// <summary>
        /// for pre_pecl(): periodics
        /// </summary>
        static double[,] pqper = new double[,]{
          {708.15, 2309, 1620, 492.2, 1183, 622, 882, 547},
          {-5486.751211, -17.127623, -617.517403, 413.44294, 78.614193, -180.732815, -87.676083, 46.140315},
          {-684.66156, 2446.28388, 399.671049, -356.652376, -186.387003, -316.80007, 198.296701, 101.135679}, /* typo in publication fixed */
          {667.66673, -2354.886252, -428.152441, 376.202861, 184.778874, 335.321713, -185.138669, -120.97283},
          {-5523.863691, -549.74745, -310.998056, 421.535876, -36.776172, -145.278396, -34.74445, 22.885731}
        };

        /// <summary>
        /// for pre_pequ(): polynomials
        /// </summary>
        static double[,] xypol = new double[,]{
          {+5453.282155, -73750.930350},
          {+0.4252841, -0.7675452},
          {-0.00037173, -0.00018725},
          {-0.000000152, +0.000000231}
        };

        /// <summary>
        /// for pre_pequ(): periodics 
        /// </summary>
        static double[,] xyper = new double[,]{
          {256.75, 708.15, 274.2, 241.45, 2309, 492.2, 396.1, 288.9, 231.1, 1610, 620, 157.87, 220.3, 1200},
          {-819.940624, -8444.676815, 2600.009459, 2755.17563, -167.659835, 871.855056, 44.769698, -512.313065, -819.415595, -538.071099, -189.793622, -402.922932, 179.516345, -9.814756},
          {75004.344875, 624.033993, 1251.136893, -1102.212834, -2660.66498, 699.291817, 153.16722, -950.865637, 499.754645, -145.18821, 558.116553, -23.923029, -165.405086, 9.344131},
          {81491.287984, 787.163481, 1251.296102, -1257.950837, -2966.79973, 639.744522, 131.600209, -445.040117, 584.522874, -89.756563, 524.42963, -13.549067, -210.157124, -44.919798},
          {1558.515853, 7774.939698, -2219.534038, -2523.969396, 247.850422, -846.485643, -1393.124055, 368.526116, 749.045012, 444.704518, 235.934465, 374.049623, -171.33018, -22.899655}
        };

        /// <summary>
        /// return dpre, deps
        /// </summary>
        internal static Tuple<double, double> LdpPeps(double tjd)
        {
            int i;
            double w, a, s, c;
            int npol = NPOL_PEPS;
            int nper = NPER_PEPS;
            double t = (tjd - SweDate.J2000) / 36525.0;
            double p = 0;
            double q = 0;
            // periodic terms
            for (i = 0; i < nper; i++)
            {
                w = D2PI * t;
                a = w / peper[0, i];
                s = Math.Sin(a);
                c = Math.Cos(a);
                p += c * peper[1, i] + s * peper[3, i];
                q += c * peper[2, i] + s * peper[4, i];
            }
            // polynomial terms
            w = 1;
            for (i = 0; i < npol; i++)
            {
                p += pepol[i, 0] * w;
                q += pepol[i, 1] * w;
                w *= t;
            }
            // both to radians
            p *= AS2R;
            q *= AS2R;
            // return
            return new Tuple<double, double>(p, q);
        }

        const double OFFSET_EPS_JPLHORIZONS = (35.95);
        const double DCOR_EPS_JPL_TJD0 = 2437846.5;
        const int NDCOR_EPS_JPL = 51;
        static double[] dcor_eps_jpl = new double[]{
            36.726, 36.627, 36.595, 36.578, 36.640, 36.659, 36.731, 36.765,
            36.662, 36.555, 36.335, 36.321, 36.354, 36.227, 36.289, 36.348, 36.257, 36.163,
            35.979, 35.896, 35.842, 35.825, 35.912, 35.950, 36.093, 36.191, 36.009, 35.943,
            35.875, 35.771, 35.788, 35.753, 35.822, 35.866, 35.771, 35.732, 35.543, 35.498,
            35.449, 35.409, 35.497, 35.556, 35.672, 35.760, 35.596, 35.565, 35.510, 35.394,
            35.385, 35.375, 35.415,
        };

        /// <summary>
        /// Obliquity of the ecliptic at Julian date jd
        /// </summary>
        /// <remarks>
        /// IAU Coefficients are from:
        /// J. H. Lieske, T. Lederle, W. Fricke, and B. Morando,
        /// "Expressions for the Precession Quantities Based upon the IAU
        /// (1976) System of Astronomical Constants,"  Astronomy and Astrophysics
        /// 58, 1-16 (1977).
        /// 
        /// Before or after 200 years from J2000, the formula used is from:
        /// J. Laskar, "Secular terms of classical planetary theories
        /// using the results of general theory," Astronomy and Astrophysics
        /// 157, 59070 (1986).
        /// 
        /// Bretagnon, P. et al.: 2003, "Expressions for Precession Consistent with 
        /// the IAU 2000A Model". A&A 400,785
        /// B03  	84381.4088  	-46.836051*t  	-1667*10-7*t2  	+199911*10-8*t3  	-523*10-9*t4  	-248*10-10*t5  	-3*10-11*t6
        /// C03   84381.406  	-46.836769*t  	-1831*10-7*t2  	+20034*10-7*t3  	-576*10-9*t4  	-434*10-10*t5
        /// 
        /// See precess and page B18 of the Astronomical Almanac.
        /// </remarks>
        internal static double Epsiln(double jd, JPL.JplHorizonMode horizons)
        {
            double eps = 0;
            double T = (jd - 2451545.0) / 36525.0;
            if ((horizons & JPL.JplHorizonMode.JplHorizons) != 0 && IncludeCodeForDpsiDepsIAU1980)
            {
                eps = (((1.813e-3 * T - 5.9e-4) * T - 46.8150) * T + 84381.448) * DEGTORAD / 3600;
            }
            else if ((horizons & JPL.JplHorizonMode.JplApproximate) != 0 && !ApproximateHorizonsAstrodienst)
            {
                eps = (((1.813e-3 * T - 5.9e-4) * T - 46.8150) * T + 84381.448) * DEGTORAD / 3600;
            }
            else if (UsePrecessionIAU == PrecessionIAU.IAU_1976 && Math.Abs(T) <= PrecessionIAU_1976_Centuries)
            {
                eps = (((1.813e-3 * T - 5.9e-4) * T - 46.8150) * T + 84381.448) * DEGTORAD / 3600;
            }
            else if (UsePrecessionIAU == PrecessionIAU.IAU_2000 && Math.Abs(T) <= PrecessionIAU_2000_Centuries)
            {
                eps = (((1.813e-3 * T - 5.9e-4) * T - 46.84024) * T + 84381.406) * DEGTORAD / 3600;
            }
            else if (UsePrecessionIAU == PrecessionIAU.IAU_2006 && Math.Abs(T) <= PrecessionIAU_2006_Centuries)
            {
                eps = (((((-4.34e-8 * T - 5.76e-7) * T + 2.0034e-3) * T - 1.831e-4) * T - 46.836769) * T + 84381.406) * DEGTORAD / 3600.0;
            }
            else if (UsePrecessionCoefficient == PrecessionCoefficients.Bretagnon2003)
            {
                eps = ((((((-3e-11 * T - 2.48e-8) * T - 5.23e-7) * T + 1.99911e-3) * T - 1.667e-4) * T - 46.836051) * T + 84381.40880) * DEGTORAD / 3600.0;
            }
            else if (UsePrecessionCoefficient == PrecessionCoefficients.Simon1994)
            {
                eps = (((((2.5e-8 * T - 5.1e-7) * T + 1.9989e-3) * T - 1.52e-4) * T - 46.80927) * T + 84381.412) * DEGTORAD / 3600.0;
            }
            else if (UsePrecessionCoefficient == PrecessionCoefficients.Williams1994)
            {
                eps = ((((-1.0e-6 * T + 2.0e-3) * T - 1.74e-4) * T - 46.833960) * T + 84381.409) * DEGTORAD / 3600.0;/* */
            }
            else if (UsePrecessionCoefficient == PrecessionCoefficients.Laskar1986)
            {
                T /= 10.0;
                eps = (((((((((2.45e-10 * T + 5.79e-9) * T + 2.787e-7) * T
                + 7.12e-7) * T - 3.905e-5) * T - 2.4967e-3) * T
                - 5.138e-3) * T + 1.99925) * T - 0.0155) * T - 468.093) * T
                + 84381.448;
                eps *= DEGTORAD / 3600.0;
            }
            else
            {
                // Vondrak2011
                var tup = LdpPeps(jd);
                eps = tup.Item2;
                if ((horizons & JPL.JplHorizonMode.JplApproximate) != 0 && !ApproximateHorizonsAstrodienst)
                {
                    double tofs = (jd - DCOR_EPS_JPL_TJD0) / 365.25;
                    double dofs = OFFSET_EPS_JPLHORIZONS;
                    if (tofs < 0)
                    {
                        tofs = 0;
                        dofs = dcor_eps_jpl[0];
                    }
                    else if (tofs >= NDCOR_EPS_JPL - 1)
                    {
                        tofs = NDCOR_EPS_JPL;
                        dofs = dcor_eps_jpl[NDCOR_EPS_JPL - 1];
                    }
                    else
                    {
                        double t0 = (int)tofs;
                        double t1 = t0 + 1;
                        dofs = dcor_eps_jpl[(int)t0];
                        dofs = (tofs - t0) * (dcor_eps_jpl[(int)t0] - dcor_eps_jpl[(int)t1]) + dcor_eps_jpl[(int)t0];
                    }
                    dofs /= (1000.0 * 3600.0);
                    eps += dofs * DEGTORAD;
                }
            }
            return eps;
        }

        #endregion

        #region Nutation in longitude and obliquity

        /// <summary>
        /// Nutation in longitude and obliquity
        /// computed at Julian date J.
        /// 
        /// References:
        /// "Summary of 1980 IAU Theory of Nutation (Final Report of the
        /// IAU Working Group on Nutation)", P. K. Seidelmann et al., in
        /// Transactions of the IAU Vol. XVIII A, Reports on Astronomy,
        /// P. A. Wayman, ed.; D. Reidel Pub. Co., 1982.
        /// 
        /// "Nutation and the Earth's Rotation",
        /// I.A.U. Symposium No. 78, May, 1977, page 256.
        /// I.A.U., 1980.
        /// 
        /// Woolard, E.W., "A redevelopment of the theory of nutation",
        /// The Astronomical Journal, 58, 1-3 (1953).
        /// 
        /// This program implements all of the 1980 IAU nutation series.
        /// Results checked at 100 points against the 1986 AA; all agreed.
        /// 
        /// 
        /// - S. L. Moshier, November 1987
        ///   October, 1992 - typo fixed in nutation matrix
        /// 
        /// - D. Koch, November 1995: small changes in structure,
        ///   Corrections to IAU 1980 Series added from Expl. Suppl. p. 116
        /// 
        /// Each term in the expansion has a trigonometric
        /// argument given by
        ///   W = i*MM + j*MS + k*FF + l*DD + m*OM
        /// where the variables are defined below.
        /// The nutation in longitude is a sum of terms of the
        /// form (a + bT) * Math.Sin(W). The terms for nutation in obliquity
        /// are of the form (c + dT) * cos(W).  The coefficients
        /// are arranged in the tabulation as follows:
        /// 
        /// Coefficient:
        /// i  j  k  l  m      a      b      c     d
        /// 0, 0, 0, 0, 1, -171996, -1742, 92025, 89,
        /// The first line of the table, above, is done separately
        /// since two of the values do not fit into 16 bit integers.
        /// The values a and c are arc seconds times 10000.  b and d
        /// are arc seconds per Julian century times 100000.  i through m
        /// are integers.  See the program for interpretation of MM, MS,
        /// etc., which are mean orbital elements of the Sun and Moon.
        /// 
        /// If terms with coefficient less than X are omitted, the peak
        /// errors will be:
        /// 
        ///   omit	error,		  omit	error,
        ///   a &lt;	longitude	  c &lt;	obliquity
        /// .0005"	.0100"		.0008"	.0094"
        /// .0046	.0492		.0095	.0481
        /// .0123	.0880		.0224	.0905
        /// .0386	.1808		.0895	.1129
        /// </summary>

        #region Nutation values

        static short[] base_nt = new short[] {  
        /* LS and OC are units of 0.0001"
         *LS2 and OC2 are units of 0.00001"
         *MM,MS,FF,DD,OM, LS, LS2,OC, OC2 */
         0, 0, 0, 0, 2,  2062,  2, -895,  5,
        -2, 0, 2, 0, 1,    46,  0,  -24,  0,
         2, 0,-2, 0, 0,    11,  0,    0,  0,
        -2, 0, 2, 0, 2,    -3,  0,    1,  0,
         1,-1, 0,-1, 0,    -3,  0,    0,  0,
         0,-2, 2,-2, 1,    -2,  0,    1,  0,
         2, 0,-2, 0, 1,     1,  0,    0,  0,
         0, 0, 2,-2, 2,-13187,-16, 5736,-31,
         0, 1, 0, 0, 0,  1426,-34,   54, -1,
         0, 1, 2,-2, 2,  -517, 12,  224, -6,
         0,-1, 2,-2, 2,   217, -5,  -95,  3,
         0, 0, 2,-2, 1,   129,  1,  -70,  0,
         2, 0, 0,-2, 0,    48,  0,    1,  0,
         0, 0, 2,-2, 0,   -22,  0,    0,  0,
         0, 2, 0, 0, 0,    17, -1,    0,  0,
         0, 1, 0, 0, 1,   -15,  0,    9,  0,
         0, 2, 2,-2, 2,   -16,  1,    7,  0,
         0,-1, 0, 0, 1,   -12,  0,    6,  0,
        -2, 0, 0, 2, 1,    -6,  0,    3,  0,
         0,-1, 2,-2, 1,    -5,  0,    3,  0,
         2, 0, 0,-2, 1,     4,  0,   -2,  0,
         0, 1, 2,-2, 1,     4,  0,   -2,  0,
         1, 0, 0,-1, 0,    -4,  0,    0,  0,
         2, 1, 0,-2, 0,     1,  0,    0,  0,
         0, 0,-2, 2, 1,     1,  0,    0,  0,
         0, 1,-2, 2, 0,    -1,  0,    0,  0,
         0, 1, 0, 0, 2,     1,  0,    0,  0,
        -1, 0, 0, 1, 1,     1,  0,    0,  0,
         0, 1, 2,-2, 0,    -1,  0,    0,  0,
         0, 0, 2, 0, 2, -2274, -2,  977, -5,
         1, 0, 0, 0, 0,   712,  1,   -7,  0,
         0, 0, 2, 0, 1,  -386, -4,  200,  0,
         1, 0, 2, 0, 2,  -301,  0,  129, -1,
         1, 0, 0,-2, 0,  -158,  0,   -1,  0,
        -1, 0, 2, 0, 2,   123,  0,  -53,  0,
         0, 0, 0, 2, 0,    63,  0,   -2,  0,
         1, 0, 0, 0, 1,    63,  1,  -33,  0,
        -1, 0, 0, 0, 1,   -58, -1,   32,  0,
        -1, 0, 2, 2, 2,   -59,  0,   26,  0,
         1, 0, 2, 0, 1,   -51,  0,   27,  0,
         0, 0, 2, 2, 2,   -38,  0,   16,  0,
         2, 0, 0, 0, 0,    29,  0,   -1,  0,
         1, 0, 2,-2, 2,    29,  0,  -12,  0,
         2, 0, 2, 0, 2,   -31,  0,   13,  0,
         0, 0, 2, 0, 0,    26,  0,   -1,  0,
        -1, 0, 2, 0, 1,    21,  0,  -10,  0,
        -1, 0, 0, 2, 1,    16,  0,   -8,  0,
         1, 0, 0,-2, 1,   -13,  0,    7,  0,
        -1, 0, 2, 2, 1,   -10,  0,    5,  0,
         1, 1, 0,-2, 0,    -7,  0,    0,  0,
         0, 1, 2, 0, 2,     7,  0,   -3,  0,
         0,-1, 2, 0, 2,    -7,  0,    3,  0,
         1, 0, 2, 2, 2,    -8,  0,    3,  0,
         1, 0, 0, 2, 0,     6,  0,    0,  0,
         2, 0, 2,-2, 2,     6,  0,   -3,  0,
         0, 0, 0, 2, 1,    -6,  0,    3,  0,
         0, 0, 2, 2, 1,    -7,  0,    3,  0,
         1, 0, 2,-2, 1,     6,  0,   -3,  0,
         0, 0, 0,-2, 1,    -5,  0,    3,  0,
         1,-1, 0, 0, 0,     5,  0,    0,  0,
         2, 0, 2, 0, 1,    -5,  0,    3,  0, 
         0, 1, 0,-2, 0,    -4,  0,    0,  0,
         1, 0,-2, 0, 0,     4,  0,    0,  0,
         0, 0, 0, 1, 0,    -4,  0,    0,  0,
         1, 1, 0, 0, 0,    -3,  0,    0,  0,
         1, 0, 2, 0, 0,     3,  0,    0,  0,
         1,-1, 2, 0, 2,    -3,  0,    1,  0,
        -1,-1, 2, 2, 2,    -3,  0,    1,  0,
        -2, 0, 0, 0, 1,    -2,  0,    1,  0,
         3, 0, 2, 0, 2,    -3,  0,    1,  0,
         0,-1, 2, 2, 2,    -3,  0,    1,  0,
         1, 1, 2, 0, 2,     2,  0,   -1,  0,
        -1, 0, 2,-2, 1,    -2,  0,    1,  0,
         2, 0, 0, 0, 1,     2,  0,   -1,  0,
         1, 0, 0, 0, 2,    -2,  0,    1,  0,
         3, 0, 0, 0, 0,     2,  0,    0,  0,
         0, 0, 2, 1, 2,     2,  0,   -1,  0,
        -1, 0, 0, 0, 2,     1,  0,   -1,  0,

         1, 0, 0,-4, 0,    -1,  0,    0,  0,
        -2, 0, 2, 2, 2,     1,  0,   -1,  0,
        -1, 0, 2, 4, 2,    -2,  0,    1,  0,
         2, 0, 0,-4, 0,    -1,  0,    0,  0,
         1, 1, 2,-2, 2,     1,  0,   -1,  0,
         1, 0, 2, 2, 1,    -1,  0,    1,  0,
        -2, 0, 2, 4, 2,    -1,  0,    1,  0,
        -1, 0, 4, 0, 2,     1,  0,    0,  0,
         1,-1, 0,-2, 0,     1,  0,    0,  0,
         2, 0, 2,-2, 1,     1,  0,   -1,  0,
         2, 0, 2, 2, 2,    -1,  0,    0,  0,
         1, 0, 0, 2, 1,    -1,  0,    0,  0,
         0, 0, 4,-2, 2,     1,  0,    0,  0,
         3, 0, 2,-2, 2,     1,  0,    0,  0,
         1, 0, 2,-2, 0,    -1,  0,    0,  0,
         0, 1, 2, 0, 1,     1,  0,    0,  0,
        -1,-1, 0, 2, 1,     1,  0,    0,  0,
         0, 0,-2, 0, 1,    -1,  0,    0,  0,
         0, 0, 2,-1, 2,    -1,  0,    0,  0,
         0, 1, 0, 2, 0,    -1,  0,    0,  0,
         1, 0,-2,-2, 0,    -1,  0,    0,  0,
         0,-1, 2, 0, 1,    -1,  0,    0,  0,
         1, 1, 0,-2, 1,    -1,  0,    0,  0,
         1, 0,-2, 2, 0,    -1,  0,    0,  0,
         2, 0, 0, 2, 0,     1,  0,    0,  0,
         0, 0, 2, 4, 2,    -1,  0,    0,  0,
         0, 1, 0, 1, 0,     1,  0,    0,  0,
        };

        static short[] corr_1987_nt = new short[] {  
        /* corrections to IAU 1980 nutation series by Herring 1987
         *             in 0.00001" !!!
         *              LS      OC      */
         101, 0, 0, 0, 1,-725, 0, 213, 0,
         101, 1, 0, 0, 0, 523, 0, 208, 0,
         101, 0, 2,-2, 2, 102, 0, -41, 0,
         101, 0, 2, 0, 2, -81, 0,  32, 0,
        /*              LC      OS !!!  */
         102, 0, 0, 0, 1, 417, 0, 224, 0,
         102, 1, 0, 0, 0,  61, 0, -24, 0,
         102, 0, 2,-2, 2,-118, 0, -47, 0,
        };

        static short[] _nt = null;
        static short[] nt
        {
            get
            {
                if (_nt == null)
                {
                    if (NutationCorrection1987)
                    {
                        _nt = base_nt.Concat(corr_1987_nt).ToArray();
                    }
                    else
                    {
                        _nt = base_nt;
                    }
                }
                return _nt;
            }
        }

        #endregion

        static double[] NutationIAU1980(double J)
        {
            // arrays to hold sines and cosines of multiple angles
            double[,] ss = new double[5, 8];
            double[,] cc = new double[5, 8];
            double[] args = new double[5];
            int[] ns = new int[5];

            // Julian centuries from 2000 January 1.5,
            // barycentric dynamical time
            double T = (J - 2451545.0) / 36525.0;
            double T2 = T * T;

            // Fundamental arguments in the FK5 reference system.
            // The coefficients, originally given to 0.001",
            // are converted here to degrees.
            // longitude of the mean ascending node of the lunar orbit
            // on the ecliptic, measured from the mean equinox of date
            double OM = -6962890.539 * T + 450160.280 + (0.008 * T + 7.455) * T2;
            OM = DegNorm(OM / 3600) * DEGTORAD;

            // mean longitude of the Sun minus the
            // mean longitude of the Sun's perigee
            double MS = 129596581.224 * T + 1287099.804 - (0.012 * T + 0.577) * T2;
            MS = DegNorm(MS / 3600) * DEGTORAD;

            // mean longitude of the Moon minus the
            // mean longitude of the Moon's perigee
            double MM = 1717915922.633 * T + 485866.733 + (0.064 * T + 31.310) * T2;
            MM = DegNorm(MM / 3600) * DEGTORAD;

            // mean longitude of the Moon minus the
            // mean longitude of the Moon's node
            double FF = 1739527263.137 * T + 335778.877 + (0.011 * T - 13.257) * T2;
            FF = DegNorm(FF / 3600) * DEGTORAD;

            // mean elongation of the Moon from the Sun.
            double DD = 1602961601.328 * T + 1072261.307 + (0.019 * T - 6.891) * T2;
            DD = DegNorm(DD / 3600) * DEGTORAD;
            args[0] = MM;
            ns[0] = 3;
            args[1] = MS;
            ns[1] = 2;
            args[2] = FF;
            ns[2] = 4;
            args[3] = DD;
            ns[3] = 4;
            args[4] = OM;
            ns[4] = 2;

            // Calculate Math.Sin( i*MM ), etc. for needed multiple angles
            for (int k = 0; k <= 4; k++)
            {
                double arg = args[k];
                int n = ns[k];
                double su = Math.Sin(arg);
                double cu = Math.Cos(arg);
                ss[k, 0] = su;			/* sin(L) */
                cc[k, 0] = cu;			/* cos(L) */
                double sv = 2.0 * su * cu;
                double cv = cu * cu - su * su;
                ss[k, 1] = sv;			/* sin(2L) */
                cc[k, 1] = cv;
                for (int i = 2; i < n; i++)
                {
                    double s = su * cv + cu * sv;
                    cv = cu * cv - su * sv;
                    sv = s;
                    ss[k, i] = sv;		/* sin( i+1 L ) */
                    cc[k, i] = cv;
                }
            }

            // first terms, not in table: */
            double C = (-0.01742 * T - 17.1996) * ss[4, 0];	/* sin(OM) */
            double D = (0.00089 * T + 9.2025) * cc[4, 0];	/* cos(OM) */
            var wnt = nt;   // Working reference to nt
            for (int p = 0, pcnt = wnt.Length; p < pcnt; p += 9)
            {
                // argument of sine and cosine
                int k1 = 0;
                double cv = 0.0;
                double sv = 0.0;
                for (int m = 0; m < 5; m++)
                {
                    int j = wnt[p + m];
                    if (j > 100)
                        j = 0; // p[0] is a flag
                    if (j != 0)
                    {
                        int k = j;
                        if (j < 0)
                            k = -k;
                        double su = ss[m, k - 1]; /* sin(k*angle) */
                        if (j < 0)
                            su = -su;
                        double cu = cc[m, k - 1];
                        if (k1 == 0)
                        { /* set first angle */
                            sv = su;
                            cv = cu;
                            k1 = 1;
                        }
                        else
                        {		/* combine angles */
                            double sw = su * cv + cu * sv;
                            cv = cu * cv - su * sv;
                            sv = sw;
                        }
                    }
                }

                // longitude coefficient, in 0.0001"
                double f = wnt[p + 5] * 0.0001;
                if (wnt[p + 6] != 0)
                    f += 0.00001 * T * wnt[p + 6];

                // obliquity coefficient, in 0.0001"
                double g = wnt[p + 7] * 0.0001;
                if (wnt[p + 8] != 0)
                    g += 0.00001 * T * wnt[p + 8];
                if (p >= 100)
                { 	/* coefficients in 0.00001" */
                    f *= 0.1;
                    g *= 0.1;
                }

                // accumulate the terms
                if (p != 102)
                {
                    C += f * sv;
                    D += g * cv;
                }
                else
                { 		/* cos for nutl and sin for nuto */
                    C += f * cv;
                    D += g * sv;
                }
            }

            // Save answers, expressed in radians
            double[] result = new double[2];
            result[0] = DEGTORAD * C / 3600.0;
            result[1] = DEGTORAD * D / 3600.0;

            //  nutlo[0] += (-0.071590 / 3600.0) * DEGTORAD;
            //  nutlo[1] += (-0.008000 / 3600.0) * DEGTORAD;

            //  nutlo[0] += (-0.047878 / 3600.0) * DEGTORAD;
            //  nutlo[1] += (-0.004035 / 3600.0) * DEGTORAD;

            return result;
        }

        /// <summary>
        /// Nutation IAU 2000A model 
        /// (MHB2000 luni-solar and planetary nutation, without free core nutation)
        ///
        /// Function returns nutation in longitude and obliquity in radians with 
        /// respect to the equinox of date. For the obliquity of the ecliptic
        /// the calculation of Lieske & al. (1977) must be used.
        ///
        /// The precision in recent years is about 0.001 arc seconds.
        ///
        /// The calculation includes luni-solar and planetary nutation.
        /// Free core nutation, which cannot be predicted, is omitted, 
        /// the error being of the order of a few 0.0001 arc seconds.
        ///
        /// References:
        ///
        /// Capitaine, N., Wallace, P.T., Chapront, J., A & A 432, 366 (2005).
        ///
        /// Chapront, J., Chapront-Touze, M. & Francou, G., A & A 387, 700 (2002).
        ///
        /// Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
        /// for the precession quantities based upon the IAU (1976) System of
        /// Astronomical Constants", A & A 58, 1-16 (1977).
        ///
        /// Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
        /// and precession   New nutation series for nonrigid Earth and
        /// insights into the Earth's interior", J.Geophys.Res., 107, B4,
        /// 2002.  
        ///
        /// Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
        /// Francou, G., Laskar, J., A & A 282, 663-683 (1994).
        ///
        /// Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M., A & A Supp.
        /// Ser. 135, 111 (1999).
        ///
        /// Wallace, P.T., "Software for Implementing the IAU 2000
        /// Resolutions", in IERS Workshop 5.1 (2002).
        ///
        /// Nutation IAU 2000A series in: 
        /// Kaplan, G.H., United States Naval Observatory Circular No. 179 (Oct. 2005)
        /// aa.usno.navy.mil/publications/docs/Circular_179.html
        ///
        /// MHB2000 code at
        /// - ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
        /// - http://www.iau-sofa.rl.ac.uk/2005_0901/Downloads.html
        /// </summary>
        static double[] NutationIAU2000ab(double J)
        {
            double dpsi = 0, deps = 0;

            double T = (J - SweDate.J2000) / 36525.0;

            // luni-solar nutation
            // Fundamental arguments, Simon & al. (1994)
            // Mean anomaly of the Moon.
            double M = DegNorm((485868.249036 +
                    T * (1717915923.2178 +
                    T * (31.8792 +
                    T * (0.051635 +
                    T * (-0.00024470))))) / 3600.0) * DEGTORAD;

            // Mean anomaly of the Sun
            double SM = DegNorm((1287104.79305 +
                    T * (129596581.0481 +
                    T * (-0.5532 +
                    T * (0.000136 +
                    T * (-0.00001149))))) / 3600.0) * DEGTORAD;

            // Mean argument of the latitude of the Moon.
            double F = DegNorm((335779.526232 +
                    T * (1739527262.8478 +
                    T * (-12.7512 +
                    T * (-0.001037 +
                    T * (0.00000417))))) / 3600.0) * DEGTORAD;

            // Mean elongation of the Moon from the Sun.
            double D = DegNorm((1072260.70369 +
                    T * (1602961601.2090 +
                    T * (-6.3706 +
                    T * (0.006593 +
                    T * (-0.00003169))))) / 3600.0) * DEGTORAD;

            //    /* Mean longitude of the ascending node of the Moon. */
            double OM = DegNorm((450160.398036 +
                    T * (-6962890.5431 +
                    T * (7.4722 +
                    T * (0.007702 +
                    T * (-0.00005939))))) / 3600.0) * DEGTORAD;

            // luni-solar nutation series, in reverse order, starting with small terms
            var nls = Lib.SweNut200a.nls;
            var cls = Lib.SweNut200a.cls;
            int inls;
            if (UseNutationModel == NutationModel.IAU_2000B)
                inls = Lib.SweNut200a.NLS_2000B;
            else
                inls = Lib.SweNut200a.NLS;
            for (int i = inls - 1; i >= 0; i--)
            {
                int j = i * 5;
                double darg = RadNorm((double)nls[j + 0] * M +
                           (double)nls[j + 1] * SM +
                           (double)nls[j + 2] * F +
                           (double)nls[j + 3] * D +
                           (double)nls[j + 4] * OM);
                double sinarg = Math.Sin(darg);
                double cosarg = Math.Cos(darg);
                int k = i * 6;
                dpsi += (cls[k + 0] + cls[k + 1] * T) * sinarg + cls[k + 2] * cosarg;
                deps += (cls[k + 3] + cls[k + 4] * T) * cosarg + cls[k + 5] * sinarg;
            }

            double[] result = new double[]{
                dpsi * Lib.SweNut200a.O1MAS2DEG,
                deps * Lib.SweNut200a.O1MAS2DEG
            };

            if (UseNutationModel == NutationModel.IAU_2000A)
            {
                // planetary nutation 
                // note: The MHB2000 code computes the luni-solar and planetary nutation
                // in different routines, using slightly different Delaunay
                // arguments in the two cases.  This behaviour is faithfully
                // reproduced here.  Use of the Simon et al. expressions for both
                // cases leads to negligible changes, well below 0.1 microarcsecond.

                // Mean anomaly of the Moon.
                double AL = RadNorm(2.35555598 + 8328.6914269554 * T);

                // Mean anomaly of the Sun.
                double ALSU = RadNorm(6.24006013 + 628.301955 * T);

                // Mean argument of the latitude of the Moon.
                double AF = RadNorm(1.627905234 + 8433.466158131 * T);

                // Mean elongation of the Moon from the Sun.
                double AD = RadNorm(5.198466741 + 7771.3771468121 * T);

                // Mean longitude of the ascending node of the Moon.
                double AOM = RadNorm(2.18243920 - 33.757045 * T);

                // Planetary longitudes, Mercury through Neptune (Souchay et al. 1999).
                double ALME = RadNorm(4.402608842 + 2608.7903141574 * T);
                double ALVE = RadNorm(3.176146697 + 1021.3285546211 * T);
                double ALEA = RadNorm(1.753470314 + 628.3075849991 * T);
                double ALMA = RadNorm(6.203480913 + 334.0612426700 * T);
                double ALJU = RadNorm(0.599546497 + 52.9690962641 * T);
                double ALSA = RadNorm(0.874016757 + 21.3299104960 * T);
                double ALUR = RadNorm(5.481293871 + 7.4781598567 * T);
                double ALNE = RadNorm(5.321159000 + 3.8127774000 * T);

                // General accumulated precession in longitude.
                double APA = (0.02438175 + 0.00000538691 * T) * T;

                //  planetary nutation series (in reverse order).
                var npl = Lib.SweNut200a.npl;
                var icpl = Lib.SweNut200a.icpl;
                dpsi = 0;
                deps = 0;
                for (int i = Lib.SweNut200a.NPL - 1; i >= 0; i--)
                {
                    int j = i * 14;
                    double darg = RadNorm((double)npl[j + 0] * AL +
                    (double)npl[j + 1] * ALSU +
                    (double)npl[j + 2] * AF +
                    (double)npl[j + 3] * AD +
                    (double)npl[j + 4] * AOM +
                    (double)npl[j + 5] * ALME +
                    (double)npl[j + 6] * ALVE +
                    (double)npl[j + 7] * ALEA +
                    (double)npl[j + 8] * ALMA +
                    (double)npl[j + 9] * ALJU +
                    (double)npl[j + 10] * ALSA +
                    (double)npl[j + 11] * ALUR +
                    (double)npl[j + 12] * ALNE +
                    (double)npl[j + 13] * APA);
                    int k = i * 4;
                    double sinarg = Math.Sin(darg);
                    double cosarg = Math.Cos(darg);
                    dpsi += (double)icpl[k + 0] * sinarg + (double)icpl[k + 1] * cosarg;
                    deps += (double)icpl[k + 2] * sinarg + (double)icpl[k + 3] * cosarg;
                }
                result[0] += dpsi * Lib.SweNut200a.O1MAS2DEG;
                result[1] += deps * Lib.SweNut200a.O1MAS2DEG;

                // changes required by adoption of P03 precession 
                // according to Capitaine et al. A & A 412, 366 (2005) = IAU 2006
                dpsi = -8.1 * Math.Sin(OM) - 0.6 * Math.Sin(2 * F - 2 * D + 2 * OM);
                dpsi += T * (47.8 * Math.Sin(OM) + 3.7 * Math.Sin(2 * F - 2 * D + 2 * OM) + 0.6 * Math.Sin(2 * F + 2 * OM) - 0.6 * Math.Sin(2 * OM));
                deps = T * (-25.6 * Math.Cos(OM) - 1.6 * Math.Cos(2 * F - 2 * D + 2 * OM));
                result[0] += dpsi / (3600.0 * 1000000.0);
                result[1] += deps / (3600.0 * 1000000.0);
            }

            // Result
            result[0] *= DEGTORAD;
            result[1] *= DEGTORAD;
            return result;
        }


        static double Bessel(double[] v, int n, double t)
        {
            //    double B; double[] d = new double[6];
            if (t <= 0)
            {
                return v[0];
            }
            if (t >= n - 1)
            {
                return v[n - 1];
            }
            double p = Math.Floor(t);
            int iy = (int)t;

            // Zeroth order estimate is value at start of year
            double ans = v[iy];
            int k = iy + 1;
            if (k >= n)
                return ans;

            // The fraction of tabulation interval
            p = t - p;
            ans += p * (v[k] - v[iy]);
            if ((iy - 1 < 0) || (iy + 2 >= n))
                return ans; // can't do second differences

            // Make table of first differences
            k = iy - 2;
            var d = new double[6];
            for (int i = 0; i < 5; i++)
            {
                if ((k < 0) || (k + 1 >= n))
                    d[i] = 0;
                else
                    d[i] = v[k + 1] - v[k];
                k += 1;
            }

            // Compute second differences
            for (int i = 0; i < 4; i++)
                d[i] = d[i + 1] - d[i];
            double B = 0.25 * p * (p - 1.0);
            ans += B * (d[1] + d[2]);
            if (iy + 2 >= n)
                return ans;

            // Compute third differences
            for (int i = 0; i < 3; i++)
                d[i] = d[i + 1] - d[i];
            B = 2.0 * B / 3.0;
            ans += (p - 0.5) * B * d[1];
            if ((iy - 2 < 0) || (iy + 3 > n))
                return ans;

            // Compute fourth differences
            for (int i = 0; i < 2; i++)
                d[i] = d[i + 1] - d[i];
            B = 0.125 * B * (p + 1.0) * (p - 2.0);
            ans += B * (d[0] + d[1]);

            return ans;
        }

        /// <summary>
        /// Nutation
        /// </summary>
        internal double[] Nutation(double J, JPL.JplHorizonMode jplMode)
        {
            double[] result = new double[2];
            if ((jplMode & JPL.JplHorizonMode.JplHorizons) != 0 && IncludeCodeForDpsiDepsIAU1980)
            {
                result = NutationIAU1980(J);
            }
            else if (UseNutationModel == NutationModel.IAU_1980)
            {
                result = NutationIAU1980(J);
            }
            else if (UseNutationModel == NutationModel.IAU_2000A || UseNutationModel == NutationModel.IAU_2000B)
            {
                result = NutationIAU2000ab(J);
                if ((jplMode & JPL.JplHorizonMode.JplApproximate) != 0 && !ApproximateHorizonsAstrodienst)
                {
                    result[0] += -41.7750 / 3600.0 / 1000.0 * DEGTORAD;
                    result[1] += -6.8192 / 3600.0 / 1000.0 * DEGTORAD;
                }
            }
            if (IncludeCodeForDpsiDepsIAU1980)
            {
                if ((jplMode & JPL.JplHorizonMode.JplHorizons) != 0)
                {
                    int n = (int)(_Swed.EopTjdEnd - _Swed.EopTjdBeg + 0.000001);
                    double J2 = J;
                    if (J < _Swed.EopTjdBeg_horizons)
                        J2 = _Swed.EopTjdBeg_horizons;
                    double dpsi = Bessel(_Swed.Dpsi, n + 1, J2 - _Swed.EopTjdBeg);
                    double deps = Bessel(_Swed.Deps, n + 1, J2 - _Swed.EopTjdBeg);
                    result[0] += dpsi / 3600.0 * DEGTORAD;
                    result[1] += deps / 3600.0 * DEGTORAD;
                }
            }
            return result;
        }

        #endregion

        #region Parameters

        /// <summary>
        /// Precession coefficients for remote past and future
        /// </summary>
        public static PrecessionCoefficients UsePrecessionCoefficient = PrecessionCoefficients.Vondrak2011;

        /// <summary>
        /// IAU precession 1976 or 2003 for recent centuries.
        /// </summary>
        public static PrecessionIAU UsePrecessionIAU = PrecessionIAU.None;

        /// <summary>
        /// Nutation model
        /// </summary>
        public static NutationModel UseNutationModel = NutationModel.IAU_2000B;

        /// <summary>
        /// You can set the latter false if you do not want to compile the
        /// code required to reproduce JPL Horizons.
        /// Keep it TRUE in order to reproduce JPL Horizons following
        /// IERS Conventions 1996 (1992), p. 22. Call swe_calc_ut() with 
        /// iflag|SEFLG_JPLHOR.  This options runs only, if the files 
        /// DPSI_DEPS_IAU1980_FILE_EOPC04 and DPSI_DEPS_IAU1980_FILE_FINALS
        /// are in the ephemeris path.
        /// </summary>
        public static bool IncludeCodeForDpsiDepsIAU1980 = true;

        /// <summary>
        /// If the above define INCLUDE_CODE_FOR_DPSI_DEPS_IAU1980 is FALSE or 
        /// the software does not find the earth orientation files (see above)
        /// in the ephemeris path, then SEFLG_JPLHOR will run as 
        /// SEFLG_JPLHOR_APPROX.
        /// The following define APPROXIMATE_HORIZONS_ASTRODIENST defines 
        /// the handling of SEFLG_JPLHOR_APPROX.
        /// With this flag, planetary positions are always calculated 
        /// using a recent precession/nutation model.  
        /// If APPROXIMATE_HORIZONS_ASTRODIENST is FALSE, then the 
        /// frame bias as recommended by IERS Conventions 2003 and 2010
        /// is *not* applied. Instead, dpsi_bias and deps_bias are added to 
        /// nutation. This procedure is found in some older astronomical software.
        /// Equatorial apparent positions will be close to JPL Horizons 
        /// (within a few mas) beetween 1962 and current years. Ecl. longitude 
        /// will be good, latitude bad.
        /// If APPROXIMATE_HORIZONS_ASTRODIENST is TRUE, the approximation of 
        /// JPL Horizons is even better. Frame bias matrix is applied with
        /// some correction to RA and another correction is added to epsilon.
        /// </summary>
        public static bool ApproximateHorizonsAstrodienst = true;

        /// <summary>
        /// The latter, if combined with SEFLG_JPLHOR provides good agreement 
        /// with JPL Horizons for 1800 - today. However, Horizons uses correct
        /// dpsi and deps only after 20-jan-1962. For all dates before that
        /// it uses dpsi and deps of 20-jan-1962, which provides a continuous 
        /// ephemeris, but does not make sense otherwise. 
        /// Before 1800, even this option does not provide agreement with Horizons,
        /// because Horizons uses a different precession model (Owen 1986) 
        /// before 1800, which is not included in the Swiss Ephemeris.
        /// If this macro is FALSE then the program defaults to SEFLG_JPLHOR_APPROX
        /// outside the time range of correction data dpsi and deps.
        /// Note that this will result in a non-continuous ephemeris near
        /// 20-jan-1962 and current years.
        /// </summary>
        /// <remarks>
        /// Horizons method before 20-jan-1962
        /// </remarks>
        public static bool UseHorizonsMethodBefore1980 = true;

        /// <summary>
        /// Set TRUE, to include Herring's (1987) corrections to IAU 1980 
        /// nutation series. AA (1996) neglects them.
        /// </summary>
        public static bool NutationCorrection1987
        {
            get { return _NutationCorrection1987; }
            set { _NutationCorrection1987 = value; _nt = null; }
        }
        static bool _NutationCorrection1987 = false;

        #endregion

    }
}
