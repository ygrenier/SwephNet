using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SwephNet
{
    partial class Sweph
    {

        ///* we always use Astronomical Almanac constants, if available */
        //public const double MOON_MEAN_DIST = 384400000.0;		/* in m, AA 1996, F2 */
        //public const double MOON_MEAN_INCL = 5.1453964;		/* AA 1996, D2 */
        //public const double MOON_MEAN_ECC = 0.054900489;		/* AA 1996, F2 */
        ///* #define SUN_EARTH_MRAT=  328900.561400;           Su/(Ea+Mo) AA 2006 K7 */
        //public const double SUN_EARTH_MRAT = 332946.050895;           /* Su / (Ea only) AA 2006 K7 */
        //public const double EARTH_MOON_MRAT = (1 / 0.0123000383);	/* AA 2006, K7 */
        ////#if 0
        ////#define EARTH_MOON_MRAT =81.30056		/* de406 */
        ////#endif

		/// <summary>
        /// au in meters, AA 2006 K6
		/// </summary>
        public const double AUNIT = 1.49597870691e+11;
		/// <summary>
        /// m/s, AA 1996 K6
		/// </summary>
        public const double CLIGHT = 2.99792458e+8;

        ////#if 0
        ////#define HELGRAVCONST   = 1.32712438e+20;		/* G * M(sun), m^3/sec^2, AA 1996 K6 */
        ////#endif
        //public const double HELGRAVCONST = 1.32712440017987e+20;	/* G * M(sun), m^3/sec^2, AA 2006 K6 */
        //public const double GEOGCONST = 3.98600448e+14; 		/* G * M(earth) m^3/sec^2, AA 1996 K6 */
        //public const double KGAUSS = 0.01720209895;		/* Gaussian gravitational constant K6 */
        //public const double SUN_RADIUS = (959.63 / 3600 * SwissEph.DEGTORAD);  /*  Meeus germ. p 391 */
        //public const double EARTH_RADIUS = 6378136.6;		/* AA 2006 K6 */
        ///*#define EARTH_OBLATENESS= (1.0/ 298.257223563)	 * AA 1998 K13 */
        //public const double EARTH_OBLATENESS = (1.0 / 298.25642);	/* AA 2006 K6 */
        //public const double EARTH_ROT_SPEED = (7.2921151467e-5 * 86400); /* in rad/day, expl. suppl., p 162 */

        //public const double LIGHTTIME_AUNIT = (499.0047838061 / 3600 / 24); 	/* 8.3167 minutes (days), AA 2006 K6 */

        ///* node of ecliptic measured on ecliptic 2000 */
        //public const double SSY_PLANE_NODE_E2000 = (107.582569 * SwissEph.DEGTORAD);
        ///* node of ecliptic measured on solar system rotation plane */
        //public const double SSY_PLANE_NODE = (107.58883388 * SwissEph.DEGTORAD);
        ///* inclination of ecliptic against solar system rotation plane */
        //public const double SSY_PLANE_INCL = (1.578701 * SwissEph.DEGTORAD);

        //public const double KM_S_TO_AU_CTY = 21.095;			/* km/s to AU/century */
        //public const double MOON_SPEED_INTV = 0.00005; 		/* 4.32 seconds (in days) */
        //public const double PLAN_SPEED_INTV = 0.0001; 	        /* 8.64 seconds (in days) */
        //public const double MEAN_NODE_SPEED_INTV = 0.001;
        //public const double NODE_CALC_INTV = 0.0001;
        //public const double NODE_CALC_INTV_MOSH = 0.1;
        //public const double NUT_SPEED_INTV = 0.0001;
        //public const double DEFL_SPEED_INTV = 0.0000005;

        //public const double SE_LAPSE_RATE = 0.0065;  /* deg K / m, for refraction */

    }
}
