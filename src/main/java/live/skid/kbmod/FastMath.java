package lib.net.jafama;

import lib.net.jafama.CmnFastMath;
import lib.net.jafama.DoubleWrapper;
import lib.net.jafama.IntWrapper;
import lib.net.jafama.NumbersUtils;

public final class FastMath
extends CmnFastMath {
    private static final boolean USE_JDK_MATH = FM_USE_JDK_MATH;
    private static final boolean USE_REDEFINED_LOG = FM_USE_REDEFINED_LOG;
    private static final boolean USE_REDEFINED_SQRT = FM_USE_REDEFINED_SQRT;
    private static final boolean USE_POWTABS_FOR_ASIN = false;

    public static double sin(double angle) {
        if (USE_JDK_MATH) {
            return Math.sin(angle);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (angle > SIN_COS_MAX_VALUE_FOR_INT_MODULO) {
            long remAndQuad = FastMath.remainderPiO2(angle);
            angle = FastMath.decodeRemainder(remAndQuad);
            int q = FastMath.decodeQuadrant(remAndQuad);
            double sin = q == 0 ? FastMath.sin(angle) : (q == 1 ? FastMath.cos(angle) : (q == 2 ? -FastMath.sin(angle) : -FastMath.cos(angle)));
            return negateResult ? -sin : sin;
        }
        int index = (int)(angle * SIN_COS_INDEXER + 0.5);
        double delta = angle - (double)index * SIN_COS_DELTA_HI - (double)index * SIN_COS_DELTA_LO;
        double indexSin = CmnFastMath.MyTSinCos.sinTab[index &= SIN_COS_TABS_SIZE - 2];
        double indexCos = CmnFastMath.MyTSinCos.cosTab[index];
        double result = indexSin + delta * (indexCos + delta * (-indexSin * 0.5 + delta * (-indexCos * 0.16666666666666666 + delta * indexSin * 0.041666666666666664)));
        return negateResult ? -result : result;
    }

    public static double sinQuick(double angle) {
        if (USE_JDK_MATH) {
            return Math.sin(angle);
        }
        return CmnFastMath.MyTSinCos.cosTab[(int)(Math.abs(angle - 1.5707963267948966) * SIN_COS_INDEXER + 0.5) & SIN_COS_TABS_SIZE - 2];
    }

    public static double cos(double angle) {
        if (USE_JDK_MATH) {
            return Math.cos(angle);
        }
        if ((angle = Math.abs(angle)) > SIN_COS_MAX_VALUE_FOR_INT_MODULO) {
            long remAndQuad = FastMath.remainderPiO2(angle);
            angle = FastMath.decodeRemainder(remAndQuad);
            int q = FastMath.decodeQuadrant(remAndQuad);
            double cos = q == 0 ? FastMath.cos(angle) : (q == 1 ? -FastMath.sin(angle) : (q == 2 ? -FastMath.cos(angle) : FastMath.sin(angle)));
            return cos;
        }
        int index = (int)(angle * SIN_COS_INDEXER + 0.5);
        double delta = angle - (double)index * SIN_COS_DELTA_HI - (double)index * SIN_COS_DELTA_LO;
        double indexCos = CmnFastMath.MyTSinCos.cosTab[index &= SIN_COS_TABS_SIZE - 2];
        double indexSin = CmnFastMath.MyTSinCos.sinTab[index];
        return indexCos + delta * (-indexSin + delta * (-indexCos * 0.5 + delta * (indexSin * 0.16666666666666666 + delta * indexCos * 0.041666666666666664)));
    }

    public static double cosQuick(double angle) {
        if (USE_JDK_MATH) {
            return Math.cos(angle);
        }
        return CmnFastMath.MyTSinCos.cosTab[(int)(Math.abs(angle) * SIN_COS_INDEXER + 0.5) & SIN_COS_TABS_SIZE - 2];
    }

    public static double sinAndCos(double angle, DoubleWrapper cosine) {
        if (USE_JDK_MATH) {
            cosine.value = Math.cos(angle);
            return Math.sin(angle);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (angle > SIN_COS_MAX_VALUE_FOR_INT_MODULO) {
            double sin;
            long remAndQuad = FastMath.remainderPiO2(angle);
            angle = FastMath.decodeRemainder(remAndQuad);
            int q = FastMath.decodeQuadrant(remAndQuad);
            if (q == 0) {
                sin = FastMath.sin(angle);
                cosine.value = FastMath.cos(angle);
            } else if (q == 1) {
                sin = FastMath.cos(angle);
                cosine.value = -FastMath.sin(angle);
            } else if (q == 2) {
                sin = -FastMath.sin(angle);
                cosine.value = -FastMath.cos(angle);
            } else {
                sin = -FastMath.cos(angle);
                cosine.value = FastMath.sin(angle);
            }
            return negateResult ? -sin : sin;
        }
        int index = (int)(angle * SIN_COS_INDEXER + 0.5);
        double delta = angle - (double)index * SIN_COS_DELTA_HI - (double)index * SIN_COS_DELTA_LO;
        double indexSin = CmnFastMath.MyTSinCos.sinTab[index &= SIN_COS_TABS_SIZE - 2];
        double indexCos = CmnFastMath.MyTSinCos.cosTab[index];
        cosine.value = indexCos + delta * (-indexSin + delta * (-indexCos * 0.5 + delta * (indexSin * 0.16666666666666666 + delta * indexCos * 0.041666666666666664)));
        double result = indexSin + delta * (indexCos + delta * (-indexSin * 0.5 + delta * (-indexCos * 0.16666666666666666 + delta * indexSin * 0.041666666666666664)));
        return negateResult ? -result : result;
    }

    public static double tan(double angle) {
        double result;
        if (USE_JDK_MATH) {
            return Math.tan(angle);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (angle > TAN_MAX_VALUE_FOR_INT_MODULO && (angle = FastMath.remainderPi(angle)) < 0.0) {
            angle = -angle;
            negateResult = !negateResult;
        }
        int index = (int)(angle * TAN_INDEXER + 0.5);
        double delta = angle - (double)index * TAN_DELTA_HI - (double)index * TAN_DELTA_LO;
        if ((index &= 2 * (TAN_VIRTUAL_TABS_SIZE - 1) - 1) > TAN_VIRTUAL_TABS_SIZE - 1) {
            index = 2 * (TAN_VIRTUAL_TABS_SIZE - 1) - index;
            delta = -delta;
            boolean bl = negateResult = !negateResult;
        }
        if (index < TAN_TABS_SIZE) {
            result = CmnFastMath.MyTTan.tanTab[index] + delta * (CmnFastMath.MyTTan.tanDer1DivF1Tab[index] + delta * (CmnFastMath.MyTTan.tanDer2DivF2Tab[index] + delta * (CmnFastMath.MyTTan.tanDer3DivF3Tab[index] + delta * CmnFastMath.MyTTan.tanDer4DivF4Tab[index])));
        } else {
            index = TAN_VIRTUAL_TABS_SIZE - 1 - index;
            result = 1.0 / (CmnFastMath.MyTTan.tanTab[index] - delta * (CmnFastMath.MyTTan.tanDer1DivF1Tab[index] - delta * (CmnFastMath.MyTTan.tanDer2DivF2Tab[index] - delta * (CmnFastMath.MyTTan.tanDer3DivF3Tab[index] - delta * CmnFastMath.MyTTan.tanDer4DivF4Tab[index]))));
        }
        return negateResult ? -result : result;
    }

    public static double asin(double value) {
        if (USE_JDK_MATH) {
            return Math.asin(value);
        }
        boolean negateResult = false;
        if (value < 0.0) {
            value = -value;
            negateResult = true;
        }
        if (value <= ASIN_MAX_VALUE_FOR_TABS) {
            int index = (int)(value * ASIN_INDEXER + 0.5);
            double delta = value - (double)index * ASIN_DELTA;
            double result = CmnFastMath.MyTAsin.asinTab[index] + delta * (CmnFastMath.MyTAsin.asinDer1DivF1Tab[index] + delta * (CmnFastMath.MyTAsin.asinDer2DivF2Tab[index] + delta * (CmnFastMath.MyTAsin.asinDer3DivF3Tab[index] + delta * CmnFastMath.MyTAsin.asinDer4DivF4Tab[index])));
            return negateResult ? -result : result;
        }
        if (value < 1.0) {
            double t = (1.0 - value) * 0.5;
            double p = t * (ASIN_PS0 + t * (ASIN_PS1 + t * (ASIN_PS2 + t * (ASIN_PS3 + t * (ASIN_PS4 + t * ASIN_PS5)))));
            double q = 1.0 + t * (ASIN_QS1 + t * (ASIN_QS2 + t * (ASIN_QS3 + t * ASIN_QS4)));
            double s = FastMath.sqrt(t);
            double z = s + s * (p / q);
            double result = ASIN_PIO2_HI - (z + z - ASIN_PIO2_LO);
            return negateResult ? -result : result;
        }
        if (value == 1.0) {
            return negateResult ? -1.5707963267948966 : 1.5707963267948966;
        }
        return Double.NaN;
    }

    public static double asinInRange(double value) {
        if (value <= -1.0) {
            return -1.5707963267948966;
        }
        if (value >= 1.0) {
            return 1.5707963267948966;
        }
        return FastMath.asin(value);
    }

    public static double acos(double value) {
        if (USE_JDK_MATH) {
            return Math.acos(value);
        }
        return 1.5707963267948966 - FastMath.asin(value);
    }

    public static double acosInRange(double value) {
        if (value <= -1.0) {
            return Math.PI;
        }
        if (value >= 1.0) {
            return 0.0;
        }
        return FastMath.acos(value);
    }

    public static double atan(double value) {
        if (USE_JDK_MATH) {
            return Math.atan(value);
        }
        boolean negateResult = false;
        if (value < 0.0) {
            value = -value;
            negateResult = true;
        }
        if (value == 1.0) {
            return negateResult ? -0.7853981633974483 : 0.7853981633974483;
        }
        if (value <= ATAN_MAX_VALUE_FOR_TABS) {
            int index = (int)(value * ATAN_INDEXER + 0.5);
            double delta = value - (double)index * ATAN_DELTA;
            double result = CmnFastMath.MyTAtan.atanTab[index] + delta * (CmnFastMath.MyTAtan.atanDer1DivF1Tab[index] + delta * (CmnFastMath.MyTAtan.atanDer2DivF2Tab[index] + delta * (CmnFastMath.MyTAtan.atanDer3DivF3Tab[index] + delta * CmnFastMath.MyTAtan.atanDer4DivF4Tab[index])));
            return negateResult ? -result : result;
        }
        if (value < TWO_POW_66) {
            double x = -1.0 / value;
            double x2 = x * x;
            double x4 = x2 * x2;
            double s1 = x2 * (ATAN_AT0 + x4 * (ATAN_AT2 + x4 * (ATAN_AT4 + x4 * (ATAN_AT6 + x4 * (ATAN_AT8 + x4 * ATAN_AT10)))));
            double s2 = x4 * (ATAN_AT1 + x4 * (ATAN_AT3 + x4 * (ATAN_AT5 + x4 * (ATAN_AT7 + x4 * ATAN_AT9))));
            double result = ATAN_HI3 - (x * (s1 + s2) - ATAN_LO3 - x);
            return negateResult ? -result : result;
        }
        if (value != value) {
            return Double.NaN;
        }
        return negateResult ? -1.5707963267948966 : 1.5707963267948966;
    }

    public static double atan2(double y, double x) {
        if (USE_JDK_MATH) {
            return Math.atan2(y, x);
        }
        if (x > 0.0) {
            if (y == 0.0) {
                return y;
            }
            if (x == Double.POSITIVE_INFINITY) {
                return FastMath.atan2_pinf_yyy(y);
            }
            return FastMath.atan(y / x);
        }
        if (x < 0.0) {
            if (y == 0.0) {
                return (double)FastMath.signFromBit(y) * Math.PI;
            }
            if (x == Double.NEGATIVE_INFINITY) {
                return FastMath.atan2_ninf_yyy(y);
            }
            if (y > 0.0) {
                return 1.5707963267948966 - FastMath.atan(x / y);
            }
            if (y < 0.0) {
                return -1.5707963267948966 - FastMath.atan(x / y);
            }
            return Double.NaN;
        }
        return FastMath.atan2_yyy_zeroOrNaN(y, x);
    }

    public static double toRadians(double angdeg) {
        if (USE_JDK_MATH) {
            return Math.toRadians(angdeg);
        }
        return angdeg * (Math.PI / 180);
    }

    public static double toDegrees(double angrad) {
        if (USE_JDK_MATH) {
            return Math.toDegrees(angrad);
        }
        return angrad * 57.29577951308232;
    }

    public static double toRadians(boolean sign, int degrees, int minutes, double seconds) {
        return FastMath.toRadians(FastMath.toDegrees(sign, degrees, minutes, seconds));
    }

    public static double toDegrees(boolean sign, int degrees, int minutes, double seconds) {
        double signFactor = sign ? 1.0 : -1.0;
        return signFactor * ((double)degrees + 0.016666666666666666 * ((double)minutes + 0.016666666666666666 * seconds));
    }

    public static boolean toDMS(double angrad, IntWrapper degrees, IntWrapper minutes, DoubleWrapper seconds) {
        boolean isNeg;
        double tmp = FastMath.toDegrees(FastMath.normalizeMinusPiPi(angrad));
        boolean bl = isNeg = tmp < 0.0;
        if (isNeg) {
            tmp = -tmp;
        }
        degrees.value = (int)tmp;
        tmp = (tmp - (double)degrees.value) * 60.0;
        minutes.value = (int)tmp;
        seconds.value = Math.min((tmp - (double)minutes.value) * 60.0, DOUBLE_BEFORE_60);
        return !isNeg;
    }

    public static boolean isInClockwiseDomain(double startAngRad, double angSpanRad, double angRad) {
        if (Math.abs(angRad) < 2.4492935982947064E-16) {
            if (angSpanRad <= Math.PI * 2) {
                double endAngRad;
                if (angSpanRad < 0.0) {
                    return false;
                }
                if ((startAngRad = FastMath.normalizeMinusPiPi(startAngRad)) <= (endAngRad = FastMath.normalizeMinusPiPi(startAngRad + angSpanRad))) {
                    return angRad >= startAngRad && angRad <= endAngRad;
                }
                return angRad >= startAngRad || angRad <= endAngRad;
            }
            return angSpanRad == angSpanRad;
        }
        return FastMath.normalizeZeroTwoPi(angRad - startAngRad) <= angSpanRad;
    }

    public static double sinh(double value) {
        double h;
        if (USE_JDK_MATH) {
            return Math.sinh(value);
        }
        if (value < 0.0) {
            value = -value;
            h = -0.5;
        } else {
            h = 0.5;
        }
        if (value < 22.0) {
            if (value < TWO_POW_N28) {
                return h < 0.0 ? -value : value;
            }
            double t = FastMath.expm1(value);
            return h * (t + t / (t + 1.0));
        }
        if (value < LOG_DOUBLE_MAX_VALUE) {
            return h * FastMath.exp(value);
        }
        double t = FastMath.exp(value * 0.5);
        return h * t * t;
    }

    public static double cosh(double value) {
        if (USE_JDK_MATH) {
            return Math.cosh(value);
        }
        if (value < 0.0) {
            value = -value;
        }
        if (value < LOG_TWO_POW_27) {
            if (value < TWO_POW_N27) {
                return 1.0;
            }
            double t = FastMath.exp(value);
            return 0.5 * (t + 1.0 / t);
        }
        if (value < LOG_DOUBLE_MAX_VALUE) {
            return 0.5 * FastMath.exp(value);
        }
        double t = FastMath.exp(value * 0.5);
        return 0.5 * t * t;
    }

    public static double coshm1(double value) {
        if (value < 0.0) {
            value = -value;
        }
        if (value < LOG_TWO_POW_27) {
            if (value < TWO_POW_N27) {
                if (value == 0.0) {
                    return value;
                }
                return 0.5 * value * value;
            }
            return 0.5 * (FastMath.expm1(value) + FastMath.expm1(-value));
        }
        if (value < LOG_DOUBLE_MAX_VALUE) {
            return 0.5 * FastMath.exp(value) - 1.0;
        }
        double t = FastMath.exp(value * 0.5);
        return 0.5 * t * t;
    }

    public static double sinhAndCosh(double value, DoubleWrapper hcosine) {
        double hsine;
        double h;
        if (USE_JDK_MATH) {
            hcosine.value = Math.cosh(value);
            return Math.sinh(value);
        }
        if (value < 0.0) {
            value = -value;
            h = -0.5;
        } else {
            h = 0.5;
        }
        if (value < LOG_TWO_POW_27) {
            double t;
            if (value < TWO_POW_N28) {
                hsine = h < 0.0 ? -value : value;
            } else {
                t = FastMath.expm1(value);
                hsine = h * (t + t / (t + 1.0));
            }
            if (value < TWO_POW_N27) {
                hcosine.value = 1.0;
            } else {
                t = FastMath.exp(value);
                hcosine.value = 0.5 * (t + 1.0 / t);
            }
        } else if (value < 22.0) {
            double t = FastMath.expm1(value);
            hsine = h * (t + t / (t + 1.0));
            hcosine.value = 0.5 * (t + 1.0);
        } else {
            if (value < LOG_DOUBLE_MAX_VALUE) {
                hsine = h * FastMath.exp(value);
            } else {
                double t = FastMath.exp(value * 0.5);
                hsine = h * t * t;
            }
            hcosine.value = Math.abs(hsine);
        }
        return hsine;
    }

    public static double tanh(double value) {
        double z;
        if (USE_JDK_MATH) {
            return Math.tanh(value);
        }
        boolean negateResult = false;
        if (value < 0.0) {
            value = -value;
            negateResult = true;
        }
        if (value < 19.061547465398498) {
            if (value < TWO_POW_N55) {
                return negateResult ? -value * (1.0 - value) : value * (1.0 + value);
            }
            if (value >= 1.0) {
                z = 1.0 - 2.0 / (FastMath.expm1(value + value) + 2.0);
            } else {
                double t = FastMath.expm1(-(value + value));
                z = -t / (t + 2.0);
            }
        } else {
            z = value != value ? Double.NaN : 1.0;
        }
        return negateResult ? -z : z;
    }

    public static double asinh(double value) {
        double result;
        boolean negateResult = false;
        if (value < 0.0) {
            value = -value;
            negateResult = true;
        }
        if (value < 0.04) {
            double x = value;
            double x2 = x * x;
            double argLog1p = x * (1.0 + 0.5 * x * (1.0 + -0.25 * x2 * (1.0 + -0.5 * x2 * (1.0 + -0.625 * x2 * (1.0 + -0.7 * x2)))));
            result = FastMath.log1p(argLog1p);
        } else {
            result = value < 1.6777216E7 ? FastMath.log(value + FastMath.sqrt(value * value + 1.0)) : LOG_2 + FastMath.log(value);
        }
        return negateResult ? -result : result;
    }

    public static double acosh(double value) {
        if (!(value > 1.0)) {
            return value < 1.0 ? Double.NaN : value - 1.0;
        }
        double result = value < 1.6777216E7 ? FastMath.log(value + FastMath.sqrt(value * value - 1.0)) : LOG_2 + FastMath.log(value);
        return result;
    }

    public static double acosh1p(double value) {
        if (!(value > 0.0)) {
            return value < 0.0 ? Double.NaN : value;
        }
        double result = value < 1.6777215E7 ? FastMath.log1p(value + FastMath.sqrt(value * (2.0 + value))) : LOG_2 + FastMath.log(1.0 + value);
        return result;
    }

    public static double atanh(double value) {
        boolean negateResult = false;
        if (value < 0.0) {
            value = -value;
            negateResult = true;
        }
        double result = !(value < 1.0) ? (value > 1.0 ? Double.NaN : Double.POSITIVE_INFINITY + value) : 0.5 * FastMath.log1p((value + value) / (1.0 - value));
        return negateResult ? -result : result;
    }

    public static double exp(double value) {
        if (USE_JDK_MATH) {
            return Math.exp(value);
        }
        if (value > EXP_OVERFLOW_LIMIT) {
            return Double.POSITIVE_INFINITY;
        }
        if (!(value >= EXP_UNDERFLOW_LIMIT)) {
            return value != value ? Double.NaN : 0.0;
        }
        int indexes = (int)(value * (double)EXP_LO_INDEXING);
        int valueInt = indexes >= 0 ? indexes >> EXP_LO_INDEXING_DIV_SHIFT : -(-indexes >> EXP_LO_INDEXING_DIV_SHIFT);
        double hiTerm = CmnFastMath.MyTExp.expHiTab[valueInt - (int)EXP_UNDERFLOW_LIMIT];
        int zIndex = indexes - (valueInt << EXP_LO_INDEXING_DIV_SHIFT);
        double y = value - (double)valueInt;
        double z = (double)zIndex * (1.0 / (double)EXP_LO_INDEXING);
        double eps = y - z;
        double expZ = CmnFastMath.MyTExp.expLoPosTab[zIndex + EXP_LO_TAB_MID_INDEX];
        double expEps = 1.0 + eps * (1.0 + eps * (0.5 + eps * (0.16666666666666666 + eps * 0.041666666666666664)));
        double loTerm = expZ * expEps;
        return hiTerm * loTerm;
    }

    public static double expQuick(double value) {
        if (USE_JDK_MATH) {
            return Math.exp(value);
        }
        return Double.longBitsToDouble((long)((int)(1512775.3952 * value + 1.0726481222E9)) << 32);
    }

    public static double expm1(double value) {
        if (USE_JDK_MATH) {
            return Math.expm1(value);
        }
        if (Math.abs(value) < 1.0) {
            int i = (int)(value * (double)EXP_LO_INDEXING);
            double delta = value - (double)i * (1.0 / (double)EXP_LO_INDEXING);
            return CmnFastMath.MyTExp.expLoPosTab[i + EXP_LO_TAB_MID_INDEX] * (CmnFastMath.MyTExp.expLoNegTab[i + EXP_LO_TAB_MID_INDEX] + delta * (1.0 + delta * (0.5 + delta * (0.16666666666666666 + delta * (0.041666666666666664 + delta * 0.008333333333333333)))));
        }
        return FastMath.exp(value) - 1.0;
    }

    public static double log(double value) {
        if (USE_JDK_MATH || !USE_REDEFINED_LOG) {
            return Math.log(value);
        }
        if (value > 0.0) {
            double h;
            if (value == Double.POSITIVE_INFINITY) {
                return Double.POSITIVE_INFINITY;
            }
            if (value > 0.95) {
                if (value < 1.14) {
                    double z = (value - 1.0) / (value + 1.0);
                    double z2 = z * z;
                    return z * (2.0 + z2 * (0.6666666666666666 + z2 * (0.4 + z2 * (0.2857142857142857 + z2 * (0.2222222222222222 + z2 * 0.18181818181818182)))));
                }
                h = 0.0;
            } else if (value < DOUBLE_MIN_NORMAL) {
                value *= TWO_POW_52;
                h = -52.0 * LOG_2;
            } else {
                h = 0.0;
            }
            int valueBitsHi = (int)(Double.doubleToRawLongBits(value) >> 32);
            int valueExp = (valueBitsHi >> 20) - 1023;
            int xIndex = valueBitsHi << 12 >>> 32 - LOG_BITS;
            double z = value * FastMath.twoPowNormalOrSubnormal(-valueExp) * CmnFastMath.MyTLog.logXInvTab[xIndex] - 1.0;
            z *= 1.0 - z * (0.5 - z * 0.3333333333333333);
            return h + (double)valueExp * LOG_2 + (CmnFastMath.MyTLog.logXLogTab[xIndex] + z);
        }
        if (value == 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        return Double.NaN;
    }

    public static double logQuick(double value) {
        double h;
        if (USE_JDK_MATH) {
            return Math.log(value);
        }
        if (value > 0.87) {
            if (value < 1.16) {
                return 2.0 * (value - 1.0) / (value + 1.0);
            }
            h = 0.0;
        } else if (value < DOUBLE_MIN_NORMAL) {
            value *= TWO_POW_52;
            h = -52.0 * LOG_2;
        } else {
            h = 0.0;
        }
        int valueBitsHi = (int)(Double.doubleToRawLongBits(value) >> 32);
        int valueExp = (valueBitsHi >> 20) - 1023;
        int xIndex = valueBitsHi << 12 >>> 32 - LOG_BITS;
        return h + (double)valueExp * LOG_2 + CmnFastMath.MyTLog.logXLogTab[xIndex];
    }

    public static double log10(double value) {
        if (USE_JDK_MATH || !USE_REDEFINED_LOG) {
            return Math.log10(value);
        }
        return FastMath.log(value) * INV_LOG_10;
    }

    public static double log1p(double value) {
        if (USE_JDK_MATH) {
            return Math.log1p(value);
        }
        if (value > -1.0) {
            if (value == Double.POSITIVE_INFINITY) {
                return Double.POSITIVE_INFINITY;
            }
            double valuePlusOne = 1.0 + value;
            if (valuePlusOne == 1.0) {
                return value;
            }
            if (Math.abs(value) < 0.15) {
                double z = value / (value + 2.0);
                double z2 = z * z;
                return z * (2.0 + z2 * (0.6666666666666666 + z2 * (0.4 + z2 * (0.2857142857142857 + z2 * (0.2222222222222222 + z2 * 0.18181818181818182)))));
            }
            int valuePlusOneBitsHi = (int)(Double.doubleToRawLongBits(valuePlusOne) >> 32) & Integer.MAX_VALUE;
            int valuePlusOneExp = (valuePlusOneBitsHi >> 20) - 1023;
            int xIndex = valuePlusOneBitsHi << 12 >>> 32 - LOG_BITS;
            double z = valuePlusOne * FastMath.twoPowNormalOrSubnormal(-valuePlusOneExp) * CmnFastMath.MyTLog.logXInvTab[xIndex] - 1.0;
            z *= 1.0 - z * (0.5 - z * 0.3333333333333333);
            return (double)valuePlusOneExp * LOG_2 + CmnFastMath.MyTLog.logXLogTab[xIndex] + (z + (value - (valuePlusOne - 1.0)) / valuePlusOne);
        }
        if (value == -1.0) {
            return Double.NEGATIVE_INFINITY;
        }
        return Double.NaN;
    }

    public static double pow(double value, double power) {
        if (USE_JDK_MATH) {
            return Math.pow(value, power);
        }
        if (power == 0.0) {
            return 1.0;
        }
        if (power == 1.0) {
            return value;
        }
        if (value <= 0.0) {
            int powerInfo;
            if (Math.abs(power) >= TWO_POW_52 * 2.0) {
                powerInfo = 1;
            } else if (Math.abs(power) <= 2.147483647E9) {
                int powerAsInt = (int)power;
                powerInfo = power == (double)powerAsInt ? ((powerAsInt & 1) == 0 ? 1 : -1) : 0;
            } else {
                long powerAsLong = (long)power;
                if (power == (double)powerAsLong) {
                    powerInfo = (powerAsLong & 1L) == 0L ? 1 : -1;
                } else {
                    if (power != power) {
                        return Double.NaN;
                    }
                    powerInfo = 0;
                }
            }
            if (value == 0.0) {
                if (power < 0.0) {
                    return powerInfo < 0 ? 1.0 / value : Double.POSITIVE_INFINITY;
                }
                return powerInfo < 0 ? value : 0.0;
            }
            if (value == Double.NEGATIVE_INFINITY) {
                if (powerInfo < 0) {
                    return power < 0.0 ? -0.0 : Double.NEGATIVE_INFINITY;
                }
                return power < 0.0 ? 0.0 : Double.POSITIVE_INFINITY;
            }
            return powerInfo == 0 ? Double.NaN : (double)powerInfo * FastMath.exp(power * FastMath.log(-value));
        }
        return FastMath.exp(power * FastMath.log(value));
    }

    public static double powQuick(double value, double power) {
        if (USE_JDK_MATH) {
            return Math.pow(value, power);
        }
        return FastMath.exp(power * FastMath.logQuick(value));
    }

    public static double powFast(double value, int power) {
        if (USE_JDK_MATH) {
            return Math.pow(value, power);
        }
        if (power < 3) {
            if (power < 0) {
                if (power == Integer.MIN_VALUE) {
                    return 1.0 / (FastMath.powFast(value, Integer.MAX_VALUE) * value);
                }
                return 1.0 / FastMath.powFast(value, -power);
            }
            if (power == 2) {
                return value * value;
            }
            if (power == 0) {
                return 1.0;
            }
            return value;
        }
        double oddRemains = 1.0;
        while (power > 5) {
            if ((power & 1) != 0) {
                oddRemains *= value;
            }
            value *= value;
            power >>= 1;
        }
        if (power == 3) {
            return oddRemains * value * value * value;
        }
        double v2 = value * value;
        if (power == 4) {
            return oddRemains * v2 * v2;
        }
        return oddRemains * v2 * v2 * value;
    }

    public static float pow2(float value) {
        return value * value;
    }

    public static double pow2(double value) {
        return value * value;
    }

    public static float pow3(float value) {
        return value * value * value;
    }

    public static double pow3(double value) {
        return value * value * value;
    }

    public static double sqrt(double value) {
        double h;
        if (USE_JDK_MATH || !USE_REDEFINED_SQRT) {
            return Math.sqrt(value);
        }
        if (!(value > 0.0)) {
            return value < 0.0 ? Double.NaN : value;
        }
        if (value == Double.POSITIVE_INFINITY) {
            return Double.POSITIVE_INFINITY;
        }
        if (value < DOUBLE_MIN_NORMAL) {
            value *= TWO_POW_52;
            h = 2.0 * TWO_POW_N26;
        } else {
            h = 2.0;
        }
        int valueBitsHi = (int)(Double.doubleToRawLongBits(value) >> 32);
        int valueExponentIndex = (valueBitsHi >> 20) + 51;
        int xIndex = valueBitsHi << 12 >>> 32 - SQRT_LO_BITS;
        double result = CmnFastMath.MyTSqrt.sqrtXSqrtHiTab[valueExponentIndex] * CmnFastMath.MyTSqrt.sqrtXSqrtLoTab[xIndex];
        double slope = CmnFastMath.MyTSqrt.sqrtSlopeHiTab[valueExponentIndex] * CmnFastMath.MyTSqrt.sqrtSlopeLoTab[xIndex];
        result += ((value *= 0.25) - result * result) * slope;
        result += (value - result * result) * slope;
        return h * (result + (value - result * result) * slope);
    }

    public static double sqrtQuick(double value) {
        if (USE_JDK_MATH) {
            return Math.sqrt(value);
        }
        long bits = Double.doubleToRawLongBits(value);
        return Double.longBitsToDouble(bits + 4606859074900000000L >>> 1);
    }

    public static double invSqrtQuick(double value) {
        if (USE_JDK_MATH) {
            return 1.0 / Math.sqrt(value);
        }
        return Double.longBitsToDouble(6910469410427058089L - (Double.doubleToRawLongBits(value) >> 1));
    }

    public static double cbrt(double value) {
        double h;
        if (USE_JDK_MATH) {
            return Math.cbrt(value);
        }
        if (value < 0.0) {
            if (value == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }
            if ((value = -value) < DOUBLE_MIN_NORMAL) {
                value *= TWO_POW_52 * TWO_POW_26;
                h = -2.0 * TWO_POW_N26;
            } else {
                h = -2.0;
            }
        } else {
            if (!(value < Double.POSITIVE_INFINITY)) {
                return value;
            }
            if (value < DOUBLE_MIN_NORMAL) {
                if (value == 0.0) {
                    return value;
                }
                value *= TWO_POW_52 * TWO_POW_26;
                h = 2.0 * TWO_POW_N26;
            } else {
                h = 2.0;
            }
        }
        int valueBitsHi = (int)(Double.doubleToRawLongBits(value) >> 32);
        int valueExponentIndex = (valueBitsHi >> 20) + 51;
        int xIndex = valueBitsHi << 12 >>> 32 - CBRT_LO_BITS;
        double result = CmnFastMath.MyTCbrt.cbrtXCbrtHiTab[valueExponentIndex] * CmnFastMath.MyTCbrt.cbrtXCbrtLoTab[xIndex];
        double slope = CmnFastMath.MyTCbrt.cbrtSlopeHiTab[valueExponentIndex] * CmnFastMath.MyTCbrt.cbrtSlopeLoTab[xIndex];
        result += ((value *= 0.125) - result * result * result) * slope;
        result += (value - result * result * result) * slope;
        return h * (result + (value - result * result * result) * slope);
    }

    public static double hypot(double x, double y) {
        double factor;
        if (USE_JDK_MATH) {
            return Math.hypot(x, y);
        }
        x = Math.abs(x);
        if ((y = Math.abs(y)) < x) {
            double a = x;
            x = y;
            y = a;
        } else if (!(y >= x)) {
            return FastMath.hypot_NaN(x, y);
        }
        if (y - x == y) {
            return y;
        }
        if (y > HYPOT_MAX_MAG) {
            x *= 1.0 / HYPOT_FACTOR;
            y *= 1.0 / HYPOT_FACTOR;
            factor = HYPOT_FACTOR;
        } else if (x < 1.0 / HYPOT_MAX_MAG) {
            x *= HYPOT_FACTOR;
            y *= HYPOT_FACTOR;
            factor = 1.0 / HYPOT_FACTOR;
        } else {
            factor = 1.0;
        }
        return factor * FastMath.sqrt(x * x + y * y);
    }

    public static double hypot(double x, double y, double z) {
        double factor;
        double a;
        if (USE_JDK_MATH) {
            // empty if block
        }
        x = Math.abs(x);
        y = Math.abs(y);
        if ((z = Math.abs(z)) > y) {
            a = z;
            z = y;
            y = a;
        } else if (!(z <= y)) {
            return FastMath.hypot_NaN(x, y, z);
        }
        if (z > x) {
            double oldZ = z;
            z = x;
            double oldY = y;
            y = oldZ;
            x = oldY;
        } else if (y > x) {
            a = y;
            y = x;
            x = a;
        } else if (x != x) {
            return FastMath.hypot_NaN(x, y, z);
        }
        if (x - y == x) {
            return x;
        }
        if (y - z == y) {
            if (x > HYPOT_MAX_MAG) {
                x *= 1.0 / HYPOT_FACTOR;
                y *= 1.0 / HYPOT_FACTOR;
                factor = HYPOT_FACTOR;
            } else if (y < 1.0 / HYPOT_MAX_MAG) {
                x *= HYPOT_FACTOR;
                y *= HYPOT_FACTOR;
                factor = 1.0 / HYPOT_FACTOR;
            } else {
                factor = 1.0;
            }
            return factor * FastMath.sqrt(x * x + y * y);
        }
        if (x > HYPOT_MAX_MAG) {
            x *= 1.0 / HYPOT_FACTOR;
            y *= 1.0 / HYPOT_FACTOR;
            z *= 1.0 / HYPOT_FACTOR;
            factor = HYPOT_FACTOR;
        } else if (z < 1.0 / HYPOT_MAX_MAG) {
            x *= HYPOT_FACTOR;
            y *= HYPOT_FACTOR;
            z *= HYPOT_FACTOR;
            factor = 1.0 / HYPOT_FACTOR;
        } else {
            factor = 1.0;
        }
        return factor * FastMath.sqrt(x * x + (y * y + z * z));
    }

    public static float floor(float value) {
        int exponent = FastMath.getExponent(value);
        if (exponent < 0) {
            if (value < 0.0f) {
                return -1.0f;
            }
            return 0.0f * value;
        }
        if (exponent < 23) {
            int bits = Float.floatToRawIntBits(value);
            int anteCommaBits = bits & -8388608 >> exponent;
            if (value < 0.0f && anteCommaBits != bits) {
                return Float.intBitsToFloat(anteCommaBits) - 1.0f;
            }
            return Float.intBitsToFloat(anteCommaBits);
        }
        return value;
    }

    public static double floor(double value) {
        if (USE_JDK_MATH) {
            return Math.floor(value);
        }
        double valueAbs = Math.abs(value);
        if (valueAbs <= 2.147483647E9) {
            if (value > 0.0) {
                return (int)value;
            }
            if (value < 0.0) {
                double anteCommaDigits = (int)value;
                if (value != anteCommaDigits) {
                    return anteCommaDigits - 1.0;
                }
                return anteCommaDigits;
            }
            return value;
        }
        if (valueAbs < TWO_POW_52) {
            double highPart = (double)((int)(value * TWO_POW_N26)) * TWO_POW_26;
            if (value > 0.0) {
                return highPart + (double)((int)(value - highPart));
            }
            double anteCommaDigits = highPart + (double)((int)(value - highPart));
            if (value != anteCommaDigits) {
                return anteCommaDigits - 1.0;
            }
            return anteCommaDigits;
        }
        return value;
    }

    public static float ceil(float value) {
        return -FastMath.floor(-value);
    }

    public static double ceil(double value) {
        if (USE_JDK_MATH) {
            return Math.ceil(value);
        }
        return -FastMath.floor(-value);
    }

    public static int round(float value) {
        int bits = Float.floatToRawIntBits(value);
        int biasedExp = bits >> 23 & 0xFF;
        int shift = 149 - biasedExp;
        if ((shift & 0xFFFFFFE0) == 0) {
            int bitsSignum = (bits >> 31 << 1) + 1;
            int extendedMantissa = (0x800000 | bits & 0x7FFFFF) * bitsSignum;
            return (extendedMantissa >> shift) + 1 >> 1;
        }
        return (int)value;
    }

    public static long round(double value) {
        long bits = Double.doubleToRawLongBits(value);
        int biasedExp = (int)(bits >> 52) & 0x7FF;
        int shift = 1074 - biasedExp;
        if ((shift & 0xFFFFFFC0) == 0) {
            long bitsSignum = (bits >> 63 << 1) + 1L;
            long extendedMantissa = (0x10000000000000L | bits & 0xFFFFFFFFFFFFFL) * bitsSignum;
            return (extendedMantissa >> shift) + 1L >> 1;
        }
        if (Math.abs(value) >= 9.223372036854776E18) {
            return value < 0.0 ? Long.MIN_VALUE : Long.MAX_VALUE;
        }
        return (long)value;
    }

    public static int roundEven(float value) {
        int sign = FastMath.signFromBit(value);
        if ((value = Math.abs(value)) < TWO_POW_23_F) {
            value = value + TWO_POW_23_F - TWO_POW_23_F;
            return sign * (int)value;
        }
        if (value < 2.14748365E9f) {
            return sign * (int)value;
        }
        return (int)((float)sign * value);
    }

    public static long roundEven(double value) {
        int sign = (int)FastMath.signFromBit(value);
        if ((value = Math.abs(value)) < TWO_POW_52) {
            value = value + TWO_POW_52 - TWO_POW_52;
        }
        if (value <= 2.147483647E9) {
            return sign * (int)value;
        }
        return (long)((double)sign * value);
    }

    public static float rint(float value) {
        int sign = FastMath.signFromBit(value);
        if ((value = Math.abs(value)) < TWO_POW_23_F) {
            value = TWO_POW_23_F + value - TWO_POW_23_F;
        }
        return (float)sign * value;
    }

    public static double rint(double value) {
        if (USE_JDK_MATH) {
            return Math.rint(value);
        }
        int sign = (int)FastMath.signFromBit(value);
        if ((value = Math.abs(value)) < TWO_POW_52) {
            value = TWO_POW_52 + value - TWO_POW_52;
        }
        return (double)sign * value;
    }

    public static int floorToInt(double value) {
        int valueInt = (int)value;
        if (value < 0.0) {
            if (value == (double)valueInt) {
                return valueInt;
            }
            if (valueInt == Integer.MIN_VALUE) {
                return valueInt;
            }
            return valueInt - 1;
        }
        return valueInt;
    }

    public static int ceilToInt(double value) {
        int valueInt = (int)value;
        if (value > 0.0) {
            if (value == (double)valueInt) {
                return valueInt;
            }
            if (valueInt == Integer.MAX_VALUE) {
                return valueInt;
            }
            return valueInt + 1;
        }
        return valueInt;
    }

    public static int roundToInt(double value) {
        return NumbersUtils.toInt(FastMath.round(value));
    }

    public static int roundEvenToInt(double value) {
        int sign = (int)FastMath.signFromBit(value);
        value = Math.abs(value);
        value = value + TWO_POW_52 - TWO_POW_52;
        return (int)((double)sign * value);
    }

    public static float toRange(float min, float max, float value) {
        return NumbersUtils.toRange(min, max, value);
    }

    public static double toRange(double min, double max, double value) {
        return NumbersUtils.toRange(min, max, value);
    }

    public static double remainder(double dividend, double divisor) {
        if (Double.isInfinite(divisor)) {
            if (Double.isInfinite(dividend)) {
                return Double.NaN;
            }
            return dividend;
        }
        double value = dividend % divisor;
        if (Math.abs(value + value) > Math.abs(divisor)) {
            return value + (value > 0.0 ? -Math.abs(divisor) : Math.abs(divisor));
        }
        return value;
    }

    public static double normalizeMinusPiPi(double angle) {
        if (angle >= -Math.PI && angle <= Math.PI) {
            return angle;
        }
        return FastMath.remainderTwoPi(angle);
    }

    public static double normalizeMinusPiPiFast(double angle) {
        if (angle >= -Math.PI && angle <= Math.PI) {
            return angle;
        }
        return FastMath.remainderTwoPiFast(angle);
    }

    public static double normalizeZeroTwoPi(double angle) {
        if (angle >= 0.0 && angle <= Math.PI * 2) {
            return angle;
        }
        if ((angle = FastMath.remainderTwoPi(angle)) < 0.0) {
            return angle + TWOPI_LO + TWOPI_HI;
        }
        return angle;
    }

    public static double normalizeZeroTwoPiFast(double angle) {
        if (angle >= 0.0 && angle <= Math.PI * 2) {
            return angle;
        }
        if ((angle = FastMath.remainderTwoPiFast(angle)) < 0.0) {
            return angle + TWOPI_LO + TWOPI_HI;
        }
        return angle;
    }

    public static double normalizeMinusHalfPiHalfPi(double angle) {
        if (angle >= -1.5707963267948966 && angle <= 1.5707963267948966) {
            return angle;
        }
        return FastMath.remainderPi(angle);
    }

    public static double normalizeMinusHalfPiHalfPiFast(double angle) {
        if (angle >= -1.5707963267948966 && angle <= 1.5707963267948966) {
            return angle;
        }
        return FastMath.remainderPiFast(angle);
    }

    public static boolean isNaNOrInfinite(float value) {
        return NumbersUtils.isNaNOrInfinite(value);
    }

    public static boolean isNaNOrInfinite(double value) {
        return NumbersUtils.isNaNOrInfinite(value);
    }

    public static int getExponent(float value) {
        return (Float.floatToRawIntBits(value) >> 23 & 0xFF) - 127;
    }

    public static int getExponent(double value) {
        return ((int)(Double.doubleToRawLongBits(value) >> 52) & 0x7FF) - 1023;
    }

    public static float signum(float value) {
        if (USE_JDK_MATH) {
            return Math.signum(value);
        }
        if (value == 0.0f || value != value) {
            return value;
        }
        return FastMath.signFromBit(value);
    }

    public static double signum(double value) {
        if (USE_JDK_MATH) {
            return Math.signum(value);
        }
        if (value == 0.0 || value != value) {
            return value;
        }
        return (int)FastMath.signFromBit(value);
    }

    public static int signFromBit(float value) {
        return Float.floatToRawIntBits(value) >> 30 | 1;
    }

    public static long signFromBit(double value) {
        return Double.doubleToRawLongBits(value) >> 62 | 1L;
    }

    public static float copySign(float magnitude, float sign) {
        return Float.intBitsToFloat(Float.floatToRawIntBits(sign) & Integer.MIN_VALUE | Float.floatToRawIntBits(magnitude) & Integer.MAX_VALUE);
    }

    public static double copySign(double magnitude, double sign) {
        return Double.longBitsToDouble(Double.doubleToRawLongBits(sign) & Long.MIN_VALUE | Double.doubleToRawLongBits(magnitude) & Long.MAX_VALUE);
    }

    public static float ulp(float value) {
        if (USE_JDK_MATH) {
            return Math.ulp(value);
        }
        int exponent = FastMath.getExponent(value);
        if (exponent >= -103) {
            if (exponent == 128) {
                return Math.abs(value);
            }
            return Float.intBitsToFloat(exponent + 104 << 23);
        }
        if (exponent == -127) {
            return Float.MIN_VALUE;
        }
        return Float.intBitsToFloat(1 << exponent - -126);
    }

    public static double ulp(double value) {
        if (USE_JDK_MATH) {
            return Math.ulp(value);
        }
        int exponent = FastMath.getExponent(value);
        if (exponent >= -970) {
            if (exponent == 1024) {
                return Math.abs(value);
            }
            return Double.longBitsToDouble((long)exponent + 971L << 52);
        }
        if (exponent == -1023) {
            return Double.MIN_VALUE;
        }
        return Double.longBitsToDouble(1L << exponent - -1022);
    }

    public static float nextAfter(float start, double direction) {
        if (direction < (double)start) {
            int bits;
            if (start == 0.0f) {
                return -1.4E-45f;
            }
            return Float.intBitsToFloat(bits + ((bits = Float.floatToRawIntBits(start)) > 0 ? -1 : 1));
        }
        if (direction > (double)start) {
            int bits;
            return Float.intBitsToFloat(bits + ((bits = Float.floatToRawIntBits(start + 0.0f)) >= 0 ? 1 : -1));
        }
        if ((double)start == direction) {
            return (float)direction;
        }
        return start + (float)direction;
    }

    public static double nextAfter(double start, double direction) {
        if (direction < start) {
            long bits;
            if (start == 0.0) {
                return -4.9E-324;
            }
            return Double.longBitsToDouble(bits + (long)((bits = Double.doubleToRawLongBits(start)) > 0L ? -1 : 1));
        }
        if (direction > start) {
            long bits;
            return Double.longBitsToDouble(bits + (long)((bits = Double.doubleToRawLongBits(start + 0.0)) >= 0L ? 1 : -1));
        }
        if (start == direction) {
            return direction;
        }
        return start + direction;
    }

    public static float nextDown(float start) {
        if (start > Float.NEGATIVE_INFINITY) {
            int bits;
            if (start == 0.0f) {
                return -1.4E-45f;
            }
            return Float.intBitsToFloat(bits + ((bits = Float.floatToRawIntBits(start)) > 0 ? -1 : 1));
        }
        if (start == Float.NEGATIVE_INFINITY) {
            return Float.NEGATIVE_INFINITY;
        }
        return start;
    }

    public static double nextDown(double start) {
        if (start > Double.NEGATIVE_INFINITY) {
            long bits;
            if (start == 0.0) {
                return -4.9E-324;
            }
            return Double.longBitsToDouble(bits + (long)((bits = Double.doubleToRawLongBits(start)) > 0L ? -1 : 1));
        }
        if (start == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }
        return start;
    }

    public static float nextUp(float start) {
        if (start < Float.POSITIVE_INFINITY) {
            int bits;
            return Float.intBitsToFloat(bits + ((bits = Float.floatToRawIntBits(start + 0.0f)) >= 0 ? 1 : -1));
        }
        if (start == Float.POSITIVE_INFINITY) {
            return Float.POSITIVE_INFINITY;
        }
        return start;
    }

    public static double nextUp(double start) {
        if (start < Double.POSITIVE_INFINITY) {
            long bits;
            return Double.longBitsToDouble(bits + (long)((bits = Double.doubleToRawLongBits(start + 0.0)) >= 0L ? 1 : -1));
        }
        if (start == Double.POSITIVE_INFINITY) {
            return Double.POSITIVE_INFINITY;
        }
        return start;
    }

    public static float scalb(float value, int scaleFactor) {
        int MAX_SCALE = 278;
        scaleFactor = Math.max(Math.min(scaleFactor, 278), -278);
        return (float)((double)value * FastMath.twoPowNormal(scaleFactor));
    }

    public static double scalb(double value, int scaleFactor) {
        double exponentDelta;
        int scaleIncrement;
        if (scaleFactor > -1023 && scaleFactor <= 1023) {
            return value * FastMath.twoPowNormal(scaleFactor);
        }
        int MAX_SCALE = 2099;
        if (scaleFactor < 0) {
            scaleFactor = Math.max(scaleFactor, -2099);
            scaleIncrement = -512;
            exponentDelta = TWO_POW_N512;
        } else {
            scaleFactor = Math.min(scaleFactor, 2099);
            scaleIncrement = 512;
            exponentDelta = TWO_POW_512;
        }
        int t = scaleFactor >> 8 >>> 23;
        int exponentAdjust = (scaleFactor + t & 0x1FF) - t;
        value *= FastMath.twoPowNormal(exponentAdjust);
        scaleFactor -= exponentAdjust;
        while (scaleFactor != 0) {
            value *= exponentDelta;
            scaleFactor -= scaleIncrement;
        }
        return value;
    }

    public static float abs(float a) {
        return Math.abs(a);
    }

    public static double abs(double a) {
        return Math.abs(a);
    }

    public static float min(float a, float b) {
        return Math.min(a, b);
    }

    public static double min(double a, double b) {
        return Math.min(a, b);
    }

    public static float max(float a, float b) {
        return Math.max(a, b);
    }

    public static double max(double a, double b) {
        return Math.max(a, b);
    }

    public static double IEEEremainder(double f1, double f2) {
        return Math.IEEEremainder(f1, f2);
    }

    public static double random() {
        return Math.random();
    }

    private FastMath() {
    }

    private static double remainderTwoPi(double angle) {
        if (USE_JDK_MATH) {
            return FastMath.jdkRemainderTwoPi(angle);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (angle <= 4.0 * NORMALIZE_ANGLE_MAX_MEDIUM_DOUBLE_PIO2) {
            double fn = (int)(angle * TWOPI_INV + 0.5);
            if ((angle = angle - fn * TWOPI_HI - fn * TWOPI_LO) < -Math.PI) {
                angle = angle + TWOPI_HI + TWOPI_LO;
            } else if (angle > Math.PI) {
                angle = angle - TWOPI_HI - TWOPI_LO;
            }
            return negateResult ? -angle : angle;
        }
        if (angle < Double.POSITIVE_INFINITY) {
            angle = FastMath.heavyRemainderTwoPi(angle);
            return negateResult ? -angle : angle;
        }
        return Double.NaN;
    }

    private static double remainderPi(double angle) {
        if (USE_JDK_MATH) {
            return FastMath.jdkRemainderPi(angle);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (angle <= 2.0 * NORMALIZE_ANGLE_MAX_MEDIUM_DOUBLE_PIO2) {
            double fn = (int)(angle * PI_INV + 0.5);
            if ((angle = angle - fn * PI_HI - fn * PI_LO) < -1.5707963267948966) {
                angle = angle + PI_HI + PI_LO;
            } else if (angle > 1.5707963267948966) {
                angle = angle - PI_HI - PI_LO;
            }
            return negateResult ? -angle : angle;
        }
        if (angle < Double.POSITIVE_INFINITY) {
            angle = FastMath.heavyRemainderPi(angle);
            return negateResult ? -angle : angle;
        }
        return Double.NaN;
    }

    private static long remainderPiO2(double angle) {
        if (USE_JDK_MATH) {
            return FastMath.jdkRemainderPiO2(angle, false);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (angle <= NORMALIZE_ANGLE_MAX_MEDIUM_DOUBLE_PIO2) {
            int n = (int)(angle * PIO2_INV + 0.5);
            double fn = n;
            if ((angle = angle - fn * PIO2_HI - fn * PIO2_LO) < -0.7853981633974483) {
                angle = angle + PIO2_HI + PIO2_LO;
                --n;
            } else if (angle > 0.7853981633974483) {
                angle = angle - PIO2_HI - PIO2_LO;
                ++n;
            }
            if (negateResult) {
                angle = -angle;
            }
            return FastMath.encodeRemainderAndQuadrant(angle, n & 3);
        }
        if (angle < Double.POSITIVE_INFINITY) {
            return FastMath.heavyRemainderPiO2(angle, negateResult);
        }
        return FastMath.encodeRemainderAndQuadrant(Double.NaN, 0);
    }

    private static double remainderTwoPiFast(double angle) {
        double fn;
        if (USE_JDK_MATH) {
            return FastMath.jdkRemainderTwoPi(angle);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (!(angle <= TWO_POW_26 * (Math.PI * 2))) {
            if (angle <= TWO_POW_52 * (Math.PI * 2)) {
                fn = (int)(angle * (TWOPI_INV / TWO_POW_26) + 0.5);
                if ((angle = angle - fn * (TWOPI_HI * TWO_POW_26) - fn * (TWOPI_LO * TWO_POW_26)) < 0.0) {
                    angle = -angle;
                    negateResult = !negateResult;
                }
            } else {
                if (angle < Double.POSITIVE_INFINITY) {
                    return 0.0;
                }
                return Double.NaN;
            }
        }
        if ((angle = angle - (fn = (double)((int)(angle * TWOPI_INV + 0.5))) * TWOPI_HI - fn * TWOPI_LO) < -Math.PI) {
            angle = angle + TWOPI_HI + TWOPI_LO;
        } else if (angle > Math.PI) {
            angle = angle - TWOPI_HI - TWOPI_LO;
        }
        return negateResult ? -angle : angle;
    }

    private static double remainderPiFast(double angle) {
        double fn;
        if (USE_JDK_MATH) {
            return FastMath.jdkRemainderPi(angle);
        }
        boolean negateResult = false;
        if (angle < 0.0) {
            angle = -angle;
            negateResult = true;
        }
        if (!(angle <= TWO_POW_26 * Math.PI)) {
            if (angle <= TWO_POW_52 * Math.PI) {
                fn = (int)(angle * (PI_INV / TWO_POW_26) + 0.5);
                if ((angle = angle - fn * (PI_HI * TWO_POW_26) - fn * (PI_LO * TWO_POW_26)) < 0.0) {
                    angle = -angle;
                    negateResult = !negateResult;
                }
            } else {
                if (angle < Double.POSITIVE_INFINITY) {
                    return 0.0;
                }
                return Double.NaN;
            }
        }
        if ((angle = angle - (fn = (double)((int)(angle * PI_INV + 0.5))) * PI_HI - fn * PI_LO) < -1.5707963267948966) {
            angle = angle + PI_HI + PI_LO;
        } else if (angle > 1.5707963267948966) {
            angle = angle - PI_HI - PI_LO;
        }
        return negateResult ? -angle : angle;
    }
}
