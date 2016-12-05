package plugin.trackmate.tram.trackanalyzer;


import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 *
 * This is a smoothing interpolation function using cubic splines. Core functionality
 * comes from Apache Commons Math.
 *
 * Written by Dan Ruderman (ruderman@usc.edu)
 * <p>
 * Date: 11/23/16
 * Time: 6:30 PM
 * <p>
 * Copyright University of Southern California.
 * All Rights Reserved.
 */
public class Smoother {
    private final PolynomialSplineFunction polynomialSplineFunction;

    /**
	 *
	 * @param xs
	 *            array of x locations
	 * @param ys
	 *            array of y values (same length as xs)
	 * @param numKnots
	 *            number of knots to use in spline (must be &lt;= length of xs
	 *            and &gt;= 4)
	 */
    public Smoother(final double[] xs, final double[] ys, final int numKnots) {
        if ( xs.length != ys.length ) throw new IllegalArgumentException("xs and ys must have the same length");
		if ( numKnots > xs.length )	throw new IllegalArgumentException( "numKnots must not be greater than length of xs" );
        if ( numKnots < 4 ) throw new IllegalArgumentException("must have at least 4 knots");

        // One knot always at each end because the polynomial spline function can't extrapolate
        // figure out which indexes to use as knots. Try to get close to both ends.
        final int spacing = (xs.length-1) / (numKnots-1);
        final int offset = (xs.length - 1 - (numKnots-1)*spacing) / 2; // try to be centered

        // copy over the data we will use
        final double[] xsToUse = new double[numKnots];
        final double[] ysToUse = new double[numKnots];
        // first and last are fixed
        xsToUse[0] = xs[0];
        ysToUse[0] = ys[0];
        xsToUse[numKnots-1] = xs[xs.length-1];
        ysToUse[numKnots-1] = ys[ys.length-1];
        // now do the ones in the middle
        for ( int i=1, idx=spacing+offset; i < numKnots-1; ++i, idx += spacing ) {
            xsToUse[i] = xs[idx];
            ysToUse[i] = ys[idx];
        }

        // create and store the interpolator
        polynomialSplineFunction = (new SplineInterpolator()).interpolate(xsToUse, ysToUse);
    }

    public double value(final double x) {
        return polynomialSplineFunction.value(x);
    }

    /**
     * Test
     * @param args
     */
    public static void main(final String[] args) {
        final double[] xs = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        final double[] ys = new double[] {  -112,  -54,  -20,   -4 ,   0 ,  -2,   -4,    0,   16,   50};

        final Smoother s = new Smoother(xs, ys, 4);

        for ( double x=1; x <= 10; x += 0.5 ) {
            System.out.println(x + " " + s.value(x));
        }
    }
}
