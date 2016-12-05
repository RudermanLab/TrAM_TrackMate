package plugin.trackmate.tram.trackanalyzer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.collections.primitives.ArrayDoubleList;
import org.apache.commons.math3.stat.descriptive.rank.Median;

/**
 * Written by Dan Ruderman (ruderman@usc.edu)
 * <p>
 * Date: 11/23/16 Time: 9:24 PM
 * <p>
 *
 *
 *     Implementation of Tracking Aberration Measure (TrAM) from:
 *     Patsch, et al. (2016). Single cell dynamic phenotyping. Scientific Reports, 6(October),
 *     34785. http://doi.org/10.1038/srep34785
 *
 *
 * Copyright (c) 2016 Dan Ruderman
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
 * Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */
public class TrAM
{
	/**
	 * Holds median absolute difference by feature name.
	 */
	private final Map< String, Double > medianAbsoluteDifferenceByFeature;

	private final int numKnots;

	private final double p;

	/**
	 * Instantiates a TrAM calculator.
	 *
	 * @param medianAbsoluteDifferenceByFeature
	 *            statistics of trajectory fluctuations (@see
	 *            computeMedianAbsoluteDifferences).
	 * @param numKnots
	 *            mumber of knots to use in time series smoothing.
	 * @param p
	 *            TrAM exponent (typically &lt; 1).
	 */
	public TrAM( final Map< String, Double > medianAbsoluteDifferenceByFeature, final int numKnots, final double p )
	{
		this.medianAbsoluteDifferenceByFeature = new HashMap< String, Double >( medianAbsoluteDifferenceByFeature );
		this.numKnots = numKnots;
		this.p = p;

		// remove any features with zero variability
		for ( final String feature : this.medianAbsoluteDifferenceByFeature.keySet().toArray( new String[ 0 ] ) )
		{
			if ( this.medianAbsoluteDifferenceByFeature.get( feature ) == 0 )
			{
				this.medianAbsoluteDifferenceByFeature.remove( feature );
			}
		}
	}

	/**
	 * Get a set of feature names for which we have fluctuation statistics.
	 * These can be used in the TrAM computation.
	 *
	 * @return a set of strings.
	 */
	public Set< String > getAvailableFeatures()
	{
		return medianAbsoluteDifferenceByFeature.keySet();
	}

	/**
	 * Computes the TrAM statistics for a set of feature values in a single
	 * track.
	 *
	 * @param timeSeriesByFeature
	 *            map containing a double[] time series for each tracked
	 *            feature.
	 * @param euclidianFeatures
	 *            map containing String[] of features which are to be handled in
	 *            Euclidian fashion. May be null.
	 * @return the TrAM statistic
	 */
	public double computeTrAM( final Map< String, double[] > timeSeriesByFeature, final Map< String, String[] > euclidianFeatures )
	{
		// length of data vectors (just use the first vector found and trust the
		// rest are all the same)
		final int timeSeriesLength = timeSeriesByFeature.get( timeSeriesByFeature.keySet().iterator().next() ).length;
		if ( timeSeriesLength < numKnots )
			return Double.NaN;

		// compute the absolute difference of each time series from its smoothed
		// version, normalized by median
		final Map< String, double[] > absDiffsByFeature = new HashMap<>();
		for ( final String key : timeSeriesByFeature.keySet() )
		{
			final double scaleFactor = medianAbsoluteDifferenceByFeature.get( key );

			final double[] vals = timeSeriesByFeature.get( key );
			final double[] times = IntStream.range( 0, vals.length ).asDoubleStream().toArray();


			final Smoother s = new Smoother( times, vals, numKnots );

			final double[] absDiff = new double[ times.length ];
			for ( int i = 0; i < times.length; ++i )
				absDiff[ i ] = Math.abs( vals[ i ] - s.value( times[ i ] ) ) / scaleFactor;

			absDiffsByFeature.put( key, absDiff );
		}

		// figure out the weightings, which depends on how many Euclidian sets
		// there are
		final Map< String, Double > weights = new HashMap<>();
		if ( euclidianFeatures == null || euclidianFeatures.size() == 0 )
		{
			// easy case - all weights the same

			final Set< String > keys = timeSeriesByFeature.keySet();

			for ( final String key : keys )
			{
				weights.put( key, 1.0 / keys.size() );
			}
		}
		else
		{
			// trickier: go through each euclidian list and compute its
			// Euclidian summary of the difference
			// Give it the right relative weighting, which is the number of
			// components the summary is made from.
			// We will be replacing absDiffsByFeature elements by their
			// Euclidian summaries

			final Set< String > originalFeaturesRemaining = timeSeriesByFeature.keySet();
			final int numOriginalFeatures = originalFeaturesRemaining.size();

			for ( final String euclidianKey : euclidianFeatures.keySet() )
			{
				final String[] featuresToCombine = euclidianFeatures.get( euclidianKey );
				final int numFeatures = featuresToCombine.length;

				// get all the time series we will combine as euclidian
				final double[][] euclidianData = new double[ numFeatures ][];
				for ( int idx = 0; idx < numFeatures; ++idx )
				{
					euclidianData[ idx ] = absDiffsByFeature.remove( featuresToCombine[ idx ] );
					// keep track of remaining keys.
					originalFeaturesRemaining.remove( featuresToCombine[ idx ] );
				}

				// perform the euclidian combination
				final double[] euclidianizedData = new double[ timeSeriesLength ];
				for ( int idx = 0; idx < euclidianizedData.length; ++idx )
				{
					double sumSquare = 0;
					for ( int i = 0; i < numFeatures; i++ )
					{
						sumSquare += Math.pow( euclidianData[ i ][ idx ], 2 );
					}
					euclidianizedData[ idx ] = Math.sqrt( sumSquare / numFeatures );
				}

				// now put in the Euclidianized version
				absDiffsByFeature.put( euclidianKey, euclidianizedData );
				// and store the relative weight
				weights.put( euclidianKey, Double.valueOf( numFeatures ) );
			}

			// Set weights for any/all remaining features
			for ( final String originalFeature : originalFeaturesRemaining )
				weights.put( originalFeature, 1.0 );

			// and divide each weight by the total
			for ( final String key : weights.keySet() )
				weights.put( key, weights.get( key ) / numOriginalFeatures );
		}

		// sanity check: make sure we have the same keys (probably shouldn't be
		// in production code)
		if ( !weights.keySet().equals( absDiffsByFeature.keySet() ) )
			throw new IllegalStateException( "Weight/data key mismatch" );

		// Now compute the TrAM statistic for each time point.
		// First extract all the data in an array format
		final String[] keys = weights.keySet().toArray( new String[ 0 ] );
		final double[] weightsArray = new double[ keys.length ];
		final double[][] absDiffArray = new double[ keys.length ][];
		for ( int idx = 0; idx < keys.length; ++idx )
		{
			weightsArray[ idx ] = weights.get( keys[ idx ] );
			absDiffArray[ idx ] = absDiffsByFeature.get( keys[ idx ] );
		}
		// then do the TrAM computation
		final double[] stat = new double[ timeSeriesLength ];
		for ( int t = 0; t < timeSeriesLength; ++t )
		{
			double sum = 0;
			for ( int idx = 0; idx < absDiffArray.length; ++idx )
				sum += weightsArray[ idx ] * Math.pow( absDiffArray[ idx ][ t ], p );
			stat[ t ] = Math.pow( sum, 1.0 / p );
		}

		/*
		 * TrAM is the maximum value. Note we could have computed this along the
		 * way instead of saving all such values, but we want the option to
		 * extract more info later (like the time point of maximal value).
		 */
		final double TrAM = Arrays.stream( stat ).max().getAsDouble();

		return TrAM;
	}

	/**
	 * Convenience method in the absence of Euclidian features.
	 *
	 * @param timeSeriesByFeature
	 *            map containing a double[] time series for each tracked
	 *            feature.
	 * @return the TrAM statistics.
	 */
	public double computeTrAM( final Map< String, double[] > timeSeriesByFeature )
	{
		return computeTrAM( timeSeriesByFeature, null );
	}

	/**
	 * Computes the median absolute difference between adjacent values within
	 * each data vector.
	 * 
	 * @param timeSeries
	 *            a staggered array whose rows are each independent time series.
	 * @return the median absolute difference between adjacent values within all
	 *         time series.
	 */
	public static double computeMedianAbsoluteDifference( final double[][] timeSeries )
	{
		if ( timeSeries.length < 1 )
			throw new IllegalArgumentException( "Must have at least one data vector." );

		// where we accumulate all the absolute differences
		final ArrayDoubleList accum = new ArrayDoubleList();

		// compute absolute differences for each time series
		for ( final double[] vals : timeSeries )
		{
			for ( int i = 0; i < vals.length - 1; ++i )
				accum.add( Math.abs( vals[ i + 1 ] - vals[ i ] ) );
		}

		final Median m = new Median();
		m.setData( accum.toArray() );

		return m.evaluate();
	}

	/**
	 * Computes the median absolute differences between adjacent time series
	 * points across all features and all tracks. Returns a map from feature
	 * name to the median absolute difference value for that feature.
	 *
	 * @param timeSeriesForAllTracks
	 *            list of maps, one map per track, with each map containing time
	 *            series named by feature.
	 * @return median absolute differences for each feature name across all
	 *         tracks.
	 */
	public static Map< String, Double > computeMedianAbsoluteDifferences( final List< Map< String, double[] > > timeSeriesForAllTracks )
	{
		// we have to pull the time series together by feature name
		// store them here
		final Map< String, List< double[] > > timeSeriesByFeature = new HashMap<>();

		// loop through all tracks, and accumulate their time series
		for ( final Map< String, double[] > timeSeriesForTrack : timeSeriesForAllTracks )
		{
			for ( final String key : timeSeriesForTrack.keySet() )
			{
				List< double[] > accum = timeSeriesByFeature.get( key );

				if ( accum == null )
				{ // first one - need to create
					accum = new ArrayList< double[] >();
					timeSeriesByFeature.put( key, accum );
				}

				accum.add( timeSeriesForTrack.get( key ) ); // add it in
			}
		}

		// now do the computation for each feature
		final Map< String, Double > resultMap = new HashMap<>();
		for ( final String key : timeSeriesByFeature.keySet() )
		{
			final double[][] allTimeSeries = timeSeriesByFeature.get( key ).toArray( new double[ 0 ][] );

			resultMap.put( key, computeMedianAbsoluteDifference( allTimeSeries ) );
		}

		return resultMap;
	}

	/**
	 * Unit test
	 */
	public static void main( final String[] args )
	{
		final int numTracks = 500;
		final int numDiscontinuousTracks = 10;
		final int numTimePoints = 50;
		final int discontinuitySafeZoneWidth = 5;
		final String[] features = { "x", "y", "f1", "f2", "f3" };
		final double[] featureFluctuationScales = { 1, 1, 5, 10, 0.1 };
		final double relativeDiscontinuityScale = 9;

		final int numKnots = 6;

		// populate with Gaussian noise. Add discontinuities for the first ones
		final Random rnd = new Random();
		final List< Map< String, double[] > > allTrackData = new ArrayList<>();

		for ( int i = 0; i < numTracks; ++i )
		{
			final Map< String, double[] > featureData = new HashMap<>( features.length );
			allTrackData.add( featureData );

			final boolean isDiscontinuous = i < numDiscontinuousTracks;
			final int discontinuityTimePoint = discontinuitySafeZoneWidth + rnd.nextInt( numTimePoints - 2 * discontinuitySafeZoneWidth );

			for ( int j = 0; j < features.length; ++j )
			{
				final double[] data = new double[ numTimePoints ];
				for ( int k = 0; k < numTimePoints; ++k )
					data[ k ] = featureFluctuationScales[ j ] * rnd.nextGaussian();

				// if discontinuous add a big jump
				if ( isDiscontinuous )
					data[ discontinuityTimePoint ] += relativeDiscontinuityScale *
							featureFluctuationScales[ j ] * rnd.nextGaussian();

				featureData.put( features[ j ], data ); // save it
			}
		}

		final Map< String, Double > mad = computeMedianAbsoluteDifferences( allTrackData );
		System.out.println( mad ); // sanity check

		final Map< String, String[] > euclidians = new HashMap<>();
		euclidians.put( "XY", new String[] { "x", "y" } );
		final TrAM tram = new TrAM( mad, numKnots, 0.5 );
		final double[] tramValues = new double[ numTracks ];
		for ( int trackNumber = 0; trackNumber < numTracks; ++trackNumber )
		{
			tramValues[ trackNumber ] = tram.computeTrAM( allTrackData.get( trackNumber ), euclidians );

			// print out some of the discontinuous and an equal number of
			// continuous track TrAM values
			if ( trackNumber < 2 * numDiscontinuousTracks )
				System.out.println( trackNumber + " " + tramValues[ trackNumber ] );
		}
	}
}
