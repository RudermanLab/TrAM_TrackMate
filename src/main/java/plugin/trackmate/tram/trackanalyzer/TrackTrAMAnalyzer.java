package plugin.trackmate.tram.trackanalyzer;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.ImageIcon;

import org.scijava.plugin.Plugin;

import fiji.plugin.trackmate.Dimension;
import fiji.plugin.trackmate.FeatureModel;
import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzerFactory;
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory;
import fiji.plugin.trackmate.features.track.TrackAnalyzer;

@Plugin( type = TrackAnalyzer.class )
public class TrackTrAMAnalyzer implements TrackAnalyzer
{

	private static final String KEY = "TRACK_TRAM_ANALYZER";

	public static final String TRAM = "TRAM";

	private static final List< String > FEATURES = new ArrayList< String >( 1 );

	private static final Map< String, Boolean > IS_INT = new HashMap< String, Boolean >( 1 );

	private static final Map< String, String > FEATURE_SHORT_NAMES = new HashMap< String, String >( 1 );

	private static final Map< String, String > FEATURE_NAMES = new HashMap< String, String >( 1 );

	private static final Map< String, Dimension > FEATURE_DIMENSIONS = new HashMap< String, Dimension >( 1 );

	/**
	 * Unmodifiable set of spot feature keys that can be used by TrAM.
	 */
	static final Set< String > DEFAULT_SELECTED_SPOT_FEATURES;
	static
	{
		final HashSet< String > sf = new HashSet<>();
		sf.add( Spot.POSITION_X );
		sf.add( Spot.POSITION_Y );
		sf.add( Spot.POSITION_Z );
		sf.add( Spot.RADIUS );
		sf.add( SpotIntensityAnalyzerFactory.MEDIAN_INTENSITY );
		sf.add( SpotContrastAndSNRAnalyzerFactory.SNR );
		DEFAULT_SELECTED_SPOT_FEATURES = Collections.unmodifiableSet( sf );
	}

	// todo: get this via interface
	private static final Map< String, String[] > euclidianMap = new HashMap<>();

	static
	{
		FEATURES.add( TRAM );

		IS_INT.put( TRAM, false );

		FEATURE_NAMES.put( TRAM, "TrAM" );

		FEATURE_SHORT_NAMES.put( TRAM, "TrAM" );

		FEATURE_DIMENSIONS.put( TRAM, Dimension.NONE );

		// TODO: get this via gui interface.
		euclidianMap.put( "XY", new String[] { "POSITION_X", "POSITION_Y" } );
	}

	static final int DEFAULT_NUM_KNOTS = 7;

	static final double DEFAULT_P = 0.5;

	private double p = DEFAULT_P;

	private int numKnots = DEFAULT_NUM_KNOTS;

	private Set< String > selectedSpotFeatures = new HashSet<>( DEFAULT_SELECTED_SPOT_FEATURES );

	private long processingTime;

	/*
	 * TrAM specific methods.
	 */

	public void setExponent( final double p )
	{
		this.p = p;
	}

	public void setNumKnots( final int numKnots )
	{
		this.numKnots = numKnots;
	}

	public void setSelectedSpotFeatures( final Collection< String > selectedSpotFeatures )
	{
		this.selectedSpotFeatures = new HashSet<>( selectedSpotFeatures );
	}

	/*
	 * TRACKANALYZER METHODS
	 */

	@Override
	public void process( final Collection< Integer > trackIDs, final Model model )
	{
		final long start = System.currentTimeMillis();

		// The feature model where we store the feature values:
		final FeatureModel fm = model.getFeatureModel();

		// First order all the tracks
		final HashMap< Integer, List< Spot > > orderedTracks = new HashMap< Integer, List< Spot > >();
		for ( final Integer trackID : trackIDs )
		{
			// The tracks are indexed by their ID. Here is how to get their
			// content:
			final Set< Spot > spots = model.getTrackModel().trackSpots( trackID );
			// Or .trackEdges( trackID ) if you want the edges.

			// This set is NOT ordered. If we want the first one and last one we
			// have to sort them:
			final Comparator< Spot > comparator = Spot.frameComparator;
			final List< Spot > sorted = new ArrayList<>( spots );
			Collections.sort( sorted, comparator );

			orderedTracks.put( trackID, sorted ); // save
		}


		// perform intersection. These are the features we want and have
		final List< String > spotFeatures = new ArrayList< String >( model.getFeatureModel().getSpotFeatures() );
		spotFeatures.retainAll( selectedSpotFeatures );

		// Now convert to the format needed by the computation. This is an array
		// of Map<String, double[]>
		final int numTracks = trackIDs.size();
		final List< Map< String, double[] > > allTrackDataList = new ArrayList<>( numTracks );
		final Integer[] orderedTrackIDs = trackIDs.toArray( new Integer[ 0 ] );
		for ( final Integer trackID : orderedTrackIDs )
		{
			// do this by order of the orderedTrackIDs
			final List< Spot > spots = orderedTracks.get( trackID );
			final int numTimePoints = spots.size();

			// pull all of the time series data into arrays
			final double[][] dataArray = new double[ spotFeatures.size() ][ numTimePoints ];
			for ( int i = 0; i < numTimePoints; ++i )
			{
				final Spot spot = spots.get( i );

				for ( int featureNumber = 0; featureNumber < spotFeatures.size(); ++featureNumber )
					dataArray[ featureNumber ][ i ] = spot.getFeature( spotFeatures.get( featureNumber ) );
			}

			// store in the map
			final Map< String, double[] > trackData = new HashMap<>( spotFeatures.size() );
			for ( int featureNumber = 0; featureNumber < spotFeatures.size(); ++featureNumber )
				trackData.put( spotFeatures.get( featureNumber ), dataArray[ featureNumber ] );

			// Collect with all the other tracks.
			allTrackDataList.add( trackData );
		}

		// compute the overall track statistics
		final Map< String, Double > medianAbsolutedifferences = TrAM.computeMedianAbsoluteDifferences( allTrackDataList );

		final TrAM tram = new TrAM( medianAbsolutedifferences, numKnots, p );

		final Set< String > availableFeatures = tram.getAvailableFeatures();

		// now for each track compute its TrAM value
		for ( int i = 0; i < allTrackDataList.size(); ++i )
		{
			final Integer trackID = orderedTrackIDs[ i ];
			// Get the track data.
			final Map< String, double[] > trackData = allTrackDataList.get( i );
			final HashSet< String > unavailableFeatures = new HashSet< String >( trackData.keySet() );
			// Remove any data without statistics.
			unavailableFeatures.removeAll( availableFeatures );
			for ( final String featureToRemove : unavailableFeatures )
				trackData.remove( featureToRemove );

			final double tramValue = tram.computeTrAM( trackData, euclidianMap );

			fm.putTrackFeature( trackID, TRAM, tramValue ); // store it
		}

		final long end = System.currentTimeMillis();
		this.processingTime = end - start;
	}

	@Override
	public boolean isLocal()
	{
		return false;
	} // not local because we amass stats across all tracks


	@Override
	public Map< String, Boolean > getIsIntFeature()
	{
		return Collections.unmodifiableMap( IS_INT );
	}

	@Override
	public boolean isManualFeature()
	{
		return true;
	}

	/*
	 * TRACKMODULE METHODS
	 */

	@Override
	public String getKey()
	{
		return KEY;
	}

	@Override
	public String getName()
	{
		return "Track aberration measure";
	}

	@Override
	public String getInfoText()
	{
		return "";
	}

	@Override
	public ImageIcon getIcon()
	{
		return null;
	}

	/*
	 * FEATUREANALYZER METHODS
	 */

	@Override
	public List< String > getFeatures()
	{
		return FEATURES;
	}

	@Override
	public Map< String, String > getFeatureShortNames()
	{
		return FEATURE_SHORT_NAMES;
	}

	@Override
	public Map< String, String > getFeatureNames()
	{
		return FEATURE_NAMES;
	}

	@Override
	public Map< String, Dimension > getFeatureDimensions()
	{
		return FEATURE_DIMENSIONS;
	}

	/*
	 * MULTITHREADED METHODS
	 */

	@Override
	public void setNumThreads()
	{}

	@Override
	public void setNumThreads( final int numThreads )
	{}

	@Override
	public int getNumThreads()
	{
		return 1;
	}

	/*
	 * BENCHMARK METHOD
	 */

	@Override
	public long getProcessingTime()
	{
		return processingTime;
	}
}
