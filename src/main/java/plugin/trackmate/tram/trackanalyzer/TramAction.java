package plugin.trackmate.tram.trackanalyzer;

import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import org.scijava.plugin.Plugin;

import fiji.plugin.trackmate.FeatureModel;
import fiji.plugin.trackmate.LoadTrackMatePlugIn_;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.action.AbstractTMAction;
import fiji.plugin.trackmate.action.TrackMateAction;
import fiji.plugin.trackmate.action.TrackMateActionFactory;
import fiji.plugin.trackmate.gui.TrackMateGUIController;
import fiji.plugin.trackmate.gui.TrackMateWizard;
import ij.ImageJ;
import ij.gui.GenericDialog;

public class TramAction extends AbstractTMAction
{

	public static final ImageIcon ICON = new ImageIcon( TrackMateWizard.class.getResource( "images/arrow_merge.png" ) );

	public static final String NAME = "Tracking Aberration Measure (TrAM)";

	public static final String KEY = "TRAM";

	public static final String INFO_TEXT = "<html>"
			+ "<b>TrAM</b>: Tracking Aberration Measure."
			+ "<p>"
			+ "See http://www.nature.com/articles/srep34785"
			+ "</html>";

	private final TrackMateGUIController controller;

	private static final List< String > PREVIOUSLY_SELECTED_FEATURES = new ArrayList<>( TrackTrAMAnalyzer.DEFAULT_SELECTED_SPOT_FEATURES );

	private static double PREVIOUS_P = TrackTrAMAnalyzer.DEFAULT_P;

	private static int PREVIOUS_NKNOTS = TrackTrAMAnalyzer.DEFAULT_NUM_KNOTS;

	public TramAction( final TrackMateGUIController controller )
	{
		this.controller = controller;
	}

	@Override
	public void execute( final TrackMate trackmate )
	{
		final GenericDialog dialog = new GenericDialog( NAME, controller.getGUI() );
		dialog.setFont( TrackMateWizard.FONT );
		dialog.addMessage( "Enter parameters for the TrAM calculation.", TrackMateWizard.BIG_FONT );

		// Feature model.
		final FeatureModel fm = trackmate.getModel().getFeatureModel();

		// Exponent and N knots.
		dialog.addNumericField( "TrAM exponent", PREVIOUS_P, 2 );
		dialog.addNumericField( "N knots", PREVIOUS_NKNOTS, 0 );

		// Spot features table.
		final List< String > spotFeatures = new ArrayList<>( fm.getSpotFeatures() );
		final String[] labels = new String[ spotFeatures.size() ];
		final boolean[] defaultValues = new boolean[ spotFeatures.size() ];
		int i = 0;
		for ( final String sf : spotFeatures )
		{
			labels[ i ] = fm.getSpotFeatureNames().get( sf );
			defaultValues[ i ] = PREVIOUSLY_SELECTED_FEATURES.contains( sf );
			i++;
		}
		final int columns = 3;
		final int rows = labels.length / columns + 1;
		dialog.addCheckboxGroup( rows, columns, labels, defaultValues );

		dialog.showDialog();
		if ( !dialog.wasOKed() )
			return;

		final double p = dialog.getNextNumber();
		final int nKnots = ( int ) dialog.getNextNumber();
		final List<String> selectedFeatures = new ArrayList<>();
		for ( final String sf : spotFeatures )
		{
			if (dialog.getNextBoolean())
				selectedFeatures.add( sf );
		}

		// For subsequent execution.
		PREVIOUSLY_SELECTED_FEATURES.clear();
		PREVIOUSLY_SELECTED_FEATURES.addAll( selectedFeatures );
		PREVIOUS_NKNOTS = nKnots;
		PREVIOUS_P = p;

		final TrackTrAMAnalyzer analyzer = new TrackTrAMAnalyzer();
		analyzer.setExponent( p );
		analyzer.setNumKnots( nKnots );
		analyzer.setSelectedSpotFeatures( selectedFeatures );
		
		logger.log( String.format( "Starting TrAM calculation with p=%.2f, N knots = %d.\n", p, nKnots ) );
		logger.log( "Spot features used for calculation:\n" );
		for ( final String sf : selectedFeatures )
			logger.log( " - " + fm.getSpotFeatureNames().get( sf ) + "\n" );

		// Declare track features.
		fm.declareTrackFeatures(
				analyzer.getFeatures(),
				analyzer.getFeatureNames(),
				analyzer.getFeatureShortNames(),
				analyzer.getFeatureDimensions(),
				analyzer.getIsIntFeature() );

		// Compute feature values.
		analyzer.process( trackmate.getModel().getTrackModel().trackIDs( true ), trackmate.getModel() );

		logger.log( "Done in " + analyzer.getProcessingTime() + " ms." );
	}

	@Plugin( type = TrackMateActionFactory.class, visible = true )
	public static class Factory implements TrackMateActionFactory
	{

		@Override
		public String getInfoText()
		{
			return INFO_TEXT;
		}

		@Override
		public String getName()
		{
			return NAME;
		}

		@Override
		public String getKey()
		{
			return KEY;
		}

		@Override
		public ImageIcon getIcon()
		{
			return ICON;
		}


		@Override
		public TrackMateAction create( final TrackMateGUIController controller )
		{
			return new TramAction( controller );
		}

	}


	/*
	 * MAIN METHOD
	 */

	public static void main( final String[] args )
	{
		ImageJ.main( args );
		new LoadTrackMatePlugIn_().run( "../TrackMate/samples/FakeTracks.xml" );
	}
}
