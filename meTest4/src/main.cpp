#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>

#include "tbb/tick_count.h"

#include "TLorentzVector.h"

using namespace std;
using namespace tbb;

bool SortJetsByPt( TLorentzVector const& i, TLorentzVector const& j )
{
	return ( i.Pt() > j.Pt() );
}

//Try to make some jets myself
int main( int argc, char * argv[] )
{
	vector< TLorentzVector > outputs;

	//Dummy TLV for conversions
	TLorentzVector dummy1;
	TLorentzVector dummy2;

	//Read the fastjet example input
	double px, py, pz, E;
	vector< double > phis, rapidities, pts, energies, etas;
	while ( cin >> px >> py >> pz >> E )
	{
		dummy1.SetPxPyPzE( px, py, pz, E );
		phis.push_back( dummy1.Phi() );
		rapidities.push_back( dummy1.Rapidity() );
		pts.push_back( dummy1.Pt() );
		energies.push_back( E );
		etas.push_back( dummy1.Eta() );
	}
	vector< double > cachedMinKts( phis.size(), 0.0 );
	vector< unsigned int > cachedMinKtPairs( phis.size(), 0 );

	//kt^2n = pt^2n //for single objects
	//kt^2n = min( pt^2n ) x delta_R^2 / D^2 //for pairs
	//Just do akt 0.6
	//double const MODE = -1; //1 for kt, -1 for akt, 0 for c/a
	double const D = 0.6;
	double const aD2 = 1.0 / ( D * D );
	double const TWO_PI = 2.0 * M_PI;

	//Make the jets
	tick_count const startTime = tick_count::now();
	tick_count::interval_t findMinKtTime;
	tick_count::interval_t updateCollectionsTime;
	unsigned int thisMinIndex = 0;
	unsigned int pairMinIndex = 0;
	unsigned int totalObjects = phis.size();
	while( totalObjects )
	{
		//Loop over all input objects
		tick_count const startKtTime = tick_count::now();
		for ( unsigned int thisObjectIndex = 0; thisObjectIndex < totalObjects; thisObjectIndex++ )
		{
			//Update cache if the indicated object was changed
			unsigned int const cachedPair = cachedMinKtPairs[ thisObjectIndex ];
			if ( cachedPair == thisMinIndex || cachedPair == pairMinIndex )
			{
				//Self values
				double const thisApt = 1.0 / pts[ thisObjectIndex ];
				double const thisApt2 = thisApt * thisApt;
				double const thisPhi = phis[ thisObjectIndex ];
				double const thisRapidity = rapidities[ thisObjectIndex ];
				double minKt2 = thisApt2;
				unsigned int minKt2Pair = thisObjectIndex;

				//Create pairs
				for ( unsigned int pairObjectIndex = thisObjectIndex + 1; pairObjectIndex < totalObjects; pairObjectIndex++ )
				{
					//Pair values
					double const pairApt = 1.0 / pts[ pairObjectIndex ];
					double const pairApt2 = pairApt * pairApt;
					double deltaPhi = fabs( thisPhi - phis[ pairObjectIndex ] );
					if ( deltaPhi > M_PI ) deltaPhi -= TWO_PI;
					double const deltaRapidity = thisRapidity - rapidities[ pairObjectIndex ];
					double kt2 = ( ( deltaPhi * deltaPhi ) + ( deltaRapidity * deltaRapidity ) ) * aD2;
					if ( pairApt2 < thisApt2 )
					{
						kt2 *= pairApt2;
					}
					else
					{
						kt2 *= thisApt2;
					}

					//Find minimum
					if ( kt2 < minKt2 )
					{
						minKt2 = kt2;
						minKt2Pair = pairObjectIndex;
					}
				}

				//Store minimum kt^2
				cachedMinKts[ thisObjectIndex ] = minKt2;
				cachedMinKtPairs[ thisObjectIndex ] = minKt2Pair;
			}
		}

		//Find min kt^2 value over all objects
		thisMinIndex = 0;
		double overallMinKt2 = cachedMinKts[ 0 ];
		for ( unsigned int thisObjectIndex = 1; thisObjectIndex < totalObjects; thisObjectIndex++ )
		{
			double cachedKt2 = cachedMinKts[ thisObjectIndex ];
			if ( cachedKt2 < overallMinKt2 )
			{
				overallMinKt2 = cachedKt2;
				thisMinIndex = thisObjectIndex;
			}
		}
		pairMinIndex = cachedMinKtPairs[ thisMinIndex ];
		findMinKtTime += tick_count::now() - startKtTime;

		//Single objects as min kt are outputs, pairs get merged
		tick_count const startUpdateTime = tick_count::now();
		if ( thisMinIndex == pairMinIndex )
		{
			//Make output TLV
			outputs.push_back( TLorentzVector() );
			outputs.back().SetPtEtaPhiE( pts[ thisMinIndex ], etas[ thisMinIndex ], phis[ thisMinIndex ], energies[ thisMinIndex ] );

			//Remove jet from active data
			phis.erase( phis.begin() + thisMinIndex );
			rapidities.erase( rapidities.begin() + thisMinIndex );
			pts.erase( pts.begin() + thisMinIndex );
			etas.erase( etas.begin() + thisMinIndex );
			energies.erase( energies.begin() + thisMinIndex );
			cachedMinKts.erase( cachedMinKts.begin() + thisMinIndex );
			cachedMinKtPairs.erase( cachedMinKtPairs.begin() + thisMinIndex );

			//Decrement cached indices > erased position
			totalObjects = phis.size();
			for ( unsigned int thisObjectIndex = 0; thisObjectIndex < totalObjects; thisObjectIndex++ )
			{
				if ( cachedMinKtPairs[ thisObjectIndex ] > thisMinIndex )
				{
					cachedMinKtPairs[ thisObjectIndex ]--;
				}
			}
		}
		else
		{
			//Make TLVs for convenience
			dummy1.SetPtEtaPhiE( pts[ thisMinIndex ], etas[ thisMinIndex ], phis[ thisMinIndex ], energies[ thisMinIndex ] );
			dummy2.SetPtEtaPhiE( pts[ pairMinIndex ], etas[ pairMinIndex ], phis[ pairMinIndex ], energies[ pairMinIndex ] );
			dummy1 += dummy2;

			//Replace the original object with the merge
			pts[ thisMinIndex ] = dummy1.Pt();
			rapidities[ thisMinIndex ] = dummy1.Rapidity();
			phis[ thisMinIndex ] = dummy1.Phi();
			etas[ thisMinIndex ] = dummy1.Eta();
			energies[ thisMinIndex ] = dummy1.E();
			cachedMinKts[ thisMinIndex ] = 0.0;
			cachedMinKtPairs[ thisMinIndex ] = thisMinIndex; //Marked to update

			//Remove the pair object
			phis.erase( phis.begin() + pairMinIndex );
			rapidities.erase( rapidities.begin() + pairMinIndex );
			pts.erase( pts.begin() + pairMinIndex );
			etas.erase( etas.begin() + pairMinIndex );
			energies.erase( energies.begin() + pairMinIndex );
			cachedMinKts.erase( cachedMinKts.begin() + pairMinIndex );
			cachedMinKtPairs.erase( cachedMinKtPairs.begin() + pairMinIndex );

			//Decrement cached indices > erased position
			totalObjects = phis.size();
			for ( unsigned int thisObjectIndex = 0; thisObjectIndex < totalObjects; thisObjectIndex++ )
			{
				if ( cachedMinKtPairs[ thisObjectIndex ] > pairMinIndex )
				{
					cachedMinKtPairs[ thisObjectIndex ]--;
				}
			}
		}
		updateCollectionsTime += tick_count::now() - startUpdateTime;
	}
	cout << "Total time: " << ( tick_count::now() - startTime ).seconds() << " sec" << endl;
	cout << "Kt finding time: " << findMinKtTime.seconds() << " sec" << endl;
	cout << "Collection update time: " << updateCollectionsTime.seconds() << " sec" << endl;

	sort( outputs.begin(), outputs.end(), SortJetsByPt );

	//Output exactly like fastjet demo
	printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
	for ( unsigned int jetIndex = 0; jetIndex < outputs.size(); jetIndex++ )
	{
		if ( outputs[ jetIndex ].Pt() < 5.0 ) break;
		double phi = outputs[ jetIndex ].Phi();
		while ( phi < 0.0 ) phi += TWO_PI;
		printf( "%5u %15.8f %15.8f %15.8f\n",
				jetIndex,
				outputs[ jetIndex ].Rapidity(),
				phi,
				outputs[ jetIndex ].Pt() );
	}

	return 0;
}
