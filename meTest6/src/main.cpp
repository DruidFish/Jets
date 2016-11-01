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
bool SortJetsByPtBack( TLorentzVector const& i, TLorentzVector const& j )
{
	return ( i.Pt() < j.Pt() );
}

//Try to make some jets myself
int main( int argc, char * argv[] )
{
	vector< TLorentzVector > inputs, outputs;

	//Dummy TLV for conversions
	TLorentzVector dummy1;
	TLorentzVector dummy2;

	//Read the fastjet example input into TLVs
	double px, py, pz, E;
	vector< double > phisIn, rapiditiesIn, ptsIn, energiesIn, etasIn;
	while ( cin >> px >> py >> pz >> E )
	{
		dummy1.SetPxPyPzE( px, py, pz, E );
		inputs.push_back( dummy1 );
	}

	//Sorting the input pT high to low gives a large speedup
	sort( inputs.begin(), inputs.end(), SortJetsByPt );

	//Sorting low to high is slightly worse than unsorted - perhaps other way around with kt algo?
	sort( inputs.begin(), inputs.end(), SortJetsByPtBack );

	//Copy input data into flat arrays
	unsigned int const totalObjects = inputs.size();
	double phis[ totalObjects ];
	double rapidities[ totalObjects ];
	double pts[ totalObjects ];
	double energies[ totalObjects ];
	double etas[ totalObjects ];
	double cachedMinKts[ totalObjects ];
	unsigned int cachedMinKtPairs[ totalObjects ];
	bool wasDeleted[ totalObjects ];
	for ( unsigned int i = 0; i < totalObjects; i++ )
	{
		TLorentzVector const& dummy = inputs[ i ];
		phis[ i ] = dummy.Phi();
		rapidities[ i ] =  dummy.Rapidity();
		pts[ i ] = dummy.Pt();
		energies[ i ] = dummy.E();
		etas[ i ] = dummy.Eta();
		cachedMinKts[ i ] = 0.0;
		cachedMinKtPairs[ i ] = 0;
		wasDeleted[ i ] = false;
	}

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
	unsigned int activeObjects = totalObjects;
	while( activeObjects )
	{
		//Loop over all input objects
		tick_count const startKtTime = tick_count::now();
		for ( unsigned int thisObjectIndex = 0; thisObjectIndex < totalObjects; thisObjectIndex++ )
		{
			if ( wasDeleted[ thisObjectIndex ] ) continue;

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
					if ( wasDeleted[ pairObjectIndex ] ) continue;

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
			wasDeleted[ thisMinIndex ] = true;
			cachedMinKts[ thisMinIndex ] = DBL_MAX;
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
			wasDeleted[ pairMinIndex ] = true;
			cachedMinKts[ pairMinIndex ] = DBL_MAX;
		}
		activeObjects--;
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
		if ( outputs[ jetIndex ].Pt() < 0.0 ) break;
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
