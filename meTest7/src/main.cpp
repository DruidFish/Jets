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

	//Read the fastjet example input into TLVs
	double px, py, pz, E;
	while ( cin >> px >> py >> pz >> E )
	{
		inputs.push_back( TLorentzVector( px, py, pz, E ) );
	}

	//Sorting the input pT high to low gives a large speedup
	sort( inputs.begin(), inputs.end(), SortJetsByPt );

	//Sorting low to high is slightly worse than unsorted - perhaps other way around with kt algo?
	//sort( inputs.begin(), inputs.end(), SortJetsByPtBack );

	//Copy input data into flat arrays
	unsigned int const totalObjects = inputs.size();
	double phis[ totalObjects ];
	double rapidities[ totalObjects ];
	double pts[ totalObjects ];
	double energies[ totalObjects ];
	double pxs[ totalObjects ];
	double pys[ totalObjects ];
	double pzs[ totalObjects ];
	double cachedMinKts[ totalObjects ];
	unsigned int cachedMinKtPairs[ totalObjects ];
	bool wasDeleted[ totalObjects ];
	for ( unsigned int i = 0; i < totalObjects; i++ )
	{
		TLorentzVector const& dummy = inputs[ i ];
		phis[ i ] = dummy.Phi();
		rapidities[ i ] =  dummy.Rapidity();
		pts[ i ] = dummy.Pt();
		pxs[ i ] = dummy.Px();
		pys[ i ] = dummy.Py();
		pzs[ i ] = dummy.Pz();
		energies[ i ] = dummy.E();
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
	unsigned int firstActive = 0;
	unsigned int lastActive = totalObjects;
	while( activeObjects )
	{
		//Loop over all input objects
		tick_count const startKtTime = tick_count::now();
		for ( unsigned int thisObjectIndex = firstActive; thisObjectIndex < lastActive; thisObjectIndex++ )
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
				for ( unsigned int pairObjectIndex = thisObjectIndex + 1; pairObjectIndex < lastActive; pairObjectIndex++ )
				{
					if ( wasDeleted[ pairObjectIndex ] ) continue;

					//Pair values
					double const pairApt = 1.0 / pts[ pairObjectIndex ];
					double const pairApt2 = pairApt * pairApt;
					double deltaPhi = thisPhi - phis[ pairObjectIndex ];
					if ( deltaPhi > M_PI ) deltaPhi -= TWO_PI;
					if ( deltaPhi < -M_PI ) deltaPhi += TWO_PI;
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
		thisMinIndex = firstActive;
		double overallMinKt2 = cachedMinKts[ firstActive ];
		for ( unsigned int thisObjectIndex = firstActive + 1; thisObjectIndex < lastActive; thisObjectIndex++ )
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
			outputs.push_back( TLorentzVector( pxs[ thisMinIndex ], pys[ thisMinIndex ], pzs[ thisMinIndex ], energies[ thisMinIndex ] ) );

			//Remove jet from active data
			wasDeleted[ thisMinIndex ] = true;
			cachedMinKts[ thisMinIndex ] = DBL_MAX;
		}
		else
		{
			//Merge the pair
			double const newPx = pxs[ thisMinIndex ] + pxs[ pairMinIndex ];
			double const newPy = pys[ thisMinIndex ] + pys[ pairMinIndex ];
			double const newPz = pzs[ thisMinIndex ] + pzs[ pairMinIndex ];
			double const newE = energies[ thisMinIndex ] + energies[ pairMinIndex ];

			//Replace the original object with the merge
			dummy1.SetPxPyPzE( newPx, newPy, newPz, newE ); //for convenience
			pts[ thisMinIndex ] = dummy1.Pt();
			pxs[ thisMinIndex ] = newPx;
			pys[ thisMinIndex ] = newPy;
			pzs[ thisMinIndex ] = newPz;
			rapidities[ thisMinIndex ] = dummy1.Rapidity();
			phis[ thisMinIndex ] = dummy1.Phi();
			energies[ thisMinIndex ] = newE;
			cachedMinKts[ thisMinIndex ] = 0.0;
			cachedMinKtPairs[ thisMinIndex ] = thisMinIndex; //Marked to update

			//Remove the pair object
			wasDeleted[ pairMinIndex ] = true;
			cachedMinKts[ pairMinIndex ] = DBL_MAX;
		}

		//Reduce the search ranges
		for ( unsigned int thisObjectIndex = firstActive; thisObjectIndex < lastActive; thisObjectIndex++ )
		{
			if ( wasDeleted[ thisObjectIndex ] ) firstActive = thisObjectIndex;
			else break;
		}
		for ( unsigned int thisObjectIndex = lastActive; thisObjectIndex > firstActive; thisObjectIndex-- )
		{
			if ( wasDeleted[ thisObjectIndex ] ) lastActive = thisObjectIndex;
			else break;
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
		//if ( outputs[ jetIndex ].Pt() < 5.0 ) break;
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
