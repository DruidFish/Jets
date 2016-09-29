#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>

#include "TLorentzVector.h"

using namespace std;

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
	vector< double > phis, rapidities, pts;
	while ( cin >> px >> py >> pz >> E )
	{
		dummy1.SetPxPyPzE( px, py, pz, E );
		phis.push_back( dummy1.Phi() );
		rapidities.push_back( dummy1.Rapidity() );
		pts.push_back( dummy1.Pt() );
	}

	//kt^2n = pt^2n //for single objects
	//kt^2n = min( pt^2n ) x delta_R^2 / D^2 //for pairs
	//Just do akt 0.6
	double const MODE = -1; //1 for kt, -1 for akt, 0 for c/a
	double const D = 0.6;
	double const aD2 = 1.0 / ( D * D );

	//Loop over all objects
	while ( phis.size() )
	{
		double kt2Min = DBL_MAX;
		unsigned int thisMinIndex = 0;
		unsigned int pairMinIndex = 0;
		unsigned int const totalObjects = phis.size();
		for ( unsigned int thisObjectIndex = 0; thisObjectIndex < totalObjects; thisObjectIndex++ )
		{
			//Single-object values
			double const thisApt = 1.0 / pts[ thisObjectIndex ];
			double const thisApt2 = thisApt * thisApt;
			double const thisPhi = phis[ thisObjectIndex ];
			double const thisRapidity = rapidities[ thisObjectIndex ];

			//Update minimum kt^2 value
			if ( thisApt2 < kt2Min )
			{
				kt2Min = thisApt2;
				thisMinIndex = thisObjectIndex;
				pairMinIndex = thisObjectIndex;
			}

			//Create pairs
			for ( unsigned int pairObjectIndex = thisObjectIndex + 1; pairObjectIndex < totalObjects; pairObjectIndex++ )
			{
				//Pair values
				double const pairApt = 1.0 / pts[ pairObjectIndex ];
				double const pairApt2 = pairApt * pairApt;
				double const deltaPhi = thisPhi - phis[ pairObjectIndex ];
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

				//Update minimum kt^2 value
				if ( kt2 < kt2Min )
				{
					kt2Min = kt2;
					thisMinIndex = thisObjectIndex;
					pairMinIndex = pairObjectIndex;
				}
			}
		}

		//Single objects as min kt are outputs, pairs get merged
		if ( thisMinIndex == pairMinIndex )
		{
			//note rapidity != eta...
			outputs.push_back( TLorentzVector() );
			outputs.back().SetPtEtaPhiM( pts[ thisMinIndex ], rapidities[ thisMinIndex ], phis[ thisMinIndex ], 0.0 );

			phis.erase( phis.begin() + thisMinIndex );
			rapidities.erase( rapidities.begin() + thisMinIndex );
			pts.erase( pts.begin() + thisMinIndex );

		}
		else
		{
			//note rapidity != eta...
			dummy1.SetPtEtaPhiM( pts[ thisMinIndex ], rapidities[ thisMinIndex ], phis[ thisMinIndex ], 0.0 );
			dummy2.SetPtEtaPhiM( pts[ pairMinIndex ], rapidities[ pairMinIndex ], phis[ pairMinIndex ], 0.0 );
			dummy1 += dummy2;
			pts[ thisMinIndex ] = dummy1.Pt();
			rapidities[ thisMinIndex ] = dummy1.Rapidity();
			phis[ thisMinIndex ] = dummy1.Phi();
			phis.erase( phis.begin() + pairMinIndex );
			rapidities.erase( rapidities.begin() + pairMinIndex );
			pts.erase( pts.begin() + pairMinIndex );
		}
	}

	sort( outputs.begin(), outputs.end(), SortJetsByPt );

	//Output exactly like fastjet demo
	printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
	for ( unsigned int jetIndex = 0; jetIndex < outputs.size(); jetIndex++ )
	{
		if ( outputs[ jetIndex ].Pt() < 5.0 ) break;
		printf( "%5u %15.8f %15.8f %15.8f\n",
				jetIndex,
				outputs[ jetIndex ].Rapidity(),
				outputs[ jetIndex ].Phi(),
				outputs[ jetIndex ].Pt() );
	}

	return 0;
}
