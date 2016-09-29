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
	vector< TLorentzVector > inputs, outputs;

	//Read the fastjet example input
	double px, py, pz, E;
	while ( cin >> px >> py >> pz >> E )
	{
		inputs.push_back( TLorentzVector() );
		inputs.back().SetPxPyPzE( px, py, pz, E );
	}

	//kt^2n = pt^2n //for single objects
	//kt^2n = min( pt^2n ) x delta_R^2 / D^2 //for pairs
	//Just do akt 0.6
	//double const MODE = -1; //1 for kt, -1 for akt, 0 for c/a
	double const D = 0.6;
	double const aD2 = 1.0 / ( D * D );

	//Loop over all objects
	tick_count const startTime = tick_count::now();
	while ( inputs.size() )
	{
		double kt2Min = DBL_MAX;
		unsigned int thisMinIndex = 0;
		unsigned int pairMinIndex = 0;
		unsigned int const totalObjects = inputs.size();

		//Loop over all input objects
		for ( unsigned int thisObjectIndex = 0; thisObjectIndex < totalObjects; thisObjectIndex++ )
		{
			//Single-object values
			TLorentzVector * thisObject = &inputs[ thisObjectIndex ];
			double const thisApt = 1.0 / thisObject->Pt();
			double const thisApt2 = thisApt * thisApt;
			double const thisPhi = thisObject->Phi();
			double const thisRapidity = thisObject->Rapidity();

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
				TLorentzVector * pairObject = &inputs[ pairObjectIndex ];
				double const pairApt = 1.0 / pairObject->Pt();
				double const pairApt2 = pairApt * pairApt;
				double const deltaPhi = thisPhi - pairObject->Phi();
				double const deltaRapidity = thisRapidity - pairObject->Rapidity();
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
			outputs.push_back( inputs[ thisMinIndex ] );
			inputs.erase( inputs.begin() + thisMinIndex );
		}
		else
		{
			inputs[ thisMinIndex ] += inputs[ pairMinIndex ];
			inputs.erase( inputs.begin() + pairMinIndex );
		}
	}
	tick_count const endTime = tick_count::now();
	cout << "Jet finding time: " << ( endTime - startTime ).seconds() << " sec" << endl;

	sort( outputs.begin(), outputs.end(), SortJetsByPt );

	//Output exactly like fastjet demo
	printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
	for ( unsigned int jetIndex = 0; jetIndex < outputs.size(); jetIndex++ )
	{
		if ( outputs[ jetIndex ].Pt() < 5.0 ) break;
		double phi = outputs[ jetIndex ].Phi();
		while ( phi < 0.0 ) phi += ( 2.0 * M_PI );
		printf( "%5u %15.8f %15.8f %15.8f\n",
				jetIndex,
				outputs[ jetIndex ].Rapidity(),
				phi,
				outputs[ jetIndex ].Pt() );
	}

	return 0;
}
