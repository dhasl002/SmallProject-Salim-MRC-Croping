#include <iostream>
#include "PDB.h"
#include "MRC.h"
#include <chrono>
#include <ctime>
#include <math.h>
#include "SSTime.h"

void readNCompareMRC();
void readMRCNFindVoxelsPercentAboveThreahold();

int main(int argc, const char * argv[]) {

    SSTime totalTimeCalculator;
    SSTime localTimer;

    PDB pdb;
    MRC mrc;

    totalTimeCalculator.startTimer();
    localTimer.startTimer();

    //Provide Input PDB File Path
    string pdbFilePath ="/Users/salimsazzad/Documents/CS-826-BioGraph/Project/Experimental_Data_Over_6_Angstrom/3izq_9.5/3izq.pdb";
    pdb.readPDBFile(pdbFilePath);
    cout<<"Time to read PDB file: "<< localTimer.getTimeElapsedSinceStart() <<" seconds \n";

    // Output File Path (Retrieved from input pdb path)
    size_t found = pdbFilePath.find_last_of(".");
    string outputFilePath =  pdbFilePath.substr(0,found) + "_output/";


    localTimer.startTimer();

    //Provide Input MRC File Path
    string mrcFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/Experimental_Data_Over_6_Angstrom/3izq_9.5/3izq_emd_1811.map";
    mrc.EMRead(mrcFilePath);
    cout<<"Time to read MRC file: "<< localTimer.getTimeElapsedSinceStart() <<" seconds \n";

    int numOfPDBChains = pdb.getNumOfChains();

    printf("num of PDB chains: %d\n", numOfPDBChains);
    printf("\nChopping of MRC file based on chain starts..........\n");

    for (int i = 0; i < numOfPDBChains; i++) {

        localTimer.startTimer();
        printf("\n\n#Chain:%d -- Start Chopping MRC file of that chain\n", i);

        string mrcOutputFilePath = outputFilePath  + pdb.getChainNameCorrespondingToNum(i)  + ".mrc";

        BoundingCoordinate coordinates = pdb.findEachAxisBoundaryPointsbasedOnChain(i);
        vector<Coordinate> atomsCoordinate = pdb.getAtomsCoordinateOfAChain(i);

        mrc.boundingCubeBasedOnChain( coordinates.minPoints.xCor, coordinates.maxPoints.xCor, coordinates.minPoints.yCor, coordinates.maxPoints.yCor, coordinates.minPoints.zCor, coordinates.maxPoints.zCor, mrcOutputFilePath, atomsCoordinate);

        printf("\n#Chain:%d Chopping of MRC Completed\n", i);
        cout<<"\nTime to chop and write MRC chain : "<< i << " = "<< localTimer.getTimeElapsedSinceStart() <<" seconds \n";
    }

    printf("Code execution complete");
    cout<<"Total time taken: "<< totalTimeCalculator.getTimeElapsedSinceStart() <<" seconds \n";

    return 0;
}


void readNCompareMRC()
{
    /*
     1HZ4:
     Number of NONZero Voxels 517
     Number of Matching Voxels 421
     Matching percentage 0.814313

     1WAB:
     Number of NONZero Voxels 441
     Number of Matching Voxels 316
     Matching percentage 0.716553


     1XQO:
     Number of NONZero Voxels 417
     Number of Matching Voxels 313
     Matching percentage 0.7506



     1Z1L:
     Number of NONZero Voxels 632
     Number of Matching Voxels 385
     Matching percentage 0.609177


     1733_H:
     Number of NONZero Voxels 140
     Number of Matching Voxels 140
     Matching percentage 1

     */



    MRC mrcThikness;
    MRC mrcTensor;

   // string mrcThiknessFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1HZ4_thr_0.38_outFiles/1HZ4_thickHLX.mrc";
     //string mrcThiknessFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1WAB_thr_0.4_outFiles/1WAB_thickHLX.mrc";
    // string mrcThiknessFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1XQO_thr_0.39_outFiles/1XQO_thickHLX.mrc";

   //  string mrcThiknessFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1Z1L_thr_0.34_outFiles/1Z1L_thickHLX.mrc";


    string mrcThiknessFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1740_A_thr_2.3_outFiles/1740_A_thickHLX.mrc";


   // string mrcThiknessFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1733_H_thr_5_outFiles/1733_H_thickHLX.mrc";

    mrcThikness.EMRead(mrcThiknessFilePath);

    ///string mrcTensorFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1HZ4_thr_0.38_outFiles///1HZ4_tensorHLX.mrc";


    // string mrcTensorFilePath = "//Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1WAB_thr_0.4_outFiles/1WAB_tensorHLX.mrc";

   // string mrcTensorFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1XQO_thr_0.39_outFiles/1XQO_tensorHLX.mrc"

    //  string mrcTensorFilePath = "//Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1Z1L_thr_0.34_outFiles/1Z1L_tensorHLX.mrc";

   //   string mrcTensorFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1733_H_thr_5_outFiles/1733_H_tensorHLX.mrc";

    string mrcTensorFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer/SSETracer/resource/1740_A_thr_2.3_outFiles/1740_A_tensorHLX.mrc";
    mrcTensor.EMRead(mrcTensorFilePath);



    int countMatchingVoxels = 0;

    int countNonZeroDensityVoxels = 0;

    for(int i = 0; i < mrcTensor.getNz(); i++) {
        for(int j = 0; j < mrcTensor.getNy(); j++) {
            for(int k = 0; k < mrcTensor.getNx(); k++) {


               // countNonZeroDensityVoxels ++;


                if (mrcTensor.cube[k][j][i] > 0 ) {

                    countNonZeroDensityVoxels ++;
                }

                if (mrcTensor.cube[k][j][i] > 0  &&    fabs(mrcThikness.cube[k][j][i] - mrcTensor.cube[k][j][i]) < 0.01) {

                    countMatchingVoxels++;
                }

            }
        }
    }

    float matchingPercentage = (float)countMatchingVoxels / countNonZeroDensityVoxels * 100;


    cout<<"Number of NONZero Voxels " << countNonZeroDensityVoxels <<"\n" ;
    cout<<"Number of Matching Voxels " << countMatchingVoxels<<"\n";

    cout<<"Matching percentage " << matchingPercentage << "\n";
}


// 1HZ4 =  0.412; Percentage: 2.83031
// 1WAB =  0.448  Percentage: 0.718611
// 1XQO =  0.45   Percentage: 1.0362
// 1Z1L =  0.438; Percentage: 1.66234




void readMRCNFindVoxelsPercentAboveThreahold()
{
    MRC mrcMap;
    string mrcFilePath = "/Users/salimsazzad/Documents/CS-826-BioGraph/Project/SSETracer\ Edited/SSETracer/resource/3iy3_output/A.mrc";
    mrcMap.EMRead(mrcFilePath);

    float threshhold = 2.93;
    int numberofNonZeroVoxels = 0;
    int numberofVoxelsAboveThreshold = 0;

    for(int i = 0; i < mrcMap.getNz(); i++) {
        for(int j = 0; j < mrcMap.getNy(); j++) {
            for(int k = 0; k < mrcMap.getNx(); k++) {
                if (mrcMap.cube[k][j][i] > 0) {
                    numberofNonZeroVoxels++;
                    if( mrcMap.cube[k][j][i] > threshhold){

                        numberofVoxelsAboveThreshold ++;

                    }

                }
            }
        }
    }

    float percentage = (float) numberofVoxelsAboveThreshold / numberofNonZeroVoxels * 100;

    cout<<"Number of NonZero Voxels: "<< numberofNonZeroVoxels<<"\n";
    cout<<"Number of NonZero Voxels: "<< numberofVoxelsAboveThreshold<<"\n";
    cout<<"Percentage: "<< percentage;


}
