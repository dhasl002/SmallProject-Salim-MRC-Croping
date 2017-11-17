#include "PDB.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <stdlib.h>
#include <sys/stat.h>
using namespace std;

#pragma mark -- File Read/Write
void PDB::readPDBFile(string fileName)
{
    printf("Reading PDB file- %s\n", fileName.c_str());

    string line;
    fstream inputFile(fileName.c_str(), ios::in);
    if(!inputFile)
    {
        cout<<"Can not open binary input file "<<fileName<<endl;
        exit(1);
    }

    ATOM atom;

    while (getline(inputFile, line )) {
        getAtoms(line);
    }

    size_t found = fileName.find_last_of(".");
    string outputFilePath = fileName.substr(0,found) + "_output/";

    writePDBAccordingToChain(outputFilePath);
}

void PDB::writeToFileInPDBFormat(string fileName)
{
    fstream output(fileName.c_str(), ios::out);
    if(!output)
    {
        cout<<"Can not open output file "<<fileName<<endl;
        cerr << "Error: " << strerror(errno);
        exit(1);
    }
}

void PDB::writePDBAccordingToChain(string filePath)
{
    printf("\nWritting PDB file based on chain.............\n");
    printf("\nnumber of chains: %d\n", numOfChains);
    printf("PDB File Output Directory: %s\n", filePath.c_str());
    mkdir(filePath.c_str(),0775);

    for (int i = 0; i < numOfChains; i++) {

        string fileName = filePath + getChainNameCorrespondingToNum(i) + ".pdb";

        fstream output(fileName.c_str(), ios::out);
        if(!output)
        {
            cout<<"Can not open output file "<<fileName<<endl;
            cerr << "Error: " << strerror(errno);
            exit(1);
        }

        size_t textSize = chainsPDBText[i].size();
        for (int j = 0; j < textSize; j++) {

            output<<chainsPDBText[i][j]<<"\n";
        }
    }

    printf("\nPDF file writting based on chain complete\n\n");
}

void PDB::getAtoms(string line)
{
    if ((line.substr(0,4) == "ATOM")) {

        numOfAtoms++;

        string chainName = line.substr(21,1);

        char chainName0 = line.at(21);

        string xPos = line.substr(30,8);
        string yPos = line.substr(38,8);
        string zPos = line.substr(46,8);

        Coordinate coordinate;
        coordinate.xCor =  stof(xPos);
        coordinate.yCor  = stof(yPos);
        coordinate.zCor =  stof(zPos);

        ATOM atom;
        atom.setValue(coordinate, chainName);
        atoms.push_back(atom);

        insertPDBTextToCorrespondingChain( chainName0,  line);
        insertAtomToCorrespondingChain(chainName0, coordinate);
    }

    else if ((line.substr(0,5) == "HELIX")) {
        char chainName = line.at(19);
        insertPDBTextToCorrespondingChain( chainName,  line);
    }

    else if ((line.substr(0,5) == "SHEET")) {
        char chainName = line.at(21);
        insertPDBTextToCorrespondingChain(chainName,  line);
    }
    else if ((line.substr(0,6) == "SEQRES")) {

        char chainName = line.at(11);

        if (chainName == ' ') {

            printf("Only 1 chain in the PDB file. So exiting.....\n");
            exit(0);
        }

        insertPDBTextToCorrespondingChain(chainName,  line);
    }
    else if ((line.substr(0,3) == "TER")) {

        char chainName = line.at(21);
        insertPDBTextToCorrespondingChain(chainName,  line);
    }
}

void PDB:: insertPDBTextToCorrespondingChain(char chainName, string line)
{
    if(chainCharToInt.find(chainName) == chainCharToInt.end())
    {
        chainCharToInt[chainName] = numOfChains;
        chainIntToChar[numOfChains] = chainName;
        chainsPDBText[numOfChains].push_back(line);
        numOfChains++;
    }
    else {

        int chainNumber = chainCharToInt[chainName];
        chainsPDBText[chainNumber].push_back(line);
    }
}

void PDB:: insertAtomToCorrespondingChain(char chainName, Coordinate coordinate)
{
    int chainNumber = chainCharToInt[chainName];
    chainsAtom[chainNumber].push_back(coordinate);
}

void PDB:: printPDBAccordingToChain()
{
    printf("\nnub of chain: %d\n", numOfChains);
    for (int i = 0; i < numOfChains; i++) {
        size_t textSize = chainsPDBText[i].size();
        for (int j = 0; j < textSize; j++) {
            cout<< i <<" " << j << " " <<chainsPDBText[i][j] <<"\n";
        }
    }
}

void PDB:: printAtomAccordingToChain()
{
    printf("\nnub of chain: %d\n", numOfChains);
    for (int i = 0; i < numOfChains; i++) {
        size_t atomNumber = chainsAtom[i].size();
        for (int j = 0; j < atomNumber; j++) {
            cout<< i <<" " << j << "  x: " <<chainsAtom[i][j].xCor<<"  y:  "<<chainsAtom[i][j].yCor <<"  z:  "<< chainsAtom[i][j].zCor <<"\n";
        }
    }
}

int PDB::getNumOfChains()
{
    return numOfChains;
}

char PDB::getChainNameCorrespondingToNum(int chainNumber)
{
    return chainIntToChar[chainNumber];
}

vector <Coordinate> PDB::getAtomsCoordinateOfAChain(int chainNumber)
{
    return chainsAtom[chainNumber];
}

// Find the smallest bounding box points that enclose the pdb structure
BoundingCoordinate PDB::findEachAxisBoundaryPoints(void)
{
    double const infinity = 10000;

    double xMin = infinity;
    double xMax = -infinity;
    double yMin = infinity;
    double yMax = -infinity;
    double zMin = infinity;
    double zMax = -infinity;


    for (int i = 0; i < numOfAtoms; i++) {

        ATOM atom = atoms[i];

        if (atom.getCoordinate().xCor < xMin) {
            xMin = atom.getCoordinate().xCor;
        }

        if (atom.getCoordinate().xCor > xMax) {
            xMax = atom.getCoordinate().xCor;
        }

        if (atom.getCoordinate().yCor < yMin) {
            yMin = atom.getCoordinate().yCor;
        }

        if (atom.getCoordinate().yCor > yMax) {
            yMax = atom.getCoordinate().yCor;
        }
        if (atom.getCoordinate().zCor < zMin) {
            zMin = atom.getCoordinate().zCor;
        }

        if (atom.getCoordinate().zCor > zMax) {
            zMax = atom.getCoordinate().zCor;
        }
    }


    BoundingCoordinate boundingCoordinates;

    boundingCoordinates.minPoints.xCor = xMin;
    boundingCoordinates.minPoints.yCor = yMin;
    boundingCoordinates.minPoints.zCor = zMin;


    boundingCoordinates.maxPoints.xCor = xMax;
    boundingCoordinates.maxPoints.yCor = yMax;
    boundingCoordinates.maxPoints.zCor = zMax;


    printf("\nPDB bounding box: xMin: %.3f xMax: %.3f  yMin= %.3f yMax= %.3f zMin= %.3f zMax= %.3f\n", xMin, xMax, yMin, yMax, zMin, zMax);

    return boundingCoordinates;
}


BoundingCoordinate PDB::findEachAxisBoundaryPointsbasedOnChain(int chainNum)
{
    double const infinity = 10000;

    double xMin = infinity;
    double xMax = -infinity;
    double yMin = infinity;
    double yMax = -infinity;
    double zMin = infinity;
    double zMax = -infinity;

    size_t numOfAtomInChain = chainsAtom[chainNum].size();

    for (int i = 0; i < numOfAtomInChain; i++) {

        Coordinate coordinate = chainsAtom[chainNum][i];

        if (coordinate.xCor < xMin) {
            xMin = coordinate.xCor;
        }

        if (coordinate.xCor > xMax) {
            xMax = coordinate.xCor;
        }

        if (coordinate.yCor < yMin) {
            yMin = coordinate.yCor;
        }

        if (coordinate.yCor > yMax) {
            yMax = coordinate.yCor;
        }
        if (coordinate.zCor < zMin) {
            zMin = coordinate.zCor;
        }

        if (coordinate.zCor > zMax) {
            zMax = coordinate.zCor;
        }
    }

    BoundingCoordinate boundingCoordinates;

    boundingCoordinates.minPoints.xCor = xMin;
    boundingCoordinates.minPoints.yCor = yMin;
    boundingCoordinates.minPoints.zCor = zMin;

    boundingCoordinates.maxPoints.xCor = xMax;
    boundingCoordinates.maxPoints.yCor = yMax;
    boundingCoordinates.maxPoints.zCor = zMax;

    printf("\nChain:%d have %d atoms\n",chainNum, (int)numOfAtomInChain);
    printf("Coordinates: xMin: %.3f xMax: %.3f  yMin= %.3f yMax= %.3f zMin= %.3f zMax= %.3f\n", xMin, xMax, yMin, yMax, zMin, zMax);

    return boundingCoordinates;
}

PDB:: ~PDB() {

}
