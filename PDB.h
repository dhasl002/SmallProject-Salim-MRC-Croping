
#ifndef ____PDB__
#define ____PDB__

#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include "ATOM.h"
#include <map>

using namespace std;


class PDB
{
private:

    int numOfAtoms;
    int numOfChains;
    vector <ATOM> atoms;
    map <char, int> chainCharToInt;
    map <int, char> chainIntToChar;
    std:: vector <std:: vector< string >> chainsPDBText;
    std:: vector <std:: vector< Coordinate >> chainsAtom;

public:
    PDB(){
        numOfAtoms = 0;
        numOfChains = 0;

        for (int i = 0; i < 30; i++) {
            vector<string> chain;
            chainsPDBText.push_back(chain);
        }

        for (int i = 0; i < 30; i++) {
            vector<Coordinate> chain;
            chainsAtom.push_back(chain);
        }
    };

    void readPDBFile(string fileName);
    void writeToFileInPDBFormat(string fileName);
    BoundingCoordinate findEachAxisBoundaryPoints(void);
    BoundingCoordinate findEachAxisBoundaryPointsbasedOnChain(int chainNum);
    void getAtoms(string line);
    void insertPDBTextToCorrespondingChain(char chainName, string line);
    void insertAtomToCorrespondingChain(char chainName, Coordinate coordinate);
    void printPDBAccordingToChain();
    void writePDBAccordingToChain(string filePath);
    void printAtomAccordingToChain();
    int  getNumOfChains();
    char getChainNameCorrespondingToNum(int chainNumber);
    vector <Coordinate> getAtomsCoordinateOfAChain(int chainNumber);
    ~PDB();
};

#endif /* defined(____PDB__) */
