import java.util.ArrayList;
import java.util.HashMap;
import java.lang.Math;

public class Functions {
    
    int totalNumberOfMolecules;
    int totalSolventAtoms;
    int totalSolventMolecules;
    int numberOfTypesOfSolutes;
    HashMap<String, Integer> typesOfSolutesWithNumbersOfEachType;
    int totalSoluteMolecules;
    ArrayList<Integer> indexOfOxygenAtoms;
    double[] f4Data;
    double[] binningOfF4DataAlongX;
    double[] binningOfF4DataAlongY;
    double[] binningOfF4DataAlongZ;
    ArrayList<ArrayList<Double>> posDataWithChangeOfCoordinates;

    public ArrayList<Integer> calcNumberOfMoleculesData(ArrayList<String> nameCodeData) {
        ArrayList<Integer> numberOfMoleculesData = new ArrayList<>();
        typesOfSolutesWithNumbersOfEachType = new HashMap<>();
        numberOfMoleculesData.add(nameCodeData.size());
        indexOfOxygenAtoms = new ArrayList<>();
        int tempIndexCountForOxygenAtoms = 0;
        int countSolventAtoms = 0;
        int countSoluteMolecules = 0;

        for(int i = 0; i < nameCodeData.size(); i = i + 1) {
            String nameCode = nameCodeData.get(i);

            if(nameCode.equals("SOL") || nameCode.equals("wat")) {
                countSolventAtoms = countSolventAtoms + 1;
                if(tempIndexCountForOxygenAtoms == 0) {
                    indexOfOxygenAtoms.add(i);
                }
                tempIndexCountForOxygenAtoms = tempIndexCountForOxygenAtoms + 1;
                if(tempIndexCountForOxygenAtoms == 4) {
                    tempIndexCountForOxygenAtoms = 0;
                }
            }
            else {
                if(!typesOfSolutesWithNumbersOfEachType.containsKey(nameCode)) {
                    typesOfSolutesWithNumbersOfEachType.put(nameCode, 1);
                    numberOfTypesOfSolutes = numberOfTypesOfSolutes + 1;
                }
                else {
                    int nameCodeValue = typesOfSolutesWithNumbersOfEachType.get(nameCode);
                    typesOfSolutesWithNumbersOfEachType.put(nameCode, nameCodeValue + 1);
                }

                countSoluteMolecules = countSoluteMolecules + 1;
            }
        }

        totalSolventAtoms = countSolventAtoms;
        totalSolventMolecules = countSolventAtoms / 4;
        totalSoluteMolecules = countSoluteMolecules;
        totalNumberOfMolecules = totalSoluteMolecules + totalSolventMolecules;

        numberOfMoleculesData.add(totalNumberOfMolecules);
        numberOfMoleculesData.add(totalSolventAtoms);
        numberOfMoleculesData.add(totalSolventMolecules);
        numberOfMoleculesData.add(totalSoluteMolecules);
        numberOfMoleculesData.add(numberOfTypesOfSolutes);

        return numberOfMoleculesData;
    }

    public NeighbourList[] neighbourListGenerator(ArrayList<ArrayList<Double>> posData, double boxX, double boxY, double boxZ) {
        NeighbourList[] neighbourList = new NeighbourList[indexOfOxygenAtoms.size()];
        int indexForNeighbourList = 0;
        int indexOfOxygenAtomsInNeighbourList = 0;
        double lowerLimit = 0.18;
        double upperLimit = 0.35;
        double dx = 0;
        double dy = 0;
        double dz = 0;
        double dist = 0;

        for(int i = 0; i < indexOfOxygenAtoms.size(); i = i + 1) {
            neighbourList[i] = new NeighbourList();
        }

        for(int i : this.indexOfOxygenAtoms) {
            neighbourList[indexForNeighbourList].parentOxygenIndex = i; 
            for(int j : this.indexOfOxygenAtoms) {
                if(i!=j) {
                    dx = Math.abs(posData.get(j).get(0) - posData.get(i).get(0));
                    dy = Math.abs(posData.get(j).get(1) - posData.get(i).get(1));
                    dz = Math.abs(posData.get(j).get(2) - posData.get(i).get(2));

                    if(dx >= boxX * 0.5) {
                        dx = boxX - dx;
                    }
                    if(dy >= boxY * 0.5) {
                        dy = boxY - dx;
                    }
                    if(dz >= boxZ * 0.5) {
                        dz = boxZ - dz;
                    }

                    dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

                    if(dist < upperLimit && dist > lowerLimit) {
                        neighbourList[indexForNeighbourList].neighboursToParent.add(j);
                    }
                }
            }
            indexForNeighbourList = indexForNeighbourList + 1;
        }

        return neighbourList;
    }

    public double calcF4Data(NeighbourList[] neighbourList, ArrayList<ArrayList<Double>> posData, double boxX, double boxY, double boxZ) {
        f4Data = new double[neighbourList.length];
        double f4 = 0;
        double f4SummationForAverage = 0;
        double averageF4 = 0;
        double numberOfPairsParticipatingInCalculation = 0;
        double dx = 0;
        double dy = 0;
        double dz = 0;
        double dist1 = 0;
        double dist2 = 0;
        double dist3 = 0;
        double dist4 = 0;
        int dihedralAtom1Index = 0;
        int dihedralAtom2Index = 0;
        int dihedralAtom3Index = 0;
        int dihedralAtom4Index = 0;
        double[] vector12 = new double[3];
        double[] vector23 = new double[3];
        double[] vector34 = new double[3];
        double[] vector123Cross = new double[3];
        double[] vector234Cross = new double[3];
        double anglePhi = 0;
        double dotProduct = 0;
        double magnitudeVector123Cross = 0;
        double magnitudeVector234Cross = 0;
        int centralOxygenIndex = 0;
        int currentNeighbourIndex = 0;

        for(int i = 0; i < neighbourList.length; i = i + 1) {
            centralOxygenIndex = neighbourList[i].parentOxygenIndex;
            for(int j : neighbourList[i].neighboursToParent) {
                currentNeighbourIndex = j;

                dx = Math.abs(posData.get(centralOxygenIndex).get(0) - posData.get(currentNeighbourIndex + 1).get(0));
                dy = Math.abs(posData.get(centralOxygenIndex).get(1) - posData.get(currentNeighbourIndex + 1).get(1));
                dz = Math.abs(posData.get(centralOxygenIndex).get(2) - posData.get(currentNeighbourIndex + 1).get(2));
                if(dx >= boxX * 0.5) {
                    dx = boxX - dx;
                }
                if(dy >= boxY * 0.5) {
                    dy = boxY - dx;
                }
                if(dz >= boxZ * 0.5) {
                    dz = boxZ - dz;
                }
                dist1 = Math.sqrt(dx * dx + dy * dy + dz * dz);

                dx = Math.abs(posData.get(centralOxygenIndex).get(0) - posData.get(currentNeighbourIndex + 2).get(0));
                dy = Math.abs(posData.get(centralOxygenIndex).get(1) - posData.get(currentNeighbourIndex + 2).get(1));
                dz = Math.abs(posData.get(centralOxygenIndex).get(2) - posData.get(currentNeighbourIndex + 2).get(2));
                if(dx >= boxX * 0.5) {
                    dx = boxX - dx;
                }
                if(dy >= boxY * 0.5) {
                    dy = boxY - dx;
                }
                if(dz >= boxZ * 0.5) {
                    dz = boxZ - dz;
                }
                dist2 = Math.sqrt(dx * dx + dy * dy + dz * dz);

                dx = Math.abs(posData.get(centralOxygenIndex + 1).get(0) - posData.get(currentNeighbourIndex).get(0));
                dy = Math.abs(posData.get(centralOxygenIndex + 1).get(1) - posData.get(currentNeighbourIndex).get(1));
                dz = Math.abs(posData.get(centralOxygenIndex + 1).get(2) - posData.get(currentNeighbourIndex).get(2));
                if(dx >= boxX * 0.5) {
                    dx = boxX - dx;
                }
                if(dy >= boxY * 0.5) {
                    dy = boxY - dx;
                }
                if(dz >= boxZ * 0.5) {
                    dz = boxZ - dz;
                }
                dist3 = Math.sqrt(dx * dx + dy * dy + dz * dz);

                dx = Math.abs(posData.get(centralOxygenIndex + 2).get(0) - posData.get(currentNeighbourIndex).get(0));
                dy = Math.abs(posData.get(centralOxygenIndex + 2).get(1) - posData.get(currentNeighbourIndex).get(1));
                dz = Math.abs(posData.get(centralOxygenIndex + 2).get(2) - posData.get(currentNeighbourIndex).get(2));
                if(dx >= boxX * 0.5) {
                    dx = boxX - dx;
                }
                if(dy >= boxY * 0.5) {
                    dy = boxY - dx;
                }
                if(dz >= boxZ * 0.5) {
                    dz = boxZ - dz;
                }
                dist4 = Math.sqrt(dx * dx + dy * dy + dz * dz);

                if(Math.min(dist1, Math.min(dist2, Math.min(dist3, dist4))) == dist4) {
                    dihedralAtom1Index = centralOxygenIndex + 1;
                    dihedralAtom2Index = centralOxygenIndex;
                    dihedralAtom3Index = currentNeighbourIndex;
                    if(Math.max(dist1, Math.max(dist2, Math.max(dist3, dist4))) == dist1) {
                        dihedralAtom4Index = currentNeighbourIndex + 1;
                    }
                    else {
                        dihedralAtom4Index = currentNeighbourIndex + 2;
                    }
                }

                if(Math.min(dist1, Math.min(dist2, Math.min(dist3, dist4))) == dist3) {
                    dihedralAtom1Index = centralOxygenIndex + 2;
                    dihedralAtom2Index = centralOxygenIndex;
                    dihedralAtom3Index = currentNeighbourIndex;
                    if(Math.max(dist1, Math.max(dist2, Math.max(dist3, dist4))) == dist1) {
                        dihedralAtom4Index = currentNeighbourIndex + 1;
                    }
                    else {
                        dihedralAtom4Index = currentNeighbourIndex + 2;
                    }
                }

                if(Math.min(dist1, Math.min(dist2, Math.min(dist3, dist4))) == dist1) {
                    dihedralAtom1Index = currentNeighbourIndex + 2;
                    dihedralAtom2Index = currentNeighbourIndex;
                    dihedralAtom3Index = centralOxygenIndex;
                    if(Math.max(dist1, Math.max(dist2, Math.max(dist3, dist4))) == dist3) {
                        dihedralAtom4Index = centralOxygenIndex + 1;
                    }
                    else {
                        dihedralAtom4Index = centralOxygenIndex + 2;
                    }
                }

                if(Math.min(dist1, Math.min(dist2, Math.min(dist3, dist4))) == dist2) {
                    dihedralAtom1Index = currentNeighbourIndex + 1;
                    dihedralAtom2Index = currentNeighbourIndex;
                    dihedralAtom3Index = centralOxygenIndex;
                    if(Math.max(dist1, Math.max(dist2, Math.max(dist3, dist4))) == dist3) {
                        dihedralAtom4Index = centralOxygenIndex + 1;
                    }
                    else {
                        dihedralAtom4Index = centralOxygenIndex + 2;
                    }
                }

                vector12[0] = posData.get(dihedralAtom2Index).get(0) - posData.get(dihedralAtom1Index).get(0);
                vector12[1] = posData.get(dihedralAtom2Index).get(1) - posData.get(dihedralAtom1Index).get(1);
                vector12[2] = posData.get(dihedralAtom2Index).get(2) - posData.get(dihedralAtom1Index).get(2);

                vector23[0] = posData.get(dihedralAtom3Index).get(0) - posData.get(dihedralAtom2Index).get(0);
                vector23[1] = posData.get(dihedralAtom3Index).get(1) - posData.get(dihedralAtom2Index).get(1);
                vector23[2] = posData.get(dihedralAtom3Index).get(2) - posData.get(dihedralAtom2Index).get(2);

                vector34[0] = posData.get(dihedralAtom4Index).get(0) - posData.get(dihedralAtom3Index).get(0);
                vector34[1] = posData.get(dihedralAtom4Index).get(1) - posData.get(dihedralAtom3Index).get(1);
                vector34[2] = posData.get(dihedralAtom4Index).get(2) - posData.get(dihedralAtom3Index).get(2);

                vector123Cross[0] = vector12[1] * vector23[2] - vector12[2] * vector23[1];
                vector123Cross[1] = vector12[2] * vector23[0] - vector12[0] * vector23[2];
                vector123Cross[2] = vector12[0] * vector23[1] - vector12[1] * vector23[0];

                vector234Cross[0] = vector23[1] * vector34[2] - vector23[2] * vector34[1];
                vector234Cross[1] = vector23[2] * vector34[0] - vector23[0] * vector34[2];
                vector234Cross[2] = vector23[0] * vector34[1] - vector23[1] * vector34[0];

                for(int k = 0; k < vector123Cross.length; k = k + 1) {
                    dotProduct = dotProduct + vector123Cross[k] * vector234Cross[k];
                }

                magnitudeVector123Cross = Math.sqrt(vector123Cross[0] * vector123Cross[0] + vector123Cross[1] * vector123Cross[1] + vector123Cross[2] * vector123Cross[2]);
                magnitudeVector234Cross = Math.sqrt(vector234Cross[0] * vector234Cross[0] + vector234Cross[1] * vector234Cross[1] + vector234Cross[2] * vector234Cross[2]);
                anglePhi = Math.acos(dotProduct / (magnitudeVector123Cross * magnitudeVector234Cross));
                f4 = f4 + Math.cos(3 * anglePhi);
                f4SummationForAverage = f4SummationForAverage + Math.cos(3 *  anglePhi);
                numberOfPairsParticipatingInCalculation = numberOfPairsParticipatingInCalculation + 1;

                dotProduct = 0;
            }

            if(neighbourList[i].neighboursToParent.size() != 0) {
                f4 = f4 / (neighbourList[i].neighboursToParent.size());
                f4Data[i] = f4;
                f4 = 0;
            }
            else {
                f4Data[i] = 10200;
            }
        }

        averageF4 = f4SummationForAverage / numberOfPairsParticipatingInCalculation;
        return averageF4;
    }
    
    public void binF4Data(double[] f4Data, NeighbourList[] neighbourList, ArrayList<ArrayList<Double>> posData, double boxX, double boxY, double boxZ) {
        int binIdentifierX;
        int binIdentifierY;
        int binIdentifierZ;
        double binSep = 0.1;
        int numberOfBinsAlongX = (int)(boxX / binSep) + 1;
        int numberOfBinsAlongY = (int)(boxY / binSep) + 1;
        int numberOfBinsAlongZ = (int)(boxZ / binSep) + 1;
        binningOfF4DataAlongX = new double[numberOfBinsAlongX];
        binningOfF4DataAlongY = new double[numberOfBinsAlongY];
        binningOfF4DataAlongZ = new double[numberOfBinsAlongZ];
        int[] numberOfMoleculesBinwiseAlongX = new int[numberOfBinsAlongX];
        int[] numberOfMoleculesBinwiseAlongY = new int[numberOfBinsAlongY];
        int[] numberOfMoleculesBinwiseAlongZ = new int[numberOfBinsAlongZ];

        posDataWithChangeOfCoordinates = new ArrayList<>();
        ArrayList<Double> tempListForPosData = new ArrayList<>();
        boolean flag = true;

        for(int i = 0; i < posData.size(); i = i + 1) {
            if(posData.get(i).get(0) < 0 || posData.get(i).get(1) < 0 || posData.get(i).get(2) < 0) {
                flag = false;
            }
        }
        
        if(!flag) {
            for(int i = 0; i < posData.size(); i = i + 1) {
                tempListForPosData.add(posData.get(i).get(0) + (boxX / 2));
                tempListForPosData.add(posData.get(i).get(1) + (boxY / 2));
                tempListForPosData.add(posData.get(i).get(2) + (boxZ / 2));
                posDataWithChangeOfCoordinates.add(tempListForPosData);
                tempListForPosData = new ArrayList<>();
            }
        }
        else {
            for(int i = 0; i < posData.size(); i = i + 1) {
                tempListForPosData.add(posData.get(i).get(0));
                tempListForPosData.add(posData.get(i).get(1));
                tempListForPosData.add(posData.get(i).get(2));
                posDataWithChangeOfCoordinates.add(tempListForPosData);
                tempListForPosData = new ArrayList<>();
            }
        }

        for(int i = 0; i < f4Data.length; i = i + 1) {
            binIdentifierX = (int)(posDataWithChangeOfCoordinates.get(neighbourList[i].parentOxygenIndex).get(0) / binSep) + 1;
            binIdentifierY = (int)(posDataWithChangeOfCoordinates.get(neighbourList[i].parentOxygenIndex).get(1) / binSep) + 1;
            binIdentifierZ = (int)(posDataWithChangeOfCoordinates.get(neighbourList[i].parentOxygenIndex).get(2) / binSep) + 1;

            if(f4Data[i] != 10200) {
                numberOfMoleculesBinwiseAlongX[binIdentifierX - 1] = numberOfMoleculesBinwiseAlongX[binIdentifierX - 1] + 1;
                numberOfMoleculesBinwiseAlongY[binIdentifierY - 1] = numberOfMoleculesBinwiseAlongY[binIdentifierY - 1] + 1;
                numberOfMoleculesBinwiseAlongZ[binIdentifierZ - 1] = numberOfMoleculesBinwiseAlongZ[binIdentifierZ - 1] + 1;

                binningOfF4DataAlongX[binIdentifierX - 1] = binningOfF4DataAlongX[binIdentifierX - 1] + f4Data[i];
                binningOfF4DataAlongY[binIdentifierY - 1] = binningOfF4DataAlongY[binIdentifierY - 1] + f4Data[i];
                binningOfF4DataAlongZ[binIdentifierZ - 1] = binningOfF4DataAlongZ[binIdentifierZ - 1] + f4Data[i];     
            }        
        }
        
        for(int i = 0; i < numberOfBinsAlongX; i = i + 1) {
            if(numberOfMoleculesBinwiseAlongX[i] != 0) {
                binningOfF4DataAlongX[i] = binningOfF4DataAlongX[i] / numberOfMoleculesBinwiseAlongX[i];
            }
            else {
                binningOfF4DataAlongX[i] = 10200;
            }
        }
        for(int i = 0; i < numberOfBinsAlongY; i = i + 1) {
            if(numberOfMoleculesBinwiseAlongY[i] != 0) {
                binningOfF4DataAlongY[i] = binningOfF4DataAlongY[i] / numberOfMoleculesBinwiseAlongY[i];
            }
            else {
                binningOfF4DataAlongY[i] = 10200;
            }
        }
        for(int i = 0; i < numberOfBinsAlongY; i = i + 1) {
            if(numberOfMoleculesBinwiseAlongZ[i] != 0) {
                binningOfF4DataAlongZ[i] = binningOfF4DataAlongZ[i] / numberOfMoleculesBinwiseAlongZ[i];
            }
            else {
                binningOfF4DataAlongZ[i] = 10200;
            }
        }

        return;
    }

}
