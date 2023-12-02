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

                    dist = Math.sqrt(dx*dx + dy*dy + dz*dz);

                    if(dist < upperLimit && dist > lowerLimit) {
                        neighbourList[indexForNeighbourList].neighboursToParent.add(j);
                    }
                }
            }
            indexForNeighbourList = indexForNeighbourList + 1;
        }

        return neighbourList;
    }

    // public double calcF4Data(NeighbourList[] neighbourList, ArrayList<ArrayList<Double>> posData, double boxX, double boxY, double boxZ) {
    //     double averageF4 = 0;
    //     double dx = 0;
    //     double dy = 0;
    //     double dz = 0;
    //     double dist = 0;

    //     for(int i = 0; i < neighbourList)

    //     return averageF4;
    // } 

}
