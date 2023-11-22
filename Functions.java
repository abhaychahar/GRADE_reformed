import java.util.ArrayList;
import java.util.HashMap;

public class Functions {
    
    int totalNumberOfMolecules;
    int totalSolventAtoms;
    int totalSolventMolecules;
    int numberOfTypesOfSolutes;
    HashMap<String, Integer> typesOfSolutesWithNumbersOfEachType;
    int totalSoluteMolecules;

    public ArrayList<Integer> calcNumberOfMoleculesData(ArrayList<String> nameCodeData) {
        ArrayList<Integer> numberOfMoleculesData = new ArrayList<>();
        typesOfSolutesWithNumbersOfEachType = new HashMap<>();
        numberOfMoleculesData.add(nameCodeData.size());
        int countSolventAtoms = 0;
        int countSoluteMolecules = 0;

        for(int i = 0; i < nameCodeData.size(); i = i + 1) {
            String nameCode = nameCodeData.get(i);

            if(nameCode.equals("SOL") || nameCode.equals("wat")) {
                countSolventAtoms = countSolventAtoms + 1;
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

}
