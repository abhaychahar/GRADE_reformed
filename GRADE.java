import java.util.Scanner;
import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;

public class GRADE {
    
    public static void main(String[] args) {

        Functions func = new Functions();

        Scanner scanFileName = new Scanner(System.in);
        String fileName = scanFileName.nextLine();
        if(!fileName.endsWith(".gro")) {
            fileName = fileName + ".gro";
        }
        scanFileName.close();
        
        ArrayList<ArrayList<Double>> posData = new ArrayList<>();
        ArrayList<String> nameCodeData = new ArrayList<>();
        ArrayList<Integer> numberOfMoleculesData = new ArrayList<>();
        double boxX = 0;
        double boxY = 0;
        double boxZ = 0;
        int lineCount = 0;
        int frameCount = 0;
        int nextFrame = 1;
        int rows = 0;

        try {
            File file = new File(fileName);
            Scanner reader = new Scanner(file);
            while(reader.hasNextLine()) {
                lineCount = lineCount + 1;

                if(lineCount == 1) {
                    reader.nextLine();
                }
                else if(lineCount == 2) {
                    rows = reader.nextInt();
                    reader.nextLine();
                }
                else {
                    String lineAsString = reader.nextLine();
                    ArrayList<Double> tempPosData = new ArrayList<>();
                    
                    if(lineAsString.length() > 30) {
                        tempPosData.add(Double.parseDouble(lineAsString.substring(23, 28)));
                        tempPosData.add(Double.parseDouble(lineAsString.substring(31, 36)));
                        tempPosData.add(Double.parseDouble(lineAsString.substring(39, 44)));

                        nameCodeData.add(lineAsString.substring(5, 8));
                        posData.add(tempPosData);
                    }
                    else {
                        boxX = Double.parseDouble(lineAsString.substring(3, 10));
                        boxY = Double.parseDouble(lineAsString.substring(13, 20));
                        boxZ = Double.parseDouble(lineAsString.substring(23, 30));
                        frameCount = frameCount + 1;
                    }                                        
                }

                if(frameCount == nextFrame) {
                    nextFrame = nextFrame + 1;

                    numberOfMoleculesData = func.calcNumberOfMoleculesData(nameCodeData);
                }
            }
            reader.close();
        }
        catch(FileNotFoundException e) {
            System.out.println("Error occurred (File not found)");
        }

        System.out.println(numberOfMoleculesData);

    }

}