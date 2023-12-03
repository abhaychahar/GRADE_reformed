import java.util.Scanner;
import java.util.ArrayList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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

        try {
            File outputFile = new File("binX.xvg");
            if(outputFile.createNewFile()) {
                System.out.println("File created: " + outputFile.getName());
            } 
            else {
                System.out.println("File already exists.");
            }
        } 
        catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

        try {
            FileWriter outputFileWriter = new FileWriter("binX.xvg");
            outputFileWriter.close();
            System.out.println("Existing file contents cleared.");
        } 
        catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        
        ArrayList<ArrayList<Double>> posData = new ArrayList<>();
        ArrayList<String> nameCodeData = new ArrayList<>();
        ArrayList<Integer> numberOfMoleculesData;
        double[] f4Data;
        double averageF4;
        NeighbourList[] neighbourList;
        BinningDataList binningDataList;
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
                    
                    if(lineCount < rows + 3) {
                        tempPosData.add(Double.parseDouble(lineAsString.substring(23, 28)));
                        tempPosData.add(Double.parseDouble(lineAsString.substring(31, 36)));
                        tempPosData.add(Double.parseDouble(lineAsString.substring(39, 44)));

                        nameCodeData.add(lineAsString.substring(5, 8));
                        posData.add(tempPosData);
                    }
                    else if(lineCount == rows + 3) {
                        boxX = Double.parseDouble(lineAsString.substring(3, 10));
                        boxY = Double.parseDouble(lineAsString.substring(13, 20));
                        boxZ = Double.parseDouble(lineAsString.substring(23, 30));
                        frameCount = frameCount + 1;
                        lineCount = 0;
                    }                                        
                }

                if(frameCount == nextFrame) {
                    nextFrame = nextFrame + 1;

                    numberOfMoleculesData = func.calcNumberOfMoleculesData(nameCodeData);
                    neighbourList = func.neighbourListGenerator(posData, boxX, boxY, boxZ);
                    averageF4 = func.calcF4Data(neighbourList, posData, boxX, boxY, boxZ);
                    f4Data = func.f4Data;
                    func.binF4Data(f4Data, neighbourList, posData, boxX, boxY, boxZ);

                    double tempSum = 0;
                    int tempNum = 0;

                    try(FileWriter fw = new FileWriter("binX.xvg", true); BufferedWriter bw = new BufferedWriter(fw); PrintWriter outputPrintWriter = new PrintWriter(bw))
                    {
                        for(int i = 0; i < func.binningOfF4DataAlongX.length; i = i + 1) {
                            if(func.binningOfF4DataAlongX[i] != 10200) {
                                outputPrintWriter.println(" " + (i+1) + "      " + func.binningOfF4DataAlongX[i]);
                                tempSum = tempSum + func.binningOfF4DataAlongX[i];
                                tempNum = tempNum + 1;
                            }
                        }
                        outputPrintWriter.println();

                        outputPrintWriter.close();
                        bw.close();
                        fw.close();
                    } 
                    catch (IOException e) {
                        System.out.println("An error occurred.");
                        e.printStackTrace();
                    }

                    System.out.println("Average f4 value calculated separately : " + averageF4);
                    System.out.println("Average f4 value calculated from binning : " + (tempSum / tempNum));
                }
            }
            reader.close();
        }
        catch(FileNotFoundException e) {
            System.out.println("Error occurred (File not found)");
        }
    }

}