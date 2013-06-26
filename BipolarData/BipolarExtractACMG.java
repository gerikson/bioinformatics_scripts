/*
 * Script that extracts the number of ACMG variants per category per genotype
 */


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author gerikson
 */

public class BipolarExtractACMG {
    
    public static String[] header;
    public static int ACMGcolumn;
    public static int notesColumn;
   
    
    /**
     *
     */
    public static Map<String, Map<String, String>> GenotypeMap = new HashMap<String, Map<String, String>>();

        
    public static int getHeaderIndex(String var) {
        int index = 0;
        for (int i = 0; i<header.length; i++) {
            if (header[i].contains(var)) {
                index = i;
                
            }
        }
    return index;
    }
     
    public static void readFile(String inputFile) throws FileNotFoundException, IOException{
            BufferedReader reader = new BufferedReader( new FileReader (inputFile));
            String         line = null;
            int lineNumber = 0;
            
            while( ( line = reader.readLine() ) != null ) {
                lineNumber ++;
                if (lineNumber == 1) {
                	System.out.println("Header is: " + line);
                	header = line.split("\t");
                	ACMGcolumn = getHeaderIndex("ACMG_Score_Clinical");
                	notesColumn = getHeaderIndex("Notes");
                        continue;
                }
                
                if (lineNumber % 10000 == 0) {
                    System.out.println("Datacount is" + lineNumber);
                }
                String[] tempLine = line.split("\t");
                countACMGperGenotype(tempLine);
            }
            
            
         WriteToFile();
    }
     
    public static void countACMGperGenotype(String[] line) {
                    if (line[ACMGcolumn].contains("1~") || line[ACMGcolumn].contains("2~") || line[ACMGcolumn].contains("2*~") ||
                    line[ACMGcolumn].contains("3~") || line[ACMGcolumn].contains("4~") || line[ACMGcolumn].contains("5~") ||
                    line[ACMGcolumn].contains("6~")){
                        //Extract the notes column and see how many genotyps are in there, put the genotypes in the hashMap
                        String notes = line[notesColumn];
                        String[] gen = notes.split("-");
                     
                        
                        for (int i =0; i < gen.length; i++) {
                            //Extract only the genotype name
                            String[] tem = gen[i].split(":");
                            String geno = tem[0];
                       //     System.out.println(geno);
                            Map<String, String> temp;
                            //check if this genotype was already present
                            if (GenotypeMap.containsKey(geno)) {
                            temp = GenotypeMap.get(geno);
                                
                                //figure out what category it is and add that to the list
                                if (line[ACMGcolumn].contains("1~")) {
                                    int t = Integer.parseInt(temp.get("1")) + 1;
                                    temp.put("1", Integer.toString(t));
                                }                                
                                else if (line[ACMGcolumn].contains("2~")) {
                                    int t = Integer.parseInt(temp.get("2")) + 1;
                                    temp.put("2", Integer.toString(t));
                                }
                                else if (line[ACMGcolumn].contains("2*~")) {
                                    int t = Integer.parseInt(temp.get("2*")) + 1;
                                    temp.put("2*", Integer.toString(t));
                                }
                                else if (line[ACMGcolumn].contains("3~")) {
                                    int t = Integer.parseInt(temp.get("3")) + 1;
                                    temp.put("3", Integer.toString(t));
                                }
                                else if (line[ACMGcolumn].contains("4~")) {
                                    try {
                                    int t = Integer.parseInt(temp.get("4")) + 1;
                                    temp.put("4", Integer.toString(t));
                                    } catch(Exception e) {
                                        System.out.println("Can't parse int" + temp.get("4"));
                                    }
                                }
                                else if (line[ACMGcolumn].contains("5~")) {
                                //    if (temp.size() < 5) {
                                  //      System.out.print("Small ArrayList size. Geno is: " + geno);
                                  //  } else {
                                    int t = Integer.parseInt(temp.get("5")) + 1;
                                    temp.put("5", Integer.toString(t));
                                  //  }
                                }
                                else if (line[ACMGcolumn].contains("6~")) {
                                    int t = Integer.parseInt(temp.get("6")) + 1;
                                    temp.put("6", Integer.toString(t));
                                }
                                
                                //insert the new counters into the map
                                GenotypeMap.put(geno, temp);
                                
                            } 
                            //This genotypes was not previously found just add it to the GenotypeMap
                           
                            else {
                                
                                Map<String, String> tempTwo = new HashMap<String, String>();
                                tempTwo.put("1", "0");
                                tempTwo.put("2", "0");
                                tempTwo.put("2*", "0");   
                                tempTwo.put("3", "0");
                                tempTwo.put("4", "0");
                                tempTwo.put("5", "0");
                                tempTwo.put("6", "0");
                                
                                
                                //figure out what category it is and add that to the list
                                if (line[ACMGcolumn].contains("1~")) {
                                    tempTwo.put("1", "1");
                                }                                
                                else if (line[ACMGcolumn].contains("2~")) {
                                    tempTwo.put("2", "1");
                                }
                                else if (line[ACMGcolumn].contains("2*~")) {
                                    tempTwo.put("2*", "1");
                                }
                                else if (line[ACMGcolumn].contains("3~")) {
                                    tempTwo.put("3", "1");
                                }
                                else if (line[ACMGcolumn].contains("4~")) {
                                    tempTwo.put("4", "1");
                                }
                                else if (line[ACMGcolumn].contains("5~")) {
                                    tempTwo.put("5", "1");
                                }
                                if (line[ACMGcolumn].contains("6~")) {
                                    tempTwo.put("6", "1");
                                }
                                
                                GenotypeMap.put(geno, tempTwo);
                            }
                        }
                        
        }
    }
    
    public static void WriteToFile() throws IOException {
                           BufferedWriter outfile = new BufferedWriter(
                                 new FileWriter("ACMG_counter_bygenotype.txt", true));     
                //Print all key values from Genotype Map     
                           outfile.write("Genotype" +'\t' + "Category 1" + '\t' + "Category 2" + '\t' + "Category 2*" + '\t');
                           outfile.write("Category 3" +  '\t' + "Category 4" + '\t' + "Category 5" + '\t' + "Category 6" + "\n");      
                for (Map.Entry<String, Map<String, String>> entry: GenotypeMap.entrySet()) {
                    String key = entry.getKey();
                    Map<String, String> value = entry.getValue();
                    
                    
                    outfile.write(key +'\t' + value.get("1") + "\t"  + value.get("2") + "\t"  + value.get("2*") + "\t"  + value.get("3") + "\t"
                             + value.get("4") + "\t" + value.get("5") + "\t" + value.get("6"));        		  
                    outfile.write("\n");
                }
      outfile.close();         
    }
    
    public static void ReadExistingCounter(String filename, String InputFile) throws FileNotFoundException, IOException {
            BufferedReader reader = new BufferedReader( new FileReader (filename));
            String         line = null;
            int lineNumber = 0;
            
            while( ( line = reader.readLine() ) != null ) {
                lineNumber ++;
                if (lineNumber == 1) {
                    System.out.println("First line counter");
                    continue;
                } 
                
                String[] spLine = line.split("\t");
                                Map<String, String> tempTwo = new HashMap<String, String>();
                                tempTwo.put("1", spLine[1]);
                                tempTwo.put("2", spLine[2]);
                                tempTwo.put("2*",spLine[3]);   
                                tempTwo.put("3", spLine[4]);
                                tempTwo.put("4", spLine[5]);
                                tempTwo.put("5", spLine[6]);
                                tempTwo.put("6", spLine[7]);
                GenotypeMap.put(spLine[0], tempTwo);
                
            }     
           readFile(InputFile);
    
   }
    
    public static void main(String[] args) throws IOException {
           System.out.println("Usage javac Stats.java full/path/to/file /path/to/existing/counter");
           String InputFile = args[0];
           String FileTwo = args[1];
           System.out.println("Input File is " + InputFile);
           //readFile(InputFile);
           ReadExistingCounter(FileTwo, InputFile);
    }    
}
