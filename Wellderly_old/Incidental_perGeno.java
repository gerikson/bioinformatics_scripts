/*
 * Per each variant calculate how many white, non related have the variant extract 
 * the number of genotypes and the name of the genotypes at the end of the file
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author gerikson
 */
public class Incidental_perGeno {
    
  //  public int[] Acmg1_count;
    public String[] header;
    public String[] headerReduced;
 //   public String genotypes;
    public String inputFile;
 //   public int ACMG_column;
    public String genotNames;
    
    public Incidental_perGeno(String path, String genoPath, String genoFull) throws FileNotFoundException, IOException {
            
            inputFile = path;
            BufferedReader reader = new BufferedReader( new FileReader (genoFull));
            String line = null;
            
            while( ( line = reader.readLine() ) != null ) {

                header = line.split("\t");
               
                System.out.println("header is: " + line + "\n");
   
            }
            
            BufferedReader reader2 = new BufferedReader( new FileReader (genoPath));
            String line2 = null;
            
            while( ( line2 = reader2.readLine() ) != null ) {

                headerReduced = line2.split("\t");
                System.out.println("Header reduced: " + line2 + "\n");
               
            }
            
            System.out.println("Full header size is: " + header.length + "\n");
            System.out.println("Reduced file size is: " + headerReduced.length + "\n");
            
    }
    
     public int getHeaderIndex(String var) {
        int index = 0;
        for (int i = 0; i<header.length; i++) {
            if (header[i].contains(var)) {
                index = i;
                
            }
        }
        
        return index;
    }
    
    public void readFile() throws FileNotFoundException, IOException{
            BufferedReader reader = new BufferedReader( new FileReader (inputFile));
            String         line = null;
            int lineNumber = 0;
            BufferedWriter outfile = new BufferedWriter(new FileWriter("ACMG_IncidentalNoRel.txt", true));    
            while( ( line = reader.readLine() ) != null ) {
                String [] l = line.split("\t");
                lineNumber++;
                if (lineNumber == 1) {
                    continue;
                }
            //    System.out.println(line);
                int t = countACMGperGenotype(l);
              //  WriteToFile(l, t);
                outfile.write(line + "\t" + t + "\t" + genotNames + "\n");
            }
            
            outfile.close();
            
    }
    
    /*
     * count how many genotypes have this mutation 
     */
    
        public int countACMGperGenotype(String[] spLine) throws IOException {
                        genotNames = null;
    			//extract the genotype
    			//String value = spLine[getHeaderIndex("Notes")];
                        String value = spLine[7];
                        System.out.println("Value is :" + value);
    			String[] notes = value.split("=");
                        String temp = notes[0];
    			String allele = temp.replace("Genotype", "");
                        System.out.println("Genotypes is: " + notes[1]);
    			String[] spValue = notes[1].split(":");
                        int count = 0;
    			for (int i =0; i<spValue.length; i++) {
    			//if this is the only allele 
    			if (allele.equals("")) {
			   if (spValue[i].contains("1")) {
                               
                               /*
                                * check if this is one is not africanAmerican or related if not, add it to the list and count
                                * the variants that have zero at the end would be deleted
                                * 
                                */
                               
                                String g = header[i];
                                for (int s = 0; s < headerReduced.length; s++) {
                                    if (g.equals(headerReduced[s])) {
                                        System.out.println("Found");
                                        genotNames = genotNames + headerReduced[s] + "_";
                                        count++;
                                      //  continue;
                                    }
                                }
  			
    			}
                        }else {
    				if (spValue[i].contains(allele)) {
                                /*
                                * check if this is one is not africanAmerican or related if not, add it to the list and count
                                * the variants that have zero at the end would be deleted
                                * 
                                */
                               
                                String g = header[i];
                                for (int s = 0; s <headerReduced.length; s++) {
                                    if (g.equals(headerReduced[s])) {
                                        System.out.println("Found");
                                        genotNames = genotNames + headerReduced[s] + "_";
                                        count ++;
                                     //   continue;
                                    }
                                }   
                            }	    			

                    }
                   
                        
        
    }
                         return count;
 }
        
    public void WriteToFile (String[] line, int count) throws IOException {
                   BufferedWriter outfile = new BufferedWriter(
                                 new FileWriter("ACMG_IncidentalNoRel.txt", true));         
                   for (int i=0; i<line.length; i++) {
                       outfile.write(line[i] + "\t");
                   }
                   outfile.write(count + "\t" + genotNames + "\n");
              outfile.close();
    }
    

    /*
    public void extractLine(String[] line) throws IOException {
    	for (int i=0; i<header.length; i++) {
            if (line[0].contains(header[i])) {
                WriteToFile(line);
                break;
            }
    	}
    }
    */
    public static void main(String[] args) throws IOException {
           System.out.println("Usage javac Stats.java full/path/to/file /Wellderly/headNorel.txt /Wellderly/headers");
           String InputFile = args[0];
           System.out.println("Input File is " + InputFile);
           String headerFile = args[1];
           System.out.println("Header File is " + headerFile);
           String headerFileFull = args[2];
           System.out.println("Header File full is " + headerFileFull);
           Incidental_perGeno ob1 = new Incidental_perGeno(InputFile, headerFile, headerFileFull);
           ob1.readFile();
    }
}
