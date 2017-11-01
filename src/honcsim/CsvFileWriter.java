package honcsim;


import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author ashraf
 * 
 */
public class CsvFileWriter {
	
	//Delimiter used in CSV file
	private static final String COMMA_DELIMITER = ",";
	private static final String NEW_LINE_SEPARATOR = "\n";
	
	//CSV file header
	private static final String FILE_HEADER = "id,firstName,lastName,gender,age";
	public static void main(String[] args){
		String fileName="testcsv.csv";
		writeCsvFile(fileName);
		
	}
	public static void writeCsvFile(String fileName) {
	

		
		
		FileWriter fileWriter = null;
				
		try {
			fileWriter = new FileWriter(fileName);

			//Write the CSV file header
			fileWriter.append(FILE_HEADER.toString());
			
			//Add a new line separator after the header
			fileWriter.append(NEW_LINE_SEPARATOR);
			  int[][] a = {
			            {1, -2, 3}, 
			            {-4, -5, 6, 9}, 
			            {7}, 
			      };
			      
			     for (int i = 0; i < a.length; ++i) {
			        for(int j = 0; j < a[i].length; ++j) {
			           System.out.println(a[i][j]);
			           fileWriter.append("hi");
			           fileWriter.append(COMMA_DELIMITER);
			        }
			        fileWriter.append(NEW_LINE_SEPARATOR);
			     }
			//Write a new student object list to the CSV file
			

			
			
			System.out.println("CSV file was created successfully !!!");
			
		} catch (Exception e) {
			System.out.println("Error in CsvFileWriter !!!");
			e.printStackTrace();
		} finally {
			
			try {
				fileWriter.flush();
				fileWriter.close();
			} catch (IOException e) {
				System.out.println("Error while flushing/closing fileWriter !!!");
                e.printStackTrace();
			}
			
		}
	}
}