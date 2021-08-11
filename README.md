# DTASelectParser
Java parser for DTASelect output files.
  
Download the latest version from [here](https://github.com/proteomicsyates/DTASelectParser/releases/latest)  

Alternatively, you can get it from our **maven repository**:  
  
**Dependency:**  
```
<dependency>  
  <groupId>edu.scripps.yates</groupId>  
  <artifactId>dtaselectparser</artifactId>  
  <version>1.1.2-SNAPSHOT</version>  
</dependency>  
```  
  
**Maven repositories:**  
 ```
<repository>  
  <id>internal</id>  
  <url>http://sealion.scripps.edu/archiva/repository/internal/</url>  
</repository>  
<snapshotRepository>  
  <id>snapshots</id>  
  <url>http://sealion.scripps.edu/archiva/repository/snapshots/</url>  
</snapshotRepository>  
```

**Usage example:**  
```
// dtaSelect File object
File dtaSelectFile = ...

// Create the parser
DTASelectParser parser = new DTASelectParser(dtaSelectFile.toURL());

// get number of proteins
final int numProteins = parser.getDTASelectProteins().size();

// get number of PSMs
final int numPSMs = parser.getDTASelectPSMsByPSMID().size();

// get number of peptides
final int numPeptides = parser.getDTASelectPSMsBySequence().size();

// get number of protein groups
final int numProteinGroups = parser.getProteinGroups().size();

// get a particular protein
final DTASelectProtein dtaSelectProtein = parser.getDTASelectProteins().get("P10120");

// print accession
System.out.println(dtaSelectProtein.getLocus());

// get PSMs associated to that protein and print the number
final List<DTASelectPSM> psMs = dtaSelectProtein.getPSMs();
System.out.println(psMs.size() + " PSMS");

// iterate over the psms
for (DTASelectPSM dtaSelectPSM : psMs) {

  // print the psm ID and the DeltaCn and Xcorr scores
	System.out.println(dtaSelectPSM.getPsmIdentifier() + " " + dtaSelectPSM.getDeltacn() + " " + dtaSelectPSM.getXcorr());
	
  // print the peptide sequence with modifications (if available), without modifications and 
	System.out.println(dtaSelectPSM.getSequence().getRawSequence() + " "	+ dtaSelectPSM.getSequence().getSequence());
	
}

```
