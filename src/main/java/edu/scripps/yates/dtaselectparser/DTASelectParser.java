/**
 * diego Sep 17, 2013
 */
package edu.scripps.yates.dtaselectparser;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import org.apache.log4j.Logger;

import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.dtaselectparser.util.DTASelectPSM;
import edu.scripps.yates.dtaselectparser.util.DTASelectProtein;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.parsers.idparser.IdentificationsParser;
import edu.scripps.yates.utilities.proteomicsmodel.MSRun;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.factories.MSRunEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.PeptideEx;
import edu.scripps.yates.utilities.proteomicsmodel.staticstorage.StaticProteomicsModelStorage;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class DTASelectParser extends IdentificationsParser {
	private static final Logger log = Logger.getLogger(DTASelectParser.class);

	public static final String PROLUCID = "ProLuCID";
	public static final String SEQUEST = "Sequest";
	private static final String UNKNOWN = "Unknown";

	public DTASelectParser(URL u) throws IOException {
		this(u.getFile(), u.openStream());
	}

	public DTASelectParser(String analysisID, RemoteSSHFileReference s) throws FileNotFoundException {
		this(analysisID, s.getRemoteInputStream());
	}

	public DTASelectParser(Map<String, RemoteSSHFileReference> s) throws FileNotFoundException {
		super(s);

	}

	public DTASelectParser(Collection<File> s) throws FileNotFoundException {
		super(s);

	}

	public DTASelectParser(File file) throws FileNotFoundException {
		super(file);

	}

	public DTASelectParser(String analysisId, InputStream f) {
		super(analysisId, f);

	}

	@Override
	protected void process(boolean checkFormat) throws IOException {
		final Set<String> psmIds = new THashSet<String>();

		ProteinGroup currentProteinGroup = null;
		final TObjectIntHashMap<String> psmHeaderPositions = new TObjectIntHashMap<String>();
		final TObjectIntHashMap<String> proteinHeaderPositions = new TObjectIntHashMap<String>();
		int numDecoy = 0;
		for (final String analysisID : fs.keySet()) {
			log.info("Reading input stream: " + analysisID + "...");
			// if we only get proteins, because the protein's msRuns are coming
			// from the PSMs, we need to assign them a fake msRun that is unique
			// for each input stream
			MSRun uniqueMSRun = null;
			if (onlyReadProteins) {
				uniqueMSRun = new MSRunEx(analysisID, null);
			}
			BufferedReader dis = null;
			try {
				final InputStream f = fs.get(analysisID);
				currentProteinGroup = new ProteinGroup();
				final BufferedInputStream bis = new BufferedInputStream(f);
				dis = new BufferedReader(new InputStreamReader(bis));

				String line;
				int numLine = 0;
				int searchEngineLine = -1;
				boolean intro = false;
				boolean conclusion = false;
				boolean isPsm = false;
				boolean locusStarted = false;
				while ((line = dis.readLine()) != null) {
					numLine++;
					if ("".equals(line.trim())) {
						continue;
					}
					// if (numLine % 10 == 0) {
					// log.info(numLine + " lines read");
					// }
					// log.info(line);
					if (numLine == 1) {
						dtaSelectVersion = line.split(" ")[1];
					}
					if (numLine == 2) {
						runPath = line.trim();
					}
					if (numLine == 3) {
						fastaPath = line.trim();
					}
					if (line.toLowerCase().startsWith(SEQUEST.toLowerCase())) {
						searchEngineLine = numLine;
						searchEngines.add(SEQUEST);
						setSearchEngineVersion(line.split(" ")[1]);
					} else if (line.toLowerCase().startsWith(PROLUCID.toLowerCase())) {
						searchEngineLine = numLine;
						searchEngines.add(PROLUCID);
						setSearchEngineVersion(line.split(" ")[1]);
					} else if (line.toLowerCase().startsWith("?")) {
						searchEngineLine = numLine;
						// searchEngines.add(UNKNOWN);
						// if not known, report as SEQUEST as it is the most
						// common
						searchEngines.add(SEQUEST);
						setSearchEngineVersion(line.split(" ")[1]);
					} else if (searchEngineLine > -1 && numLine >= searchEngineLine + 1 && !locusStarted
							&& !line.startsWith("Locus")) {
						commandLineParameterStrings.add(line);
						commandLineParameter = new DTASelectCommandLineParameters(line);
						if (onlyReadParameters) {
							return;
						}
					} else if (line.startsWith("DTASelect")) {
						// DTASelectProtein p = new DTASelectProtein("DTA", 0,
						// 0, 0,
						// line);
						// ptable.put("DTA", p);
						intro = true;
						continue;
					} else if (line.startsWith("Locus")) {
						locusStarted = true;
						// parse psm header positions
						final String[] splitted = line.split("\t");
						for (int position = 0; position < splitted.length; position++) {
							final String header = splitted[position];
							proteinHeaderPositions.put(header, position);
						}

						// DTASelectProtein p = new DTASelectProtein("Title", 0,
						// 0,
						// 0,
						// line);
						// ptable.put("Title", p);
						intro = false;
						continue;
					} else if (line.startsWith("Unique")) {
						// parse psm header positions
						final String[] splitted = line.split("\t");
						for (int position = 0; position < splitted.length; position++) {
							final String header = splitted[position];
							psmHeaderPositions.put(header, position);
						}

						// DTASelectProtein p = new DTASelectProtein("Unique",
						// 0, 0,
						// 0, line);
						// ptable.put("Unique", p);
						continue;
					}
					if (intro) {
						if (line.contains("Remove subset proteins")) {
							final String[] splitted = line.split("\t");
							if ("TRUE".equals(splitted[0])) {
							}
						}
						// DTASelectProtein p = ptable.get("DTA");
						// p.addExtra(line);
						continue;
					}
					if (conclusion) {
						// DTASelectProtein p = ptable.get("Conclusion");
						// p.addExtra(line);
						continue;
					}

					final String[] elements = line.split("\t");
					// if (elements[1].equals("DTASelectProteins")) {
					if (elements[1].equals("Proteins")) {
						conclusion = true;
						// DTASelectProtein p = new
						// DTASelectProtein("Conclusion",
						// 0,
						// 0, 0, line);
						// ptable.put("Conclusion", p);
						continue;
					}

					// this is the case of a protein
					if (isNumeric(elements[1])) {
						// log.info(line);
						// if comes from a psm line, clear the current group of
						// proteins
						if (isPsm) {
							addProteinGroup(currentProteinGroup);

							// restart the protein group
							currentProteinGroup = new ProteinGroup();

						}
						String acc = line.split("\t")[proteinHeaderPositions.get(DTASelectProtein.LOCUS)];
						if (!ignoreACCFormat) {
							acc = FastaParser.getACC(acc).getAccession();
						}
						boolean skip = false;
						if (decoyPattern != null) {
							final Matcher matcher = decoyPattern.matcher(acc);
							if (matcher.find()) {
								numDecoy++;
								skip = true;
							}
						}
						if (!skip) {
							Protein p;
							if (containsProteinByAccession(acc)) {
								p = getProteinByAccession(acc);
							} else {
								p = new DTASelectProtein(line, proteinHeaderPositions, ignoreACCFormat);
							}
							if (!searchEngines.isEmpty()) {
								p.setSearchEngine(searchEngines.iterator().next());
							}

							addProteinToMapAndList(p);
							currentProteinGroup.add(p);
						}

						isPsm = false;
					} else {
						// log.info(line);

						// this is the case of a psm
						isPsm = true;
						if (onlyReadProteins) {
							// if I only need proteins, because they get the
							// msRuns from their PSMs, we have to add them a
							// fake msRun for the entire collection of proteins
							// here
							for (final GroupableProtein protein : currentProteinGroup) {
								((Protein) protein).addMSRun(uniqueMSRun);
							}
							continue;
						}

						PSM psm = new DTASelectPSM(line, psmHeaderPositions, runPath, isDistinguishModifiedSequences(),
								isChargeSensible());
						if (super.containsPSMByPSMID(psm.getIdentifier())) {
							psm = super.getPSMByPSMID(psm.getIdentifier());
						}
						addPSMToMaps(psm);
						if (!searchEngines.isEmpty()) {
							psm.setSearchEngine(searchEngines.iterator().next());
						}
						if (!spectraFileNames.contains(psm.getMSRun().getRunId())) {
							spectraFileNames.add(psm.getMSRun().getRunId());
							log.debug(psm.getMSRun().getRunId() + " added to a set of " + spectraFileNames.size()
									+ " spectra file names in total");
						}
						final String spectraFileFullPath = new File(analysisID).getParent() + File.separator
								+ psm.getMSRun().getRunId() + ".ms2";
						if (!spectraFileFullPaths.contains(spectraFileFullPath)) {
							spectraFileFullPaths.add(spectraFileFullPath);
							log.debug(spectraFileFullPath + " added to a set of " + spectraFileFullPaths.size()
									+ " spectra paths in total");
						}
						// create the peptide
						Peptide peptide = null;
						final String peptideKey = KeyUtils.getInstance().getSequenceChargeKey(psm,
								isDistinguishModifiedSequences(), isChargeSensible());
						if (StaticProteomicsModelStorage.containsPeptide(psm.getMSRun(), null, peptideKey)) {
							peptide = StaticProteomicsModelStorage.getSinglePeptide(psm.getMSRun(), null, peptideKey);
						} else {
							peptide = new PeptideEx(psm.getFullSequence());
							StaticProteomicsModelStorage.addPeptide(peptide, psm.getMSRun(), null);
							peptide.setSearchEngine(psm.getSearchEngine());
							peptide.addMSRun(psm.getMSRun());
						}

						psm.setPeptide(peptide, true);

						if (dbIndex != null) {
							final Set<IndexedProtein> indexedProteins = dbIndex.getProteins(psm.getSequence());
							if (indexedProteins.isEmpty()) {
								if (!ignoreNotFoundPeptidesInDB) {
									throw new PeptideNotFoundInDBIndexException("The peptide " + psm.getSequence()
											+ " is not found in Fasta DB.\nReview the default indexing parameters such as the number of allowed misscleavages.");
								}
							}
							if (indexedProteins != null) {
								log.debug(indexedProteins.size() + " proteins contains " + psm.getSequence()
										+ " on fasta file");
								for (final IndexedProtein indexedProtein : indexedProteins) {
									final String indexedAccession = indexedProtein.getAccession();
									// we should take into account that in
									// the
									// indexed database you may have decoy
									// hits
									// that you want to avoid
									if (decoyPattern != null) {
										final Matcher matcher = decoyPattern.matcher(indexedAccession);
										if (matcher.find()) {
											numDecoy++;
											continue;
										}
									}
									Protein protein = null;
									if (super.containsProteinByAccession(indexedAccession)) {
										protein = super.getProteinByAccession(indexedAccession);
									} else {
										protein = new DTASelectProtein(indexedProtein, ignoreACCFormat);
										addProteinToMapAndList(protein);
									}
									// add the psm to the protein and
									// vice-versa
									psm.addProtein(protein, true);
								}
							}
						}

						if (dbIndex == null || psm.getProteins().isEmpty()) {
							// add the PSM to all the proteins in the current
							// group and all the proteins to the psm
							for (final GroupableProtein prot : currentProteinGroup) {
								((Protein) prot).addPSM(psm, true);
							}
						}
						if (checkFormat) {
							// just return with no errors
							return;
						}
					}

				}
				if (dbIndex == null) {
					super.addProteinGroup(currentProteinGroup);
					log.info(super.getProteinGroupsNumber() + " proteins groups read in " + fs.size()
							+ " DTASelect file(s).");
				}
			} catch (final IOException e) {
				e.printStackTrace();
				log.error(e.getMessage());
				throw e;
			} catch (final DBIndexStoreException e) {
				e.printStackTrace();
				log.error(e.getMessage());
				throw new IOException(e);
			} finally {
				log.info(numDecoy + " proteins discarded as decoy.");
				try {
					if (dis != null) {
						log.info("Closing input stream");
						dis.close();
						log.info("Input stream closed");

					}
				} catch (final IOException e) {
					e.printStackTrace();
					log.error(e.getMessage());
				}
			}
		}

		processed = true;

	}

	// /**
	// * Gets the identifier of the PSM, being: filename-scan
	// *
	// * @param line
	// * @param positions
	// * @return
	// */
	// public static String getPSMIdentifier(String line,
	// TObjectIntHashMap<String> positions) {
	//
	// final String[] elements = line.split("\t");
	//
	// final String psmIdentifier =
	// elements[positions.get(DTASelectPSM.PSM_ID)];
	// // log.info("PSM id: " + psmIdentifier);
	// final String scan = FastaParser.getScanFromPSMIdentifier(psmIdentifier);
	// // log.info("scan number: " + scan);
	// final String fileName =
	// FastaParser.getFileNameFromPSMIdentifier(psmIdentifier);
	//
	// final Integer charge =
	// FastaParser.getChargeStateFromPSMIdentifier(psmIdentifier);
	// // log.info("file name: " + fileName);
	// final String fullSequence =
	// elements[positions.get(DTASelectPSM.SEQUENCE)];
	// return fileName + "-" + scan + "-" + charge + "-" + fullSequence;
	// }

}
