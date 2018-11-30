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
import java.util.List;
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
import edu.scripps.yates.utilities.parsers.idparser.IdentificationsParser;
import edu.scripps.yates.utilities.parsers.idparser.IdentifiedProteinGroup;
import edu.scripps.yates.utilities.parsers.idparser.IdentifiedProteinInterface;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class DTASelectParser extends IdentificationsParser<IdentifiedProteinGroup, DTASelectProtein, DTASelectPSM> {
	private static final Logger log = Logger.getLogger(DTASelectParser.class);

	public static final String PROLUCID = "ProLuCID";
	public static final String SEQUEST = "Sequest";
	private static final String UNKNOWN = "Unknown";

	public DTASelectParser(URL u) throws IOException {
		this(u.getFile(), u.openStream());
	}

	public DTASelectParser(String runid, RemoteSSHFileReference s) throws FileNotFoundException {
		this(runid, s.getRemoteInputStream());
	}

	public DTASelectParser(Map<String, RemoteSSHFileReference> s) throws FileNotFoundException {
		super(s);

	}

	public DTASelectParser(List<File> s) throws FileNotFoundException {
		super(s);

	}

	public DTASelectParser(File file) throws FileNotFoundException {
		super(file);

	}

	public DTASelectParser(String runId, InputStream f) {
		super(runId, f);

	}

	@Override
	protected void process(boolean checkFormat) throws IOException {
		final Set<String> psmIds = new THashSet<String>();

		IdentifiedProteinGroup currentProteinGroup = null;
		final TObjectIntHashMap<String> psmHeaderPositions = new TObjectIntHashMap<String>();
		final TObjectIntHashMap<String> proteinHeaderPositions = new TObjectIntHashMap<String>();
		int numDecoy = 0;
		for (final String runId : fs.keySet()) {
			log.info("Reading input stream: " + runId + "...");

			InputStream f = null;
			try {
				f = fs.get(runId);
				currentProteinGroup = new IdentifiedProteinGroup();
				final BufferedInputStream bis = new BufferedInputStream(f);
				final BufferedReader dis = new BufferedReader(new InputStreamReader(bis));

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
							proteinGroups.add(currentProteinGroup);

							// restart the protein group
							currentProteinGroup = new IdentifiedProteinGroup();

						}
						// if (dbIndex == null) {
						DTASelectProtein p = new DTASelectProtein(line, proteinHeaderPositions, ignoreACCFormat);
						if (proteinsByAccession.containsKey(p.getAccession())) {
							final DTASelectProtein p2 = proteinsByAccession.get(p.getAccession());
							p = mergeProteins(p, p2);
						}
						if (!searchEngines.isEmpty()) {
							p.setSearchEngine(searchEngines.iterator().next());
						}
						boolean skip = false;
						if (decoyPattern != null) {
							final Matcher matcher = decoyPattern.matcher(p.getLocus());
							if (matcher.find()) {
								numDecoy++;
								skip = true;
							}
						}
						if (!skip) {
							proteinsByAccession.put(p.getAccession(), p);
							currentProteinGroup.add(p);
						}
						// }
						isPsm = false;
					} else {
						// log.info(line);

						// this is the case of a psm
						isPsm = true;
						if (onlyReadProteins) {
							continue;
						}
						final String psmID = getPSMIdentifier(line, psmHeaderPositions);
						psmIds.add(psmID);
						DTASelectPSM psm;
						if (psmTableByPSMID.containsKey(psmID)) {
							psm = psmTableByPSMID.get(psmID);
						} else {
							psm = new DTASelectPSM(line, psmHeaderPositions, runPath);

							if (!searchEngines.isEmpty()) {
								psm.setSearchEngine(searchEngines.iterator().next());
							}
							if (!spectraFileNames.contains(psm.getRawFileName())) {
								spectraFileNames.add(psm.getRawFileName());
								log.debug(psm.getRawFileName() + " added to a set of " + spectraFileNames.size()
										+ " spectra file names in total");
							}
							final String spectraFileFullPath = new File(runId).getParent() + File.separator
									+ psm.getRawFileName() + ".ms2";
							if (!spectraFileFullPaths.contains(spectraFileFullPath)) {
								spectraFileFullPaths.add(spectraFileFullPath);
								log.debug(spectraFileFullPath + " added to a set of " + spectraFileFullPaths.size()
										+ " spectra paths in total");
							}
							psmTableByPSMID.put(psm.getPsmIdentifier(), psm);
							addToPSMTable(psm);
							// if there is an indexed Fasta, look into it to get
							// the
							// proteins
							if (dbIndex != null) {
								final Set<IndexedProtein> indexedProteins = dbIndex
										.getProteins(psm.getSequence().getSequence());
								if (indexedProteins.isEmpty()) {
									if (!ignoreNotFoundPeptidesInDB) {
										throw new PeptideNotFoundInDBIndexException("The peptide "
												+ psm.getSequence().getSequence()
												+ " is not found in Fasta DB.\nReview the default indexing parameters such as the number of allowed misscleavages.");
									}
								}
								if (indexedProteins != null) {
									log.debug(indexedProteins.size() + " proteins contains "
											+ psm.getSequence().getSequence() + " on fasta file");
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
										DTASelectProtein protein = null;
										if (proteinsByAccession.containsKey(indexedAccession)) {
											protein = proteinsByAccession.get(indexedAccession);
										} else {
											protein = new DTASelectProtein(indexedProtein, ignoreACCFormat);
											proteinsByAccession.put(indexedAccession, protein);
										}
										// add the psm to the protein and
										// vice-versa
										psm.addProtein(protein);
										protein.addPSM(psm);
									}
								}
							}
						}
						if (dbIndex == null || psm.getProteins().isEmpty()) {
							// add the PSM to all the proteins in the current
							// group and all the proteins to the psm
							for (final IdentifiedProteinInterface prot : currentProteinGroup) {
								prot.addPSM(psm);
								psm.addProtein(prot);
							}
						}
						if (checkFormat) {
							// just return with no errors
							return;
						}
					}

				}
				if (dbIndex == null) {
					proteinGroups.add(currentProteinGroup);
					log.info(proteinGroups.size() + " proteins groups read in " + fs.size() + " DTASelect file(s).");
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
					if (f != null) {
						log.info("Closing input stream");
						f.close();
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

	@Override
	protected DTASelectProtein mergeProteins(DTASelectProtein destination, DTASelectProtein origin) {
		if (destination == null)
			return null;

		destination.mergeWithProtein(origin);

		return destination;
	}

	/**
	 * Gets the identifier of the PSM, being: filename-scan
	 *
	 * @param line
	 * @param positions
	 * @return
	 */
	public static String getPSMIdentifier(String line, TObjectIntHashMap<String> positions) {

		final String[] elements = line.split("\t");

		final String psmIdentifier = elements[positions.get(DTASelectPSM.PSM_ID)];
		// log.info("PSM id: " + psmIdentifier);
		final String scan = FastaParser.getScanFromPSMIdentifier(psmIdentifier);
		// log.info("scan number: " + scan);
		final String fileName = FastaParser.getFileNameFromPSMIdentifier(psmIdentifier);

		final Integer charge = FastaParser.getChargeStateFromPSMIdentifier(psmIdentifier);
		// log.info("file name: " + fileName);
		final String fullSequence = elements[positions.get(DTASelectPSM.SEQUENCE)];
		return fileName + "-" + scan + "-" + charge + "-" + fullSequence;
	}

	@Override
	protected String getPSMFullSequence(DTASelectPSM psm) {
		return psm.getFullSequence();
	}

	@Override
	public String getProteinAcc(DTASelectProtein protein) {
		return protein.getLocus();
	}

	@Override
	protected void setPrimaryAccession(String primaryAccession, DTASelectProtein protein) {
		protein.setLocus(primaryAccession);
	}

}
