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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.lang.math.NumberUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.IndexedProtein;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.dtaselectparser.util.DTASelectPSM;
import edu.scripps.yates.dtaselectparser.util.DTASelectProtein;
import edu.scripps.yates.dtaselectparser.util.DTASelectProteinGroup;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.files.Parser;
import edu.scripps.yates.utilities.ipi.IPI2UniprotACCMap;
import edu.scripps.yates.utilities.model.enums.AccessionType;
import edu.scripps.yates.utilities.model.factories.AccessionEx;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class DTASelectParser implements Parser {
	private static final Logger log = Logger.getLogger(DTASelectParser.class);
	private final Map<String, DTASelectProtein> proteinsByAccession = new THashMap<String, DTASelectProtein>();
	private final Map<String, DTASelectPSM> psmTableByPSMID = new THashMap<String, DTASelectPSM>();

	private final Map<String, Set<DTASelectPSM>> psmTableByFullSequence = new THashMap<String, Set<DTASelectPSM>>();
	private final Set<DTASelectProteinGroup> dtaSelectProteinGroups = new THashSet<DTASelectProteinGroup>();
	private final List<String> keys = new ArrayList<String>();
	public static final String PROLUCID = "ProLuCID";
	public static final String SEQUEST = "Sequest";
	private static final String UNKNOWN = "Unknown";
	private String runPath;
	private DBIndexInterface dbIndex;
	private boolean processed = false;
	private final Map<String, InputStream> fs;
	private Pattern decoyPattern;
	private final List<String> commandLineParameterStrings = new ArrayList<String>();
	private final Set<String> searchEngines = new THashSet<String>();
	private DTASelectCommandLineParameters commandLineParameter;
	private String searchEngineVersion;
	private String fastaPath;
	private final Set<String> spectraFileNames = new THashSet<String>();
	private String dtaSelectVersion;
	private final Set<String> spectraFileFullPaths = new THashSet<String>();
	private UniprotProteinLocalRetriever uplr;
	private String uniprotVersion;
	private boolean ignoreNotFoundPeptidesInDB;

	public DTASelectParser(URL u) throws IOException {
		this(u.getFile(), u.openStream());
	}

	public DTASelectParser(String runid, RemoteSSHFileReference s) throws FileNotFoundException {
		this(runid, s.getRemoteInputStream());
	}

	public DTASelectParser(Map<String, RemoteSSHFileReference> s) throws FileNotFoundException {
		fs = new THashMap<String, InputStream>();
		for (final String key : s.keySet()) {
			final RemoteSSHFileReference server = s.get(key);
			// final File remoteFile = server.getRemoteFile();
			fs.put(key, server.getRemoteInputStream());
		}

	}

	public DTASelectParser(List<File> s) throws FileNotFoundException {
		fs = new THashMap<String, InputStream>();
		for (final File remoteFile : s) {
			fs.put(remoteFile.getAbsolutePath(), FileUtils.getInputStream(remoteFile));
		}

	}

	public DTASelectParser(File file) throws FileNotFoundException {
		log.debug("Beggining of constructor with file " + file.getAbsolutePath());
		fs = new THashMap<String, InputStream>();
		fs.put(file.getAbsolutePath(), FileUtils.getInputStream(file));
		log.debug("end of constructor");
	}

	public DTASelectParser(String runId, InputStream f) {
		fs = new THashMap<String, InputStream>();
		fs.put(runId, f);

	}

	private void process() throws IOException {
		process(false);
	}

	private void process(boolean checkFormat) throws IOException {
		final Set<String> psmIds = new THashSet<String>();

		DTASelectProteinGroup currentProteinGroup = null;
		final TObjectIntHashMap<String> psmHeaderPositions = new TObjectIntHashMap<String>();
		final TObjectIntHashMap<String> proteinHeaderPositions = new TObjectIntHashMap<String>();
		int numDecoy = 0;
		for (final String runId : fs.keySet()) {
			log.info("Reading input stream: " + runId + "...");

			InputStream f = null;
			try {
				f = fs.get(runId);
				currentProteinGroup = new DTASelectProteinGroup();
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
						if (isPsm && dbIndex == null) {
							dtaSelectProteinGroups.add(currentProteinGroup);

							// restart the protein group
							currentProteinGroup = new DTASelectProteinGroup();

						}
						// if (dbIndex == null) {
						DTASelectProtein p = new DTASelectProtein(line, proteinHeaderPositions);
						if (proteinsByAccession.containsKey(p.getAccession())) {
							final DTASelectProtein p2 = proteinsByAccession.get(p.getAccession());
							p = mergeProteins(p, p2);
						}
						if (!searchEngines.isEmpty())
							p.setSearchEngine(searchEngines.iterator().next());
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
							keys.add(p.getAccession());
							currentProteinGroup.add(p);
						}
						// }
						isPsm = false;
					} else {
						// log.info(line);

						// this is the case of a psm
						isPsm = true;
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
										keys.add(indexedAccession);
										DTASelectProtein protein = null;
										if (proteinsByAccession.containsKey(indexedAccession)) {
											protein = proteinsByAccession.get(indexedAccession);
										} else {
											protein = new DTASelectProtein(indexedProtein);
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
						if (dbIndex == null) {
							// add the PSM to all the proteins in the current
							// group
							// and all the proteins to the psm
							for (final DTASelectProtein prot : currentProteinGroup) {
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
					dtaSelectProteinGroups.add(currentProteinGroup);
					log.info(dtaSelectProteinGroups.size() + " proteins groups read in " + fs.size()
							+ " DTASelect file(s).");
				}
			} catch (final Exception e) {
				e.printStackTrace();
				log.error(e.getMessage());
				throw e;
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

	private DTASelectProtein mergeProteins(DTASelectProtein destination, DTASelectProtein origin) {
		if (destination == null)
			return null;

		destination.mergeWithProtein(origin);

		return destination;
	}

	private void addToPSMTable(DTASelectPSM psm) {
		final String psmKey = psm.getFullSequence();
		if (psmTableByFullSequence.containsKey(psmKey)) {
			psmTableByFullSequence.get(psmKey).add(psm);
		} else {
			final Set<DTASelectPSM> set = new THashSet<DTASelectPSM>();
			set.add(psm);
			psmTableByFullSequence.put(psmKey, set);
		}

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

	public Map<String, DTASelectProtein> getDTASelectProteins() throws IOException {
		if (!processed)
			startProcess();

		return proteinsByAccession;
	}

	public Set<DTASelectProteinGroup> getProteinGroups() throws IOException {
		if (dbIndex != null)
			throw new IllegalArgumentException(
					"Reading proteins with a FASTA database, will not result in Protein groups");
		if (!processed)
			startProcess();
		return dtaSelectProteinGroups;
	}

	public Map<String, DTASelectPSM> getDTASelectPSMsByPSMID() throws IOException {
		if (!processed)
			startProcess();
		return psmTableByPSMID;
	}

	public Map<String, Set<DTASelectPSM>> getDTASelectPSMsByFullSequence() throws IOException {
		if (!processed)
			startProcess();
		return psmTableByFullSequence;
	}

	private boolean isNumeric(String string) {
		return NumberUtils.isNumber(string);

	}

	/**
	 * @param dbIndex
	 *            the dbIndex to set
	 */
	public void setDbIndex(DBIndexInterface dbIndex) {
		this.dbIndex = dbIndex;
	}

	public void setDecoyPattern(String patternString) throws PatternSyntaxException {
		if (patternString != null) {
			decoyPattern = Pattern.compile(patternString);
		}
	}

	/**
	 * @return the runPath
	 * @throws IOException
	 */
	public String getRunPath() throws IOException {
		if (!processed)
			startProcess();
		return runPath;
	}

	/**
	 * @return the commandLineParameterStrings
	 * @throws IOException
	 */
	public List<String> getCommandLineParameterStrings() throws IOException {
		if (!processed)
			startProcess();
		return commandLineParameterStrings;
	}

	/**
	 * @return the searchEngines
	 * @throws IOException
	 */
	public Set<String> getSearchEngines() throws IOException {
		if (!processed)
			startProcess();
		return searchEngines;
	}

	/**
	 * @return the commandLineParameter
	 * @throws IOException
	 */
	public DTASelectCommandLineParameters getCommandLineParameter() throws IOException {
		if (!processed)
			startProcess();
		return commandLineParameter;
	}

	/**
	 * @return the searchEngineVersion
	 * @throws IOException
	 */
	public String getSearchEngineVersion() throws IOException {
		if (!processed)
			startProcess();

		return searchEngineVersion;
	}

	/**
	 * @param searchEngineVersion
	 *            the searchEngineVersion to set
	 */
	public void setSearchEngineVersion(String searchEngineVersion) {
		this.searchEngineVersion = searchEngineVersion;
	}

	public String getFastaPath() throws IOException {
		if (!processed)
			startProcess();

		return fastaPath;
	}

	public Set<String> getSpectraFileNames() throws IOException {
		if (!processed)
			startProcess();

		return spectraFileNames;
	}

	public String getDecoyPattern() throws IOException {
		if (!processed)
			startProcess();

		if (decoyPattern != null) {
			return decoyPattern.toString();
		}
		return null;
	}

	public String getDTASelectVersion() throws IOException {
		if (!processed)
			startProcess();

		return dtaSelectVersion;
	}

	public Set<String> getSpectraFileFullPaths() throws IOException {
		if (!processed)
			startProcess();
		return spectraFileFullPaths;
	}

	private void mergeProteinsWithSecondaryAccessionsInParser() throws IOException {
		if (uplr == null) {
			return;
		}
		final Set<String> accessions = new THashSet<String>();
		final Map<String, String> accToLocus = new THashMap<String, String>();
		final Map<String, DTASelectProtein> dtaSelectProteins = getDTASelectProteins();
		for (final DTASelectProtein protein : dtaSelectProteins.values()) {
			final String accession = FastaParser.getACC(protein.getLocus()).getFirstelement();
			accessions.add(accession);
			accToLocus.put(accession, protein.getLocus());
		}
		String latestVersion = "latestVersion";
		if (uniprotVersion != null) {
			latestVersion = "version " + uniprotVersion;
		}
		// split into chunks of 500 accessions in order to show progress
		final int chunckSize = 500;
		final List<Set<String>> listOfSets = new ArrayList<Set<String>>();
		Set<String> set = new THashSet<String>();
		for (final String accession : accessions) {
			set.add(accession);
			if (set.size() == chunckSize) {
				listOfSets.add(set);
				set = new THashSet<String>();
			}
		}
		listOfSets.add(set);

		int numObsoletes = 0;
		log.info("Merging proteins that have secondary accessions according to Uniprot " + latestVersion + "...");

		final int initialSize = accessions.size();
		final ProgressCounter counter = new ProgressCounter(initialSize, ProgressPrintingType.PERCENTAGE_STEPS, 0);
		for (final Set<String> accessionSet : listOfSets) {
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(uniprotVersion, accessionSet);
			for (final String accession : accessionSet) {
				counter.increment();
				final String progress = counter.printIfNecessary();
				if (!"".contentEquals(progress)) {
					log.info(progress);
				}
				final DTASelectProtein protein = dtaSelectProteins.get(accToLocus.get(accession));
				final Entry entry = annotatedProteins.get(accession);
				if (entry != null && entry.getAccession() != null && !entry.getAccession().isEmpty()) {
					final String primaryAccession = entry.getAccession().get(0);
					if (!accession.equals(primaryAccession) && !accession.contains(primaryAccession)) {
						log.info("Replacing Uniprot accession " + accession + " by " + primaryAccession);
						protein.setLocus(primaryAccession);
						if (dtaSelectProteins.containsKey(primaryAccession)) {
							// there was already a protein with that
							// primaryAccession
							final DTASelectProtein quantifiedProtein2 = dtaSelectProteins.get(primaryAccession);
							// merge quantifiedPRotein and quantifiedPRotein2
							mergeProteins(protein, quantifiedProtein2);

						} else {
							numObsoletes++;
						}
						// remove old/secondary accession
						dtaSelectProteins.remove(accession);
						dtaSelectProteins.put(primaryAccession, protein);

					}
				} else {
					// // remove the protein because is obsolete
					// log.info(quantifiedProtein.getAccession());
					// parser.getProteinMap().remove(accession);
				}
			}
		}
		final int finalSize = dtaSelectProteins.size();
		if (initialSize != finalSize) {
			log.info(finalSize - initialSize
					+ " proteins with secondary accessions were merged with the corresponding protein with primary accession");
		}
		log.info("Obsolete accessions from " + numObsoletes + " proteins were changed to primary ones");
	}

	/**
	 * To be called after process().<br>
	 * If proteins have IPI accessions, look for the mapping from IPI 2 Uniprot.
	 * It adds new entries to the map, but it doesn't create any new
	 * {@link QuantifiedProteinInterface}
	 */
	private void mapIPI2Uniprot() {
		if (!proteinsByAccession.isEmpty()) {
			final int originalNumberOfEntries = proteinsByAccession.size();
			final Map<String, DTASelectProtein> newMap = new THashMap<String, DTASelectProtein>();
			for (final String accession : proteinsByAccession.keySet()) {

				final Pair<String, String> acc = FastaParser.getACC(accession);
				if (acc.getSecondElement().equals("IPI")) {
					final DTASelectProtein protein = proteinsByAccession.get(accession);
					Accession primaryAccession = new AccessionEx(accession, AccessionType.IPI);
					final Pair<Accession, Set<Accession>> pair = IPI2UniprotACCMap.getInstance()
							.getPrimaryAndSecondaryAccessionsFromIPI(primaryAccession);
					if (pair.getFirstelement() != null) {
						primaryAccession = pair.getFirstelement();
						if (!newMap.containsKey(primaryAccession)) {
							newMap.put(primaryAccession.getAccession(), protein);
						}
					}
					final Set<Accession> secondaryAccs = pair.getSecondElement();
					if (secondaryAccs != null) {
						for (final Accession secondaryAcc : secondaryAccs) {
							if (!newMap.containsKey(secondaryAcc.getAccession())) {
								newMap.put(secondaryAcc.getAccession(), protein);
							}
						}

					}
				}
			}
			for (final String acc : newMap.keySet()) {
				if (!proteinsByAccession.containsKey(acc)) {
					proteinsByAccession.put(acc, newMap.get(acc));
				}
			}
			if (originalNumberOfEntries != proteinsByAccession.size()) {
				log.info("Protein Map expanded from " + originalNumberOfEntries + " to " + proteinsByAccession.size());
			}
		}
	}

	private void startProcess() throws IOException {
		// first process
		process();
		// remove psms assigned to decoy proteins that were discarded
		removeDecoyPSMs();
		// second expand protein map
		mapIPI2Uniprot();
		// third merge proteins with secondary accessions
		mergeProteinsWithSecondaryAccessionsInParser();

		log.info(proteinsByAccession.size() + " proteins read in " + fs.size() + " DTASelect file(s).");
		log.info(psmTableByFullSequence.size() + " peptides read in " + fs.size() + " DTASelect file(s).");
		log.info(psmTableByPSMID.size() + " psms read in " + fs.size() + " DTASelect file(s).");

	}

	private void removeDecoyPSMs() {
		if (decoyPattern != null) {
			// in case of decoyPattern is enabled, we may have some PSMs
			// assigned to
			// those decoy proteins that have not been saved,
			// so we need to discard them
			// We iterate over the psms, and we will remove the ones with no
			// proteins
			final Set<String> psmIdsToDelete = new THashSet<String>();
			for (final String psmID : psmTableByPSMID.keySet()) {
				if (psmTableByPSMID.get(psmID).getProteins().isEmpty()) {
					psmIdsToDelete.add(psmID);
				}
			}
			log.info("Removing " + psmIdsToDelete.size() + " PSMs assigned to decoy discarded proteins");
			for (final String psmID : psmIdsToDelete) {
				final DTASelectPSM dtaSelectPSM = psmTableByPSMID.get(psmID);
				if (!dtaSelectPSM.getProteins().isEmpty()) {
					throw new IllegalArgumentException("This should not happen");
				}
				final Set<DTASelectPSM> set = psmTableByFullSequence.get(dtaSelectPSM.getFullSequence());
				final boolean removed = set.remove(dtaSelectPSM);
				if (!removed) {
					throw new IllegalArgumentException("This should not happen");
				}
				if (set.isEmpty()) {
					// remove the entry by full sequence
					psmTableByFullSequence.remove(dtaSelectPSM.getFullSequence());
				}
				// remove psmTableByPsmID
				psmTableByPSMID.remove(psmID);
			}

			log.info(psmIdsToDelete.size() + " PSMs discarded as decoy");
		}
	}

	/**
	 * @return the ignoreNotFoundPeptidesInDB
	 */
	public boolean isIgnoreNotFoundPeptidesInDB() {
		return ignoreNotFoundPeptidesInDB;
	}

	/**
	 * @param ignoreNotFoundPeptidesInDB
	 *            the ignoreNotFoundPeptidesInDB to set
	 */
	public void setIgnoreNotFoundPeptidesInDB(boolean ignoreNotFoundPeptidesInDB) {
		this.ignoreNotFoundPeptidesInDB = ignoreNotFoundPeptidesInDB;
	}

	public void enableProteinMergingBySecondaryAccessions(UniprotProteinLocalRetriever uplr, String uniprotVersion) {
		this.uplr = uplr;
		this.uniprotVersion = uniprotVersion;
	}

	public Set<String> getInputFilePathes() {
		return fs.keySet();
	}

	public DBIndexInterface getDBIndex() {
		return dbIndex;
	}

	public Set<String> getUniprotAccSet() {
		final Set<String> ret = new THashSet<String>();
		try {
			final Set<String> keySet = getDTASelectProteins().keySet();
			for (final String acc : keySet) {
				final String uniProtACC = FastaParser.getUniProtACC(acc);
				if (uniProtACC != null) {
					ret.add(uniProtACC);
				}
			}
		} catch (final IOException e) {
			e.printStackTrace();
		}

		return ret;
	}

	@Override
	public boolean canRead(File file) {
		try {
			fs.put("TMP***", FileUtils.getInputStream(file));
			if (!processed) {
				process(true);
			}
			if (!proteinsByAccession.isEmpty()) {
				if (!psmTableByPSMID.isEmpty()) {
					return true;
				}
			}
			return false;
		} catch (final Exception e) {
			return false;
		} finally {
			fs.remove("TMP***");
		}
	}
}
