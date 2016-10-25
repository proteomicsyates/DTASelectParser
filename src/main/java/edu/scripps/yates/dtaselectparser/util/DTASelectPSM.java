package edu.scripps.yates.dtaselectparser.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.utilities.fasta.FastaParser;

public class DTASelectPSM {
	private static final Logger log = Logger.getLogger(DTASelectParser.class);
	public static final Map<String, DTASelectPSM> map = new HashMap<String, DTASelectPSM>();
	public static final String PSM_ID = "FileName";
	private static final String XCORR = "XCorr";
	private static final String DELTACN = "DeltCN";
	private static final String PROB = "Prob%";
	private static final String CONF = "Conf%";
	private static final String MH = "M+H+";
	private static final String CALC_MH = "CalcM+H+";
	private static final String TOTAL_INTENSITY = "TotalIntensity";
	private static final String SPR = "SpR";
	private static final String PPM = "PPM";
	private static final String PI = "pI";
	private static final String PROB_SCORE = "Prob Score";
	private static final String RT = "RT";
	private static final String ION_PROPORTION = "IonProportion";
	private static final String REDUNDANCY = "Redundancy";
	public static final String SEQUENCE = "Sequence";
	private final String rawPSMIdentifier;
	private final Double xcorr;
	private final Double deltacn;
	private final Double mh;
	private final Double prob;
	private final Double conf;
	private final Double calcMh;
	private final Double totalIntensity;
	private final Integer spr;
	private final Double prob_score;
	private final Double ionProportion;
	private final Double ppmError;
	private final Integer redundancy;
	private final DTASelectPeptideSequence sequence;
	private final Double pi;
	private final String fullSequence;
	private final String scan;
	private final Set<DTASelectProtein> proteins = new HashSet<DTASelectProtein>();
	private final String runID;
	private final String runPath;
	private String searchEngine;
	private final String fileName;
	private final Integer chargeState;
	private String psmIdentifier;
	private String msRunId;
	private Double rtInMinutes;

	public DTASelectPSM(String dtaSelectRow, HashMap<String, Integer> positions, String runPath) {
		// log.info("Creating PSM from line: " + dtaSelectRow);

		this.runPath = runPath;
		// parse the headerRow
		String[] elements = dtaSelectRow.split("\t");
		rawPSMIdentifier = elements[positions.get(PSM_ID)];
		scan = FastaParser.getScanFromPSMIdentifier(rawPSMIdentifier);
		fileName = FastaParser.getFileNameFromPSMIdentifier(rawPSMIdentifier);
		chargeState = FastaParser.getChargeStateFromPSMIdentifier(rawPSMIdentifier);
		runID = fileName;
		// store by scan number in the map
		map.put(scan, this);

		xcorr = Double.parseDouble(elements[positions.get(XCORR)]);
		deltacn = Double.parseDouble(elements[positions.get(DELTACN)]);
		conf = Double.parseDouble(elements[positions.get(CONF)]);
		mh = Double.parseDouble(elements[positions.get(MH)]);
		calcMh = Double.parseDouble(elements[positions.get(CALC_MH)]);
		totalIntensity = Double.valueOf(elements[positions.get(TOTAL_INTENSITY)]);
		spr = Integer.valueOf(elements[positions.get(SPR)]);

		ionProportion = Double.parseDouble(elements[positions.get(ION_PROPORTION)]);
		redundancy = Integer.valueOf(elements[positions.get(REDUNDANCY)]);

		fullSequence = elements[positions.get(SEQUENCE)];
		sequence = new DTASelectPeptideSequence(fullSequence);

		if (positions.containsKey(PROB))
			prob = Double.valueOf(elements[positions.get(PROB)]);
		else
			prob = null;
		if (positions.containsKey(PPM))
			ppmError = Double.parseDouble(elements[positions.get(PPM)]);
		else
			ppmError = null;
		if (positions.containsKey(PI))
			pi = Double.parseDouble(elements[positions.get(PI)]);
		else
			pi = null;
		if (positions.containsKey(PROB_SCORE))
			prob_score = Double.parseDouble(elements[positions.get(PROB_SCORE)]);
		else
			prob_score = null;

		if (positions.containsKey(RT)) {
			rtInMinutes = Double.parseDouble(elements[positions.get(RT)]);
		} else {
			rtInMinutes = null;
		}
	}

	public String getSpectraFileName() {
		return fileName;

	}

	/**
	 * @return the psmIdentifier
	 */
	public String getPsmIdentifier() {
		if (psmIdentifier == null) {
			psmIdentifier = KeyUtils.getSpectrumKey(this, true);
		}
		return psmIdentifier;
	}

	public String getRawPSMID() {
		return rawPSMIdentifier;
	}

	/**
	 * @return the xcorr
	 */
	public Double getXcorr() {
		return xcorr;
	}

	/**
	 * @return the deltacn
	 */
	public Double getDeltacn() {
		return deltacn;
	}

	/**
	 * @return the mh
	 */
	public Double getMh() {
		return mh;
	}

	/**
	 * @return the prob
	 */
	public Double getProb() {
		return prob;
	}

	/**
	 * @return the conf
	 */
	public Double getConf() {
		return conf;
	}

	/**
	 * @return the calcMh
	 */
	public Double getCalcMh() {
		return calcMh;
	}

	/**
	 * @return the totalIntensity
	 */
	public Double getTotalIntensity() {
		return totalIntensity;
	}

	/**
	 * @return the spr
	 */
	public Integer getSpr() {
		return spr;
	}

	/**
	 * @return the prob_score
	 */
	public Double getProb_score() {
		return prob_score;
	}

	/**
	 * @return the ionProportion
	 */
	public Double getIonProportion() {
		return ionProportion;
	}

	/**
	 * @return the redundancy
	 */
	public Integer getRedundancy() {
		return redundancy;
	}

	/**
	 * @return the sequence
	 */
	public DTASelectPeptideSequence getSequence() {
		return sequence;
	}

	public List<DTASelectModification> getModifications() {
		return sequence.getModifications();
	}

	/**
	 * @return the ppmError
	 */
	public Double getPpmError() {
		return ppmError;
	}

	/**
	 * @return the pi
	 */
	public Double getPi() {
		return pi;
	}

	public String getFullSequence() {
		return fullSequence;
	}

	/**
	 * @return the scan
	 */
	public String getScan() {
		return scan;
	}

	/**
	 * @return the proteins
	 */
	public Set<DTASelectProtein> getProteins() {
		return proteins;
	}

	public void addProtein(DTASelectProtein protein) {
		proteins.add(protein);
	}

	/**
	 * @return the runID
	 */
	public String getRunID() {
		return runID;
	}

	/**
	 * @return the runPath
	 */
	public String getRunPath() {
		return runPath;
	}

	public void setSearchEngine(String searchEngine) {
		this.searchEngine = searchEngine;
	}

	/**
	 * @return the searchEngine
	 */
	public String getSearchEngine() {
		return searchEngine;
	}

	/**
	 *
	 * Parse the scan number to get the last number that whould be the charge:
	 * <br>
	 * brain2dayAHA092613s09.1849.1849.2 is charge 2
	 *
	 *
	 * @return the chargeState
	 */
	public Integer getChargeState() {
		return chargeState;
	}

	public String getMsRunId() {
		if (msRunId == null) {
			return getSpectraFileName();
		}
		return msRunId;
	}

	/**
	 * @param msRunId
	 *            the msRunId to set
	 */
	public void setMsRunId(String msRunId) {
		this.msRunId = msRunId;
	}

	public Double getRTInMin() {
		return rtInMinutes;
	}
}
