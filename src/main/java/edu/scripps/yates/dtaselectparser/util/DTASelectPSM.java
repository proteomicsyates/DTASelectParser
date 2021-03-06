package edu.scripps.yates.dtaselectparser.util;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.parsers.idparser.PeptideSequence;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractPSM;
import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.factories.MSRunEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.ScoreEx;
import edu.scripps.yates.utilities.proteomicsmodel.staticstorage.StaticProteomicsModelStorage;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.staticstorage.StaticStrings;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class DTASelectPSM extends AbstractPSM {
	private static final Logger log = Logger.getLogger(DTASelectParser.class);
	public static final Map<String, DTASelectPSM> map = new THashMap<String, DTASelectPSM>();
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
	public static final String XIC = "XIC";
	public static final String ESTIMATED_XIC = "Estimated_XIC";
	public static final String CCS = "CCS";
	private final Float prob;
	private final Float conf;
	private final Float prob_score;
	private final Integer redundancy;
	private final PeptideSequence peptideSequence;
	private final String rawFileName;
	private Float xic;
	private Float estimatedXIC;
	private Float ccs;
	private final static String scoreType = "PSM-level identification statistic";

	public DTASelectPSM(String dtaSelectRow, TObjectIntHashMap<String> positions, String runPath,
			boolean distinguishModifiedSequence, boolean chargeStateSensible) {
		// log.info("Creating PSM from line: " + dtaSelectRow);
		super(distinguishModifiedSequence, chargeStateSensible);
		// parse the headerRow
		final String[] elements = dtaSelectRow.split("\t");
		final String rawPSMIdentifier = elements[positions.get(PSM_ID)];
		setScanNumber(FastaParser.getScanFromPSMIdentifier(rawPSMIdentifier));
		rawFileName = FastaParser.getFileNameFromPSMIdentifier(rawPSMIdentifier);

		if (StaticProteomicsModelStorage.containsMSRun(rawFileName)) {
			setMSRun(StaticProteomicsModelStorage.getMSRun(rawFileName));
			getMSRun().setPath(runPath);
		} else {
			setMSRun(new MSRunEx(rawFileName, runPath));
			StaticProteomicsModelStorage.addMSRun(getMSRun());
		}
		setChargeState(FastaParser.getChargeStateFromPSMIdentifier(rawPSMIdentifier));

		// store by scan number in the map
		map.put(getScanNumber(), this);

		setXCorr(Float.parseFloat(elements[positions.get(XCORR)]));
		setDeltaCn(Float.parseFloat(elements[positions.get(DELTACN)]));
		conf = Float.parseFloat(elements[positions.get(CONF)]);
		setExperimentalMH(Float.parseFloat(elements[positions.get(MH)]));
		setCalcMH(Float.parseFloat(elements[positions.get(CALC_MH)]));
		super.setTotalIntensity(Float.valueOf(elements[positions.get(TOTAL_INTENSITY)]));
		setSpr(Integer.valueOf(elements[positions.get(SPR)]));

		setIonProportion(Float.parseFloat(elements[positions.get(ION_PROPORTION)]));
		redundancy = Integer.valueOf(elements[positions.get(REDUNDANCY)]);
		peptideSequence = new PeptideSequence(StaticStrings.getUniqueInstance(elements[positions.get(SEQUENCE)]), true);

		if (positions.containsKey(PROB)) {
			prob = Float.valueOf(elements[positions.get(PROB)]);
		} else {
			prob = null;
		}
		if (positions.containsKey(PPM)) {
			setMassErrorPPM(Float.parseFloat(elements[positions.get(PPM)]));
		}
		if (positions.containsKey(PI)) {
			setPi(Float.parseFloat(elements[positions.get(PI)]));
		}
		if (positions.containsKey(PROB_SCORE)) {
			prob_score = Float.parseFloat(elements[positions.get(PROB_SCORE)]);
		} else {
			prob_score = null;
		}

		if (positions.containsKey(RT)) {
			setRtInMinutes(Float.parseFloat(elements[positions.get(RT)]));
		}
		// XIC added on Dec2020
		if (positions.containsKey(XIC)) {
			xic = Float.parseFloat(elements[positions.get(XIC)]);
		} else {
			xic = null;
		}
		if (positions.containsKey(ESTIMATED_XIC)) {
			estimatedXIC = Float.parseFloat(elements[positions.get(ESTIMATED_XIC)]);
		} else {
			estimatedXIC = null;
		}
		if (positions.containsKey(CCS)) {
			ccs = Float.parseFloat(elements[positions.get(CCS)]);
		} else {
			ccs = null;
		}

		setIdentifier(
				KeyUtils.getInstance().getSpectrumKey(this, isDistinguishModifiedSequence(), isChargeStateSensible()));
		setKey(getIdentifier());
	}

	/**
	 * @return the prob
	 */
	public Float getProb() {
		return prob;
	}

	/**
	 * @return the conf
	 */
	public Float getConf() {
		return conf;
	}

	/**
	 * @return the prob_score
	 */
	public Float getProb_score() {
		return prob_score;
	}

	/**
	 * @return the redundancy
	 */
	public Integer getRedundancy() {
		return redundancy;
	}

	@Override
	public List<PTM> getPTMs() {
		return peptideSequence.getModifications();
	}

	@Override
	public String getFullSequence() {
		return peptideSequence.getFullSequence();
	}

	@Override
	public String getSequence() {
		return peptideSequence.getSequence();
	}

	@Override
	public String getAfterSeq() {
		return String.valueOf(peptideSequence.getAfterSeq());
	}

	@Override
	public String getBeforeSeq() {
		return String.valueOf(peptideSequence.getBeforeSeq());
	}

	@Override
	public Set<Score> getScores() {
		if (super.getScores().isEmpty()) {
			// add xcorr and deltacn
			addScore(new ScoreEx(String.valueOf(getXCorr()), "XCorr", scoreType, "XCorr"));
			addScore(new ScoreEx(String.valueOf(getDeltaCn()), "DeltaCN", scoreType, "DeltaCN"));
			// add CCS and XIC and Estimated_XIC as scores
			if (xic != null) {
				addScore(new ScoreEx(String.valueOf(xic), XIC, "quantification datatype",
						"Area of the extracted ion chromatogram"));
			}
			if (estimatedXIC != null) {
				addScore(new ScoreEx(String.valueOf(estimatedXIC), ESTIMATED_XIC, "quantification datatype",
						"Estimated area of the extracted ion chromatogram"));
			}
			if (ccs != null) {
				addScore(new ScoreEx(String.valueOf(ccs), CCS, "collisional cross sectional area", ""));
			}
		}
		return super.getScores();
	}

}
