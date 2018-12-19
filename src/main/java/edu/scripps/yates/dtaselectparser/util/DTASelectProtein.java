package edu.scripps.yates.dtaselectparser.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.proteomicsmodel.AbstractProtein;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.MSRun;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.GeneEx;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class DTASelectProtein extends AbstractProtein {
	private static final Logger log = Logger.getLogger(DTASelectProtein.class);
	public static final String LOCUS = "Locus";
	private static final String SP_COUNT = "Spectrum Count";
	private static final String COVERAGE = "Sequence Coverage";
	private static final String LENGTH = "Length";
	private static final String MW = "MolWt";
	private static final String PI = "pI";
	private static final String DESCRIPTION = "Descriptive Name";
	private static final String NSAF = "NSAF";
	private static final String EMPAI = "EMPAI";
	private String locus;

	// whether to use locus or to parse the locus from input dtaselect file to
	// get the accession (which may be more time consuming)
	private final boolean ignoreACCFormat;

	public DTASelectProtein(String lineToParse, TObjectIntHashMap<String> positions, boolean ignoreACCFormat) {
		// log.info("Creating protein from line: " + lineToParse);
		this.ignoreACCFormat = ignoreACCFormat;
		final String[] elements = lineToParse.split("\t");
		locus = elements[positions.get(LOCUS)];

		spectrumCount = Integer.parseInt(elements[positions.get(SP_COUNT)]);
		setCoverage(Float.parseFloat(elements[positions.get(COVERAGE)].replace("%", "")));

		setLength(Integer.parseInt(elements[positions.get(LENGTH)]));
		setMw(Float.parseFloat(elements[positions.get(MW)]));
		setPi(Float.parseFloat(elements[positions.get(PI)]));
		final String description = elements[positions.get(DESCRIPTION)];
		getPrimaryAccession().setDescription(description);
		// get gene name from description
		final String gene = FastaParser.getGeneFromFastaHeader(description);
		if (gene != null) {
			addGene(new GeneEx(gene));
		}
		setNsaf_norm(Float.valueOf(spectrumCount) / Float.valueOf(length));

		if (positions.containsKey(NSAF)) {
			setNsaf(Float.valueOf(elements[positions.get(NSAF)]));
		} else {
			setNsaf(null);
		}
		if (positions.containsKey(EMPAI)) {
			setEmpai(Float.parseFloat(elements[positions.get(EMPAI)]));
		} else {
			setEmpai(null);
		}
	}

	public DTASelectProtein(IndexedProtein indexedProtein, boolean ignoreACCFormat) {
		this.ignoreACCFormat = ignoreACCFormat;
		locus = indexedProtein.getAccession();
		final Accession accession2 = FastaParser.getACC(locus);
		super.setPrimaryAccession(accession2);
		final String description = indexedProtein.getFastaDefLine();
		accession2.setDescription(description);
		spectrumCount = null;
		setCoverage(null);
		length = null;
		setMw(null);
		setPi(null);
		setNsaf_norm(null);
		setNsaf(null);
		setEmpai(null);
	}

	@Override
	public Integer getLength() {
		return length;
	}

	@Override
	public Integer getSpectrumCount() {
		return spectrumCount;
	}

	/**
	 * Gets the empai number calculated as ((Math.pow(10, coverage / 100)) - 1)
	 *
	 * @return
	 */
	public double getEmpaiCov() {
		return ((Math.pow(10, getCoverage() / 100)) - 1);
	}

	public String getLocus() {
		return locus;
	}

	public void setLocus(String locus) {
		this.locus = locus;
		// set acc to null, since it is a value derived from id
		final Accession acc = null;
		this.setPrimaryAccession(acc);
	}

	/**
	 *
	 * @return
	 */
	public List<Protein> getSibilingProteinsInGroup() {
		final List<Protein> ret = new ArrayList<Protein>();

		for (final GroupableProtein protein : getProteinGroup()) {
			if (protein != this && !protein.getAccession().equals(getAccession()))
				ret.add((Protein) protein);
		}

		// if (!ret.isEmpty())
		// log.info("Protein with " + ret.size() + " other proteins in "
		// + proteinGroups.size() + " groups");
		return ret;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return getLocus();
	}

	@Override
	public void mergeWithProtein(Protein protein) {

		super.mergeWithProtein(protein);
		if (protein instanceof DTASelectProtein) {
			if (((DTASelectProtein) protein).getLocus() != null) {
				setLocus(((DTASelectProtein) protein).getLocus());
			}
		}

	}

	public Set<String> getPeptideSequences() {
		final Set<String> peptideSequences = new THashSet<String>();
		final List<PSM> psMs2 = getPSMs();
		for (final PSM psm : psMs2) {
			final DTASelectPSM dtaSelectPSM = (DTASelectPSM) psm;
			peptideSequences.add(dtaSelectPSM.getSequence());
		}
		return peptideSequences;
	}

	/**
	 * returns the result of calling to FastaParser.getACC(getLocus))
	 * 
	 * @return
	 */
	@Override
	public Accession getPrimaryAccession() {
		if (super.getPrimaryAccession() == null) {
			if (ignoreACCFormat) {
				setPrimaryAccession(AccessionType.UNIPROT, getLocus());
			} else {
				setPrimaryAccession(FastaParser.getACC(getLocus()));
			}
		}
		return super.getPrimaryAccession();
	}

	@Override
	public Set<MSRun> getMSRuns() {
		if (super.getMSRuns() == null || super.getMSRuns().isEmpty()) {
			// get msruns from psms
			if (getPSMs() != null) {
				for (final PSM psm : getPSMs()) {
					addMSRun(psm.getMSRun());
				}
			}
			// getPSMs().stream().forEach(psm -> addMSRun(psm.getMSRun()));
		}
		return super.getMSRuns();
	}

	@Override
	public boolean addPSM(PSM psm, boolean recursively) {
		return super.addPSM(psm, recursively);
	}

}
