package edu.scripps.yates.dtaselectparser.util;

import edu.scripps.yates.utilities.fasta.FastaParser;

public class KeyUtils {
	/**
	 *
	 * @param psm
	 * @param chargeSensible
	 *            if true, then, the charge will be considered for
	 *            differentiating peptides with different charge states. If
	 *            false, peptides with different charge states will have the
	 *            same key
	 * @return
	 */
	public static String getSpectrumKey(DTASelectPSM psm, boolean chargeSensible) {

		StringBuilder sb = new StringBuilder();
		if (psm.getRawFileName() != null) {
			sb.append(psm.getRawFileName());
		}
		if (!"".equals(sb.toString())) {
			sb.append("-");
		}
		if (psm.getScan() != null) {
			sb.append(psm.getScan());
		}
		if (!"".equals(sb.toString())) {
			sb.append("-");
		}
		if (psm.getFullSequence() != null) {
			sb.append(FastaParser.getSequenceInBetween(psm.getFullSequence()));
		}

		if (chargeSensible) {
			if (!"".equals(sb.toString())) {
				sb.append("-");
			}
			sb.append(psm.getChargeState());
		}
		return sb.toString();
	}
}
