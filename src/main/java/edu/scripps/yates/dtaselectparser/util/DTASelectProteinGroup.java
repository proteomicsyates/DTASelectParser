package edu.scripps.yates.dtaselectparser.util;

import java.util.ArrayList;

public class DTASelectProteinGroup extends ArrayList<DTASelectProtein> {
	/**
	 * 
	 */
	private static final long serialVersionUID = -1051891453675105485L;

	/**
	 * This assumes that the {@link DTASelectProtein} only belongs to one group
	 */
	@Override
	public boolean add(DTASelectProtein e) {
		if (super.contains(e)) {
			return false;
		}
		e.addProteinGroup(this);
		return super.add(e);
	}

}
