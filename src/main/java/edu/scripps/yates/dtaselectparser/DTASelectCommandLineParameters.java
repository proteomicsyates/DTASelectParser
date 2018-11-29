package edu.scripps.yates.dtaselectparser;

import java.util.Map;

import edu.scripps.yates.utilities.parsers.idparser.CommandLineParameters;
import gnu.trove.map.hash.THashMap;

/**
 * Parse command line parameters of DTASElect like -p 1 -y 2 --trypstat --fp
 * 0.01 --modstat --extra --pI -DM 5 --DB --dm -m 0 -S 3.5 --quiet
 * 
 * @author Salva
 * 
 */
public class DTASelectCommandLineParameters implements CommandLineParameters {

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		for (final String key : parameters.keySet()) {
			String value = "";
			if (parameters.get(key) != null) {
				value = parameters.get(key);
			}
			if (!"".equals(sb.toString())) {
				sb.append(" ");
			}
			sb.append(key);
			if (!"".equals(value)) {
				sb.append(" " + value);
			}
		}
		return sb.toString();
	}

	private final Map<String, String> parameters = new THashMap<String, String>();

	public DTASelectCommandLineParameters(String commandLineParameters) {
		final String[] split = commandLineParameters.trim().split(" ");

		String parameter = null;
		String value = null;
		for (final String string : split) {
			if (string.startsWith("-")) {
				if (parameter != null) {
					parameters.put(parameter, null);
				}
				parameter = string;
				value = null;
			} else {
				value = string;
				if (parameter != null) {
					parameters.put(parameter, value);
					parameter = null;
				}
			}
		}
		// insert the last one
		if (parameter != null) {
			parameters.put(parameter, value);
			parameter = null;
		}
	}

	/**
	 * @return the parameters
	 */
	@Override
	public Map<String, String> getParametersMap() {
		return parameters;
	}

	public static void main(String[] args) {
		String paramString = " -p 1 -y 2 --trypstat --fp 0.01 --modstat --extra --pI -DM 5 --DB --dm -m 0 -S 3.5 --quiet ";
		DTASelectCommandLineParameters d = new DTASelectCommandLineParameters(paramString);
		System.out.println(paramString);
		Map<String, String> parameters2 = d.getParametersMap();
		for (final String param : parameters2.keySet()) {
			System.out.println(param + "\t" + parameters2.get(param));
		}

		String parameterValue = d.getParameterValue("t");
		System.out.println("PARAMETER t= '" + parameterValue + "'");

		paramString = " -p 1 -y 2 --trypstat --fp 0.01 --modstat --extra --pI -DM 5 --DB --dm -m 0 -S 3.5 --quiet -t 0";
		d = new DTASelectCommandLineParameters(paramString);
		System.out.println(paramString);
		parameters2 = d.getParametersMap();
		for (final String param : parameters2.keySet()) {
			System.out.println(param + "\t" + parameters2.get(param));
		}

		parameterValue = d.getParameterValue("t");
		System.out.println("PARAMETER t= '" + parameterValue + "'");
	}

	/**
	 * Search for the value of that parameter, avoiding if necessary the "-" or
	 * "--" in front of the parameters
	 * 
	 * @param parameterName
	 * @return
	 */
	@Override
	public String getParameterValue(String parameterName) {
		if (parameters.containsKey(parameterName)) {
			return parameters.get(parameterName);
		} else {
			for (final String parameter : parameters.keySet()) {
				final String parameterWithoutDash = parameter.replace("-", "");
				if (parameterWithoutDash.equals(parameterName))
					return parameters.get(parameter);
			}
		}
		return null;
	}
}
