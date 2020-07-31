package org.ohdsi.metaAnalysis;

import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.operators.*;

import java.util.List;

public interface Analysis {

	Likelihood getJoint();
	OperatorSchedule getSchedule();
	List<Loggable> getLoggerColumns();
}
