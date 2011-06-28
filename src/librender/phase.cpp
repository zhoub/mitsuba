/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>

MTS_NAMESPACE_BEGIN

std::string PhaseFunctionQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "PhaseFunctionQueryRecord[" << std::endl
		<< "  mRec = " << indent(mRec.toString()) << "," << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  wo = " << wo.toString() << "," << std::endl
		<< "  quantity = " << quantity << std::endl
		<< "]";
	return oss.str();
}

Float PhaseFunction::pdf(const PhaseFunctionQueryRecord &pRec) const {
	return f(pRec);
}
	
bool PhaseFunction::needsDirectionallyVaryingCoefficients() const {
	return false;
}
	
Float PhaseFunction::sigmaDir(Float cosTheta) const {
	Log(EError, "sigmaDir(): Not implemented! (this is not"
		" an anisotropic medium)");
	return 0.0f;
}
	
Float PhaseFunction::sigmaDirMax() const {
	Log(EError, "sigmaDirMax(): Not implemented! (this is not"
		" an anisotropic medium)");
	return 0.0f;
}

MTS_IMPLEMENT_CLASS(PhaseFunction, true, ConfigurableObject)
MTS_NAMESPACE_END
