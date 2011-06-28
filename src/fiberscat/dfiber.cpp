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

#include <mitsuba/render/fiberscat.h>

MTS_NAMESPACE_BEGIN

/**
 * Diffuse fiber scattering function
 */
class DiffuseFiber : public FiberScatteringFunction {
public:
	DiffuseFiber(const Properties &props) 
		: FiberScatteringFunction(props) {
	}

	DiffuseFiber(Stream *stream, InstanceManager *manager) 
		: FiberScatteringFunction(stream, manager) {
	}

	virtual ~DiffuseFiber() { }


	void serialize(Stream *stream, InstanceManager *manager) const {
		FiberScatteringFunction::serialize(stream, manager);
	}

	Spectrum sample(FiberScatteringRecord &pRec, 
			Sampler *sampler) const {
		Log(EError, "To be implemented");
		return Spectrum(0.0f);
	}

	Spectrum f(const FiberScatteringRecord &pRec) const {
		Log(EError, "To be implemented");
		return Spectrum(0.0f);
	}
	
	Float pdf(const FiberScatteringRecord &pRec) const {
		Log(EError, "To be implemented");
		return 0.0f;
	}

	std::string toString() const {
		return "DiffuseFiber[]";
	}

	MTS_DECLARE_CLASS()
};


MTS_IMPLEMENT_CLASS_S(DiffuseFiber, false, FiberScatteringFunction)
MTS_EXPORT_PLUGIN(DiffuseFiber, "Diffuse fiber scattering function");
MTS_NAMESPACE_END
