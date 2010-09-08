/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

Luminaire::Luminaire(const Properties &props)
 : ConfigurableObject(props), m_surfaceArea(0.0f) {
	m_type = 0;
	m_luminaireToWorld = props.getTransform("toWorld", Transform());
	m_worldToLuminaire = m_luminaireToWorld.inverse();
	AssertEx(!m_luminaireToWorld.hasScale(), "toWorld transformation can't have scale factors!");
	m_intersectable = false;
}

Luminaire::Luminaire(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
	m_surfaceArea = stream->readFloat();
	m_type = (EType) stream->readInt();
	m_intersectable = stream->readBool();
	m_worldToLuminaire = Transform(stream);
	m_luminaireToWorld = m_worldToLuminaire.inverse();
}

Luminaire::~Luminaire() {
}

void Luminaire::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);

	stream->writeFloat(m_surfaceArea);
	stream->writeInt(m_type);
	stream->writeBool(m_intersectable);
	m_worldToLuminaire.serialize(stream);
}

void Luminaire::preprocess(const Scene *scene) {
}

bool Luminaire::isBackgroundLuminaire() const {
	return false;
}

Spectrum Luminaire::Le(const Ray &ray) const {
	return Spectrum(0.0f);
}

MTS_IMPLEMENT_CLASS(Luminaire, true, ConfigurableObject)
MTS_NAMESPACE_END
