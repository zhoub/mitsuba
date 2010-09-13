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

#if !defined(__FILM_H)
#define __FILM_H

#include <mitsuba/render/sampler.h>
#include <mitsuba/render/renderproc_wr.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

MTS_NAMESPACE_BEGIN

class Bitmap;

/** \brief Abstract Film base class - used to store samples
 * generated by the Integrator.
 */
class MTS_EXPORT_RENDER Film : public ConfigurableObject {
public:
	/// Add an image block to the film
	virtual void putImageBlock(const ImageBlock *block) = 0;

	/**
	 * Get a pixel value while the rendering is still in progress
	 * - coordinates are specified in (X, Y) sensor pixel positions
	 */
	virtual Spectrum getValue(int x, int y) = 0;

	/// Clear the film
	virtual void clear() = 0;

	/**
	 * Overwrite the film with the given bitmap (of equal size).
	 * It is assumed that the bitmap uses a linear color space
	 */
	virtual void fromBitmap(const Bitmap *bitmap) = 0;

	/**
	 * Inverse of the above -- store the film contents to a bitmap
	 */
	virtual void toBitmap(Bitmap *bitmap) const = 0;

	/// Develop the film and write the result to the specified filename
	virtual void develop(const fs::path &fileName) = 0;

	/// Ignoring the crop window, return the resolution of the underlying sensor
	inline const Vector2i &getSize() const { return m_size; }
	
	/// Return the size of the crop window
	inline const Vector2i &getCropSize() const { return m_cropSize; }

	/// Return the offset of the crop window
	inline const Point2i &getCropOffset() const { return m_cropOffset; }

	/** 
	 * Should regions slightly outside the image plane be sampled to improve 
	 * the quality of the reconstruction at the edges? This only makes
	 * sense when reconstruction filters other than the box filter are used.
	 */
	inline bool hasHighQualityEdges() const { return m_highQualityEdges; }

	/// Return the tabulated image reconstruction filter
	inline const TabulatedFilter *getTabulatedFilter() const { return m_tabulatedFilter.get(); }
	
	/// Return the original image reconstruction filter
	inline const ReconstructionFilter *getReconstructionFilter() const { return m_filter.get(); }

	/// Add a child node
	virtual void addChild(const std::string &name, ConfigurableObject *child);

	/// Configure the film
	virtual void configure();
	
	/// Serialize this film to disk
	virtual void serialize(Stream *stream, InstanceManager *manager) const;
	
	/// Does the destination already exist?
	virtual bool destinationExists(const fs::path &baseName) const = 0;

	/// Return the properties of this film
	inline const Properties &getProperties() const { return m_properties; }

	MTS_DECLARE_CLASS()
protected:
    /// Create a film
    Film(const Properties &props);
    
	/// Unserialize a film
    Film(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~Film();
protected:
	Point2i m_cropOffset;
	Vector2i m_size, m_cropSize;
	bool m_highQualityEdges;
	ref<ReconstructionFilter> m_filter;
	ref<TabulatedFilter> m_tabulatedFilter;
	Properties m_properties;
};

MTS_NAMESPACE_END

#endif /* __FILM_H */
