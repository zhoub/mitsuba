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

#if !defined(__SPIRAL_H)
#define __SPIRAL_H

#include <mitsuba/render/block.h>
#include <mitsuba/render/film.h>

#define MTS_BLOCK_SIZE 32

MTS_NAMESPACE_BEGIN

/// Block listener callback for use with the Spiral class
class MTS_EXPORT_RENDER BlockListener {
public:
	/// Called whenever an image block is acquired
	virtual void acquireBlockEvent(const ImageBlock *block) = 0;

	/// Called whenever an image block is released
	virtual void releaseBlockEvent(const ImageBlock *block) = 0;

	/// Called when the whole film has changed
	virtual void filmChangedEvent() = 0;

	/// Called when rendering is done
	virtual void finishEvent() = 0;
protected:
	virtual ~BlockListener() {}
};

/**
 * Class wraps code that generates the spiral of image blocks
 * for rendering.
 *
 * \author Adam Arbree
 * Aug 25, 2005
 * RayTracer.java
 * Copyright 2005 Program of Computer Graphics, Cornell University
  */

class MTS_EXPORT_RENDER Spiral : public Object {
public:
	/**
	 * Create a new spiral generator for an image
	 * of the given width and height.
	 */
	Spiral(const Film *film);

	/**
	 * Reset the spiral to its initial state
	 */
	void reset();

	/**
	 * Add a block listener, which will be notified
	 * whenever a block has been acquired or released.
	 */
	void addBlockListener(BlockListener *listener);

	/// Remove a block listener
	void removeBlockListener(BlockListener *listener);

	/// Acquire an image block from the spiral (thread-safe)
	bool acquireBlock(ImageBlock *block);

	/// Release a finished image block (thread-safe)
	void releaseBlock(ImageBlock *block);

	/// Send the finished event even if not all blocks have been processed yet
	void finish();
	
	/// Send events notifying all listeners that the film has changed
	void notifyFilmChanged();

	/// Return the list blocks, which are currently being worked on
	std::vector<ImageBlock *> getActiveBlocks() const;

	/// Return the maximum block size
	int getMaxBlockSize() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Spiral();
private:
	/* variables controling ordering of sub-blocks */
	int numSubX;             //number of subblocks in x direction
	int numSubY;             //number of subblocks in y direction
	int totalSubblocks;    //total number of subblocks to be computed in all passes
	int subblocksLeft;     //remaining number of subblocks to be computed
	int curSubX, curSubY;  //current subblock to request
	int curDir;              //direction to next subblock
	int movesLeft;           //moves left in this direction
	int moveLength;          //total moves before next direction change
	int width, height;
	int m_activeBlocks;
	int negOffsetX, negOffsetY;
	mutable ref<Mutex> m_mutex;
	std::vector<ImageBlock *> m_openBlocks;
	std::vector<BlockListener *> m_listeners;
	ProgressReporter *m_progress;
};

MTS_NAMESPACE_END

#endif
