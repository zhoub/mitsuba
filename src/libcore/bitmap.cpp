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

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

#if defined(WIN32)
#undef _CRT_SECURE_NO_WARNINGS
#define _MATH_DEFINES_DEFINED
#endif

#include <ImfRgba.h>
#include <ImfRgbaFile.h>
#include <ImfIO.h>
#include <ImathBox.h>

#include <png.h>
extern "C" {
	#include <jpeglib.h>
	#include <jerror.h>
};

MTS_NAMESPACE_BEGIN

/* ========================== *
 *     EXR helper classes     *
 * ========================== */

class EXRIStream : public Imf::IStream {
public:
	EXRIStream(Stream *stream) : IStream(stream->toString().c_str()),
		m_stream(stream) {
		m_offset = stream->getPos();
	}

	bool read(char *c, int n) {
		m_stream->read(c, n);
		return m_stream->isEOF();
	}

	Imf::Int64 tellg() {
		return m_stream->getPos()-m_offset;
	}

	void seekg(Imf::Int64 pos) {
		m_stream->setPos((size_t) pos + m_offset);
	}

	void clear() {
	}
private:
	ref<Stream> m_stream;
	size_t m_offset;
};

class EXROStream : public Imf::OStream {
public:
	EXROStream(Stream *stream) : OStream(stream->toString().c_str()),
		m_stream(stream) {
	}

	void write(const char *c, int n) {
		m_stream->write(c, n);
	}

	Imf::Int64 tellp() {
		return m_stream->getPos();
	}

	void seekp(Imf::Int64 pos) {
		m_stream->setPos((size_t) pos);
	}

	void clear() {
	}
private:
	ref<Stream> m_stream;
};

/* ========================== *
 *    PNG helper functions    *
 * ========================== */

static void png_flush_data(png_structp png_ptr) {
	voidp flush_io_ptr = png_get_io_ptr(png_ptr);
	((Stream *) flush_io_ptr)->flush();
}

static void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length) {
	voidp read_io_ptr = png_get_io_ptr(png_ptr);
	((Stream *) read_io_ptr)->read(data, length);
}

static void png_write_data(png_structp png_ptr, png_bytep data, png_size_t length) {
	voidp write_io_ptr = png_get_io_ptr(png_ptr);
	((Stream *) write_io_ptr)->write(data, length);
}

static void png_error_func(png_structp png_ptr, png_const_charp msg) {
	SLog(EError, "Fatal libpng error: %s\n", msg);
	exit(-1);
}

/* ========================== *
 *   JPEG helper functions    *
 * ========================== */

extern "C" {
	typedef struct {
		struct jpeg_source_mgr pub;
		JOCTET * buffer;
		size_t buflen;
	} jbuf_t;

	METHODDEF(void) dsm_init_source(j_decompress_ptr cinfo) {
	}


	METHODDEF(boolean) dsm_fill_input_buffer (j_decompress_ptr cinfo) {
		ERREXIT(cinfo, JERR_INPUT_EOF);
		return TRUE;
	}
  
	METHODDEF(void) dsm_skip_input_data (j_decompress_ptr cinfo, long num_bytes) {
		jbuf_t *p = (jbuf_t *)cinfo->src;
		long i = (long) (p->pub.bytes_in_buffer - num_bytes);
		if (i < 0) i = 0;
		p->pub.bytes_in_buffer = i;
		p->pub.next_input_byte += num_bytes;
	}

	METHODDEF(void) dsm_term_source (j_decompress_ptr cinfo) {
	}

	METHODDEF(void) dsm_error_exit (j_common_ptr cinfo) {
		char msg[JMSG_LENGTH_MAX];
		(*cinfo->err->format_message) (cinfo, msg);
		SLog(EError, "Critcal libjpeg error: %s", msg);
	}
};

/* ========================== *
 *        Bitmap class        *
 * ========================== */

Bitmap::Bitmap(int width, int height, int bpp)
 : m_width(width), m_height(height), m_bpp(bpp), m_data(NULL) {
	AssertEx(m_bpp == 1 || m_bpp == 8 || m_bpp == 16 || m_bpp == 24 || m_bpp == 32
		|| m_bpp == 96 || m_bpp == 128, "Invalid number of bits per pixel");
	AssertEx(width > 0 && height > 0, "Invalid bitmap size");

	if (bpp == 96 || bpp == 128)
		m_gamma = 1.0f;
	else
		m_gamma = -1.0f; // sRGB

	// 1-bit masks are stored in a packed format. 
	m_size = (size_t) std::ceil(((double) m_width * m_height * m_bpp) / 8.0);
	m_data = static_cast<unsigned char *>(allocAligned(m_size));
}

Bitmap::Bitmap(EFileFormat format, Stream *stream) : m_data(NULL) {
	if (format == EPNG)
		loadPNG(stream);
	else if (format == EJPEG)
		loadJPEG(stream);
	else if (format == EEXR)
		loadEXR(stream);
	else if (format == ETGA)
		loadTGA(stream);
	else if (format == EBMP)
		loadBMP(stream);
	else
		Log(EError, "Bitmap: Invalid file format!");
}
	
void Bitmap::loadEXR(Stream *stream) {
	EXRIStream istr(stream);
	Imf::RgbaInputFile file(istr);

	/* Determine dimensions and allocate space */
	Imath::Box2i dw = file.dataWindow();
	m_width = dw.max.x - dw.min.x + 1;
	m_height = dw.max.y - dw.min.y + 1;
	m_size = m_width * m_height * 16;
	m_bpp = 4*4*8;
	m_gamma = 1.0f;
	m_data = static_cast<unsigned char *>(allocAligned(m_size));
	Imf::Rgba *rgba = new Imf::Rgba[m_width*m_height];
	Log(ETrace, "Reading %ix%i EXR file", m_width, m_height);

	/* Convert to 32-bit floating point per channel */
	file.setFrameBuffer(rgba, 1, m_width);
	file.readPixels(dw.min.y, dw.max.y);
	float *m_buffer = getFloatData();
	for (int i=0; i<m_width*m_height; i++) {
		*m_buffer = (float) rgba[i].r; m_buffer++;
		*m_buffer = (float) rgba[i].g; m_buffer++;
		*m_buffer = (float) rgba[i].b; m_buffer++;
		*m_buffer = (float) rgba[i].a; m_buffer++;
	}
	delete[] rgba;
}

void Bitmap::loadTGA(Stream *stream) {
	int headerSize = stream->readUChar();
	if (stream->readUChar() != 0)
		Log(EError, "Invalid TGA format -- only raw (non-RLE encoded) RGB is supported for now");
	if (stream->readUChar() != 2)
		Log(EError, "Invalid TGA format -- only raw (non-RLE encoded) RGB is supported for now");
	stream->setPos(8);
	int x1 = stream->readShort();
	int y1 = stream->readShort();
	int x2 = stream->readShort();
	int y2 = stream->readShort();
	m_width = x2-x1;
	m_height = y2-y1;
	Log(EInfo, "Reading %ix%i TGA file", m_width, m_height);

	stream->setPos(16);
	m_bpp = stream->readUChar();
	if (m_bpp != 24 && m_bpp != 32)
		Log(EError, "Invalid TGA format -- only 24 or 32 bpp images are supported for now");

	m_gamma = -1;
	int channels = m_bpp / 8;
	m_size = m_width * m_height * channels;
	m_data = static_cast<unsigned char *>(allocAligned(m_size));
	stream->setPos(18 + headerSize);
	stream->read(m_data, m_size);

	/* Convert BGR to RGB */
	for (size_t i=0; i<m_size; i += channels) {
		uint8_t tmp = m_data[i];
		m_data[i] = m_data[i+2];
		m_data[i+2] = tmp;
	}
}

void Bitmap::loadBMP(Stream *stream) {
#if defined(WIN32)
#pragma pack(push, 1)
#endif
	struct 
#if !defined(WIN32)
		__attribute__((__packed__))
#endif
	{
		uint16_t magic;
		uint32_t size;
		uint16_t creator1, creator2;
		uint32_t offset;
	} BMPFileHeader;

	struct 
#if !defined(WIN32)
	__attribute__((__packed__))
#endif
	{
		uint32_t header_sz;
		uint32_t width;
		int32_t height;
		uint16_t nplanes;
		uint16_t bitspp;
		uint32_t compress_type;
		uint32_t bmp_bytesz;
		uint32_t hres;
		uint32_t vres;
		uint32_t ncolors;
		uint32_t nimpcolors;
	} DIBHeader;
#if defined(WIN32)
#pragma pack(pop)
#endif

	Assert(sizeof(BMPFileHeader) == 14);
	Assert(sizeof(DIBHeader) == 40);

	stream->read(&BMPFileHeader, sizeof(BMPFileHeader));

	if (memcmp(&BMPFileHeader.magic, "BM", 2) != 0)
		Log(EError, "Unsupported BMP format encountered (invalid file header)!");

	stream->read(&DIBHeader, sizeof(DIBHeader));

	if (DIBHeader.header_sz != 40 || DIBHeader.nplanes != 1)
		Log(EError, "Unsupported BMP format encountered (strange DIB header)!");

	if (DIBHeader.ncolors != 0)
		Log(EError, "Only BMP images without a palette are supported for now");

	if (DIBHeader.bitspp != 8 && DIBHeader.bitspp != 24)
		Log(EError, "Only 8- and 24-bit BMP images are supported for now");

	if (DIBHeader.compress_type != 0)
		Log(EError, "Only uncompressed BMP images are supported for now");

	m_width = DIBHeader.width;
	m_height = std::abs(DIBHeader.height);
	m_bpp = DIBHeader.bitspp;

	m_size = m_width * m_height * (m_bpp / 8);
	m_data = static_cast<unsigned char *>(allocAligned(m_size));
	Log(ETrace, "Reading %ix%ix%i BMP file", m_width, m_height, m_bpp);

	int nChannels = m_bpp / 8;
	int lineWidth = m_width * nChannels; 
	int paddedLineWidth = lineWidth + lineWidth % 4;

	uint8_t *line = new uint8_t[paddedLineWidth];

	for (int y=0; y<m_height; ++y) {
		stream->read(line, paddedLineWidth);

		int targetY = y;

		if (DIBHeader.height > 0)
			targetY = m_height - 1 - y; // inverted rows

		memcpy(&m_data[targetY * m_width * nChannels], line, lineWidth);
	
		if (nChannels == 3) {
			for (int x=0; x<m_width; ++x)
				std::swap(m_data[(targetY * m_width + x) * nChannels], 
						  m_data[(targetY * m_width + x) * nChannels + 2]);
		}
	}

	delete[] line;
}

void Bitmap::loadPNG(Stream *stream) {
	png_structp png_ptr;
	png_infop info_ptr;
	volatile png_bytepp rows = NULL;

	/* Create buffers */
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, &png_error_func, NULL);
	if (png_ptr == NULL) {
		Log(EError, "Error while creating PNG data structure");
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		png_destroy_read_struct(&png_ptr, (png_infopp) NULL, (png_infopp) NULL);
		Log(EError, "Error while creating PNG information structure");
	}

	/* Error handling */
	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);
		if (rows)
			delete[] rows;
		Log(EError, "Error reading the PNG file");
	}

	/* Set read helper function */
	png_set_read_fn(png_ptr, stream, (png_rw_ptr) png_read_data);

	int bitdepth, colortype, interlacetype, compressiontype, filtertype;
	png_read_info(png_ptr, info_ptr);
	png_uint_32 width=0, height=0;
	png_get_IHDR(png_ptr, info_ptr, &width, &height, &bitdepth,
			&colortype, &interlacetype, &compressiontype, &filtertype);
	int newDepth = bitdepth;

	if (bitdepth == 1) {
		png_set_packing(png_ptr); // Unpack and later re-pack
	} else if (colortype == PNG_COLOR_TYPE_PALETTE) {
		png_set_expand(png_ptr); // expand indexed files
		newDepth = 8;
	} else if (colortype == PNG_COLOR_TYPE_GRAY && bitdepth < 8) {

		png_set_expand_gray_1_2_4_to_8(png_ptr); // convert grayscale to 8bit
		newDepth = 8;
	} else if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) {
		png_set_expand(png_ptr); // transparency
	} else if (bitdepth < 8) {
		newDepth = 8;
		png_set_expand(png_ptr);
	}

	// handle interlacing
	if (interlacetype != PNG_INTERLACE_NONE)
		png_set_interlace_handling(png_ptr);

	if (bitdepth == 1) {
		m_bpp = 1; 
	} else if (colortype == PNG_COLOR_TYPE_GRAY) {
		m_bpp = newDepth;
	} else if (colortype == PNG_COLOR_TYPE_GRAY_ALPHA) {
		m_bpp = newDepth*2;
	} else {
		m_bpp = newDepth * ((colortype & PNG_COLOR_MASK_ALPHA) ? 4 : 3);
	}

	int intent; double gamma;
	if (png_get_sRGB(png_ptr, info_ptr, &intent)) {
		m_gamma = -1;
	} else if (png_get_gAMA(png_ptr, info_ptr, &gamma)) {
		m_gamma = (Float) gamma;
	} else {
		m_gamma = 1.0f/2.2f;
	}

	Log(ETrace, "Reading %ix%ix%i PNG file", width, height, m_bpp);

	/* Update the information */
	png_read_update_info(png_ptr, info_ptr);

	/* re-read */
	png_get_IHDR(png_ptr, info_ptr, &width, &height, &bitdepth,
			&colortype, &interlacetype, &compressiontype, &filtertype);
	m_width = width; m_height = height;

	m_size = (size_t) std::ceil((m_width * m_height * m_bpp) / 8.0);
	m_data = static_cast<unsigned char *>(allocAligned(m_size));

	rows = new png_bytep[m_height];

	if (m_bpp == 1) {
		for (int i=0; i<m_height; i++)
			rows[i] = new unsigned char[m_width];
	} else {
		int rowBytes = m_width * (m_bpp / 8);
		for (int i=0; i<m_height; i++)
			rows[i] = m_data + i * rowBytes;
	}

	png_read_image(png_ptr, rows);
	png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);

	if (m_bpp == 1) {
		// The bitmask has been unpacked by the decoding, now re-pack it
		memset(m_data, 0, m_size);
		for (int i=0; i<m_height; i++) {
			for (int o=0; o<m_width; o++) {
				int pos = i * m_width + o;
				m_data[pos / 8] |= rows[i][o] * (1 << (pos % 8));
			}
			delete[] rows[i];
		}
	}

	delete[] rows;
}

void Bitmap::loadJPEG(Stream *stream) {
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	jbuf_t jbuf;
	size_t length = stream->getSize();

	memset(&jbuf, 0, sizeof(jbuf_t));
	cinfo.err = jpeg_std_error(&jerr);
	jerr.error_exit = dsm_error_exit;

	jpeg_create_decompress(&cinfo);
	cinfo.src = (struct jpeg_source_mgr *) &jbuf;
	jbuf.buffer = new JOCTET[length];
	jbuf.pub.init_source = dsm_init_source;
	jbuf.pub.fill_input_buffer = dsm_fill_input_buffer;
	jbuf.pub.skip_input_data = dsm_skip_input_data;

 	/* Use default method (in libjpeg) */
	jbuf.pub.resync_to_restart = jpeg_resync_to_restart;
	jbuf.pub.term_source = dsm_term_source;
	jbuf.buflen = jbuf.pub.bytes_in_buffer = length;
	jbuf.pub.next_input_byte = jbuf.buffer;
  
	stream->read(jbuf.buffer, length);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	m_width = cinfo.output_width;
	m_height = cinfo.output_height;
	m_bpp = cinfo.output_components*8;
	m_gamma = 1.0f/2.2f;
	m_size = m_width * m_height * cinfo.output_components;
	Log(ETrace, "Reading %ix%ix%i JPG file", m_width, m_height, m_bpp);

	int row_stride = cinfo.output_width * cinfo.output_components;
	unsigned char **buffer = new unsigned char *[m_height];
	m_data = static_cast<unsigned char *>(allocAligned(m_size));
	for (int i=0; i<m_height; ++i) 
		buffer[i] = m_data + row_stride*i;

	/* Process scanline by scanline */	
	int counter = 0;
	while (cinfo.output_scanline < (unsigned int) m_height) 
		counter += jpeg_read_scanlines(&cinfo, &buffer[counter], m_height);

	/* Release the libjpeg data structures */
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);

	delete[] buffer;
	delete[] jbuf.buffer;
}

Bitmap::~Bitmap() {
	if (m_data)
		freeAligned(m_data);
}

void Bitmap::clear() {
	memset(m_data, 0, m_size);
}

Bitmap *Bitmap::clone() const {
	Bitmap *bitmap = new Bitmap(m_width, m_height, m_bpp);
	memcpy(bitmap->m_data, m_data, m_size);
	return bitmap;
}

bool Bitmap::operator==(const Bitmap &bitmap) const {
	if (bitmap.m_width == m_width &&
		bitmap.m_height == m_height &&
		bitmap.m_bpp == m_bpp) {
		return memcmp(bitmap.m_data, m_data, m_size) == 0;
	}
	return false;
}

void Bitmap::save(EFileFormat format, Stream *stream, int compression) const {
	if (m_bpp == 96 || m_bpp == 128)
		AssertEx(format == EEXR, "Bitmap: 96/128 bpp images can only be stored "
			"using the EXR file format");
	else
		AssertEx(format == EPNG, "Bitmap: 1-32 bpp images can only be stored "
			"using the PNG file format");

	if (format == EEXR)
		saveEXR(stream);
	else if (format == EPNG)
		savePNG(stream, compression);
	else
		Log(EError, "Bitmap: Invalid file format!");
}

void Bitmap::saveEXR(Stream *stream) const {
	Log(EDebug, "Writing %ix%i EXR file", m_width, m_height);
	EXROStream ostr(stream);
	Imf::RgbaOutputFile file(ostr, Imf::Header(m_width, m_height), 
		Imf::WRITE_RGBA);

	Imf::Rgba *rgba = new Imf::Rgba[m_width*m_height];
	const float *m_buffer = getFloatData();
	for (int i=0; i<m_width*m_height; i++) {
		rgba[i].r = *m_buffer; m_buffer++;
		rgba[i].g = *m_buffer; m_buffer++;
		rgba[i].b = *m_buffer; m_buffer++;
		if (m_bpp == 128) {
			rgba[i].a = *m_buffer; m_buffer++;
		} else {
			rgba[i].a = 1;
		}
	}
	file.setFrameBuffer(rgba, 1, m_width);
	file.writePixels(m_height);
	delete[] rgba;
}

void Bitmap::savePNG(Stream *stream, int compression) const {
	png_structp png_ptr;
	png_infop info_ptr;
	png_text text[4];
	volatile png_bytepp rows = NULL;

	if (m_gamma == -1)
		Log(EDebug, "Writing %ix%ix%i sRGB PNG file", m_width, m_height, m_bpp);
	else
		Log(EDebug, "Writing %ix%ix%i PNG file (gamma=%.1f)", m_width, m_height, m_bpp, m_gamma);

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, &png_error_func, NULL);
	if (png_ptr == NULL) {
		Log(EError, "Error while creating PNG data structure");
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
		Log(EError, "Error while creating PNG information structure");
	}

	/* Error handling */
	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		Log(EError, "Error writing the PNG file");
	}

	png_set_write_fn(png_ptr, stream, (png_rw_ptr) png_write_data, (png_flush_ptr) png_flush_data);
	png_set_compression_level(png_ptr, compression);

	memset(text, 0, sizeof(png_text)*4);
	text[0].key = (char *) "Generated by";
	text[0].text =  (char *) "Vitsuba version " MTS_VERSION;
	text[0].compression = PNG_TEXT_COMPRESSION_NONE;
	text[1].key = (char *) "Title";
	text[1].text = (char *) m_title.c_str();
	text[1].compression = PNG_TEXT_COMPRESSION_NONE;
	text[2].key = (char *) "Author";
	text[2].text = (char *) m_author.c_str();
	text[2].compression = PNG_TEXT_COMPRESSION_NONE;
	text[3].key = (char *) "Comment";
	text[3].text = (char *) m_comment.c_str();
	text[3].compression = PNG_TEXT_COMPRESSION_zTXt;
	png_set_text(png_ptr, info_ptr, text, 4);

	if (m_gamma == -1)
		png_set_sRGB_gAMA_and_cHRM(png_ptr, info_ptr, PNG_sRGB_INTENT_ABSOLUTE);
	else
		png_set_gAMA(png_ptr, info_ptr, m_gamma);

	int colortype;
	switch (m_bpp) {
		case 32: colortype = PNG_COLOR_TYPE_RGBA; break;
		case 24: colortype = PNG_COLOR_TYPE_RGB; break;
		case 16: colortype = PNG_COLOR_TYPE_GRAY_ALPHA; break; 
		case 8: colortype = PNG_COLOR_TYPE_GRAY; break; 
		default: colortype = PNG_COLOR_TYPE_PALETTE;
	}

	/* Simple 8 bit/color RGB/RGBA format or 1-bit mask */
	png_set_IHDR(png_ptr, info_ptr, m_width, m_height, m_bpp == 1 ? 1 : 8,
			colortype, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
			PNG_FILTER_TYPE_BASE);

	if (m_bpp == 1) {
		png_color palette[2];
		palette[0].blue = palette[0].red = palette[0].green = 0;
		palette[1].blue = palette[1].red = palette[1].green = 0xFF;
		png_set_PLTE(png_ptr, info_ptr, palette, 2);
	}
	
	png_write_info(png_ptr, info_ptr);
	png_set_packing(png_ptr);

	rows = new png_bytep[m_height];

	if (m_bpp == 1) {
		/* Convert to 8 bit */
		for (int i=0; i<m_height; i++) {
			rows[i] = new unsigned char[m_width];
			for (int o=0; o<m_width; o++) {
				int pos = i * m_width + o;
				rows[i][o] = (m_data[pos / 8] & (1 << (pos % 8))) == false ? 0 : 1;
			}
		}
	} else {
		int rowBytes = m_width * (m_bpp / 8);
		for (int i=0; i<m_height; i++) 
			rows[i] = &m_data[rowBytes * i];
	}

	png_write_image(png_ptr, rows);
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);

	if (m_bpp == 1) {
		for (int i=0; i<m_height; i++)
			delete[] rows[i];
	}

	delete[] rows;
}

std::string Bitmap::toString() const {
	std::ostringstream oss;
	oss << "Bitmap[width=" << m_width << ", height=" << m_height << ", bpp=" << m_bpp << "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(Bitmap, false, Object)
MTS_NAMESPACE_END
