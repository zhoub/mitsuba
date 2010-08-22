#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/**
 * Gamma-corrected bitmap texture using the JPG, TGA or PNG file format
 */
class LDRTexture : public Texture {
public:
	LDRTexture(const Properties &props) : Texture(props) {
		m_filename = props.getString("filename");
		m_filename = FileResolver::getInstance()->resolve(m_filename);
		m_gamma = props.getFloat("gamma", -1); /* -1 means sRGB */
		Log(EInfo, "Loading texture \"%s\"", m_filename.c_str());

		ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
		std::string lower = toLowerCase(m_filename);
	
		if (endsWith(lower, ".jpg") || endsWith(lower, ".jpeg"))
			m_format = Bitmap::EJPEG;
		else if (endsWith(lower, ".png"))
			m_format = Bitmap::EPNG;
		else if (endsWith(lower, ".tga"))
			m_format = Bitmap::ETGA;
		else
			Log(EError, "Cannot deduce the file type of '%s'!", m_filename.c_str());

		ref<Bitmap> bitmap = new Bitmap(m_format, fs);
		initializeFrom(bitmap);
	}

	LDRTexture(Stream *stream, InstanceManager *manager) 
	 : Texture(stream, manager) {
		m_filename = stream->readString();
		Log(EInfo, "Unserializing texture \"%s\"", m_filename.c_str());
		m_gamma = stream->readFloat();
		m_format = static_cast<Bitmap::EFileFormat>(stream->readInt());
		unsigned int size = stream->readUInt();
		ref<MemoryStream> mStream = new MemoryStream(size);
		stream->copyTo(mStream, size);
		mStream->setPos(0);
		ref<Bitmap> bitmap = new Bitmap(m_format, mStream);
		initializeFrom(bitmap);

		if (Scheduler::getInstance()->hasRemoteWorkers()
			&& !FileStream::exists(m_filename)) {
			/* This code is running on a machine different from
			   the one that created the stream. Because we might
			   later have to handle a call to serialize(), the
			   whole bitmap must be kept in memory */
			m_stream = mStream;
		}
	}

	inline Float fromSRGBComponent(Float value) {
		if (value <= (Float) 0.04045)
			return value / (Float) 12.92;
		return std::pow((value + (Float) 0.055)
			/ (Float) (1.0 + 0.055), (Float) 2.4);
	}

	void initializeFrom(Bitmap *bitmap) {
		ref<Bitmap> corrected = new Bitmap(bitmap->getWidth(), bitmap->getHeight(), 128);

		unsigned char *data = bitmap->getData();
		float *flData = corrected->getFloatData();
		Spectrum spec;
		if (bitmap->getBitsPerPixel() == 32) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					if (m_gamma == -1) {
						Float r = (*data++)/255.0f,
							  g = (*data++)/255.0f,
							  b = (*data++)/255.0f,
							  a = (*data++)/255.0f;
						spec.fromSRGB(r, g, b);
						spec.toLinearRGB(r, g, b);

						*flData++ = r;
						*flData++ = g;
						*flData++ = b;
						*flData++ = a;
					} else {
						Float r = std::pow((Float) (*data++)/255.0f, m_gamma),
							  g = std::pow((Float) (*data++)/255.0f, m_gamma),
							  b = std::pow((Float) (*data++)/255.0f, m_gamma),
							  a = (*data++)/255.0f;

						*flData++ = r;
						*flData++ = g;
						*flData++ = b;
						*flData++ = a;
					}
				}
			}
		} else if (bitmap->getBitsPerPixel() == 24) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					if (m_gamma == -1) {
						Float r = (*data++)/255.0f,
							  g = (*data++)/255.0f,
							  b = (*data++)/255.0f;
						spec.fromSRGB(r, g, b);
						spec.toLinearRGB(r, g, b);

						*flData++ = r;
						*flData++ = g;
						*flData++ = b;
						*flData++ = 1.0f;
					} else {
						Float r = std::pow((Float) (*data++)/255.0f, m_gamma),
							  g = std::pow((Float) (*data++)/255.0f, m_gamma),
							  b = std::pow((Float) (*data++)/255.0f, m_gamma);

						*flData++ = r;
						*flData++ = g;
						*flData++ = b;
						*flData++ = 1.0f;
					}
				}
			}
		} else if (bitmap->getBitsPerPixel() == 16) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					if (m_gamma == -1) {
						Float col = (*data++)/255.0f,
							  a = (*data++)/255.0f;
						col = fromSRGBComponent(col);

						*flData++ = col;
						*flData++ = col;
						*flData++ = col;
						*flData++ = a;
					} else {
						Float col = std::pow((Float) (*data++)/255.0f, m_gamma),
							  a = (*data++)/255.0f;

						*flData++ = col;
						*flData++ = col;
						*flData++ = col;
						*flData++ = a;
					}
				}
			}
		} else if (bitmap->getBitsPerPixel() == 8) {
			for (int y=0; y<bitmap->getHeight(); ++y) {
				for (int x=0; x<bitmap->getWidth(); ++x) {
					if (m_gamma == -1) {
						Float col = (*data++)/255.0f;
						col = fromSRGBComponent(col);

						*flData++ = col;
						*flData++ = col;
						*flData++ = col;
						*flData++ = 1.0f;
					} else {
						Float col = std::pow((Float) (*data++)/255.0f, m_gamma);
						*flData++ = col;
						*flData++ = col;
						*flData++ = col;
						*flData++ = 1.0f;
					}
				}
			}
		} else {
			Log(EError, "%i bpp JPG/PNGs are currently not supported!", bitmap->getBitsPerPixel());
		}
	
		m_mipmap = MIPMap::fromBitmap(corrected);
		m_average = m_mipmap->triangle(m_mipmap->getLevels()-1, 0, 0);
		m_maximum = m_mipmap->getMaximum();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture::serialize(stream, manager);
		stream->writeString(m_filename);
		stream->writeFloat(m_gamma);
		stream->writeInt(m_format);
		if (m_stream.get()) {
			stream->writeUInt((unsigned int) m_stream->getSize());
			stream->write(m_stream->getData(), m_stream->getSize());
		} else {
			ref<Stream> mStream = new MemoryStream();
			ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
			stream->writeUInt((unsigned int) is->getSize());
			is->copyTo(stream);
		}
	}

	Spectrum getValue(const Intersection &its) const {
		return m_mipmap->getValue(its);
	}
	
	Spectrum getAverage() const {
		return m_average;
	}

	Spectrum getMaximum() const {
		return m_maximum;
	}

	bool usesRayDifferentials() const {
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "LDRTexture[filename=\"" << m_filename << "\", gamma=" << m_gamma << "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap> m_mipmap;
	ref<MemoryStream> m_stream;
	std::string m_filename;
	Bitmap::EFileFormat m_format;
	Spectrum m_average, m_maximum;
	Float m_gamma;
};

// ================ Hardware shader implementation ================ 

class LDRTextureShader : public Shader {
public:
	LDRTextureShader(Renderer *renderer, std::string filename, ref<Bitmap> bitmap) 
		: Shader(renderer, ETextureShader) {
		m_gpuTexture = renderer->createGPUTexture(filename, bitmap);
		m_gpuTexture->setWrapType(GPUTexture::ERepeat);
		m_gpuTexture->setMaxAnisotropy(8);
		m_gpuTexture->init();
		/* Release the memory on the host side */
		m_gpuTexture->setBitmap(0, NULL);
	}

	void cleanup(Renderer *renderer) {
		m_gpuTexture->cleanup();
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform sampler2D " << evalName << "_texture;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return texture2D(" << evalName << "_texture, uv).rgb;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_texture", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, 
		int &textureUnitOffset) const {
		m_gpuTexture->bind(textureUnitOffset++);
		program->setParameter(parameterIDs[0], m_gpuTexture.get());
	}

	void unbind() const {
		m_gpuTexture->unbind();
	}

	MTS_DECLARE_CLASS()
private:
	ref<GPUTexture> m_gpuTexture;
};

Shader *LDRTexture::createShader(Renderer *renderer) const {
	return new LDRTextureShader(renderer, m_filename, m_mipmap->getLDRBitmap());
}

MTS_IMPLEMENT_CLASS_S(LDRTexture, false, Texture)
MTS_IMPLEMENT_CLASS(LDRTextureShader, false, Shader)
MTS_EXPORT_PLUGIN(LDRTexture, "LDR texture (JPG/TGA/PNG)");
MTS_NAMESPACE_END
