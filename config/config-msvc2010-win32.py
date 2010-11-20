CXX            = 'cl'
GCC            = 'cl'
# /Ox=optimize for speed, global optimizations, intrinsic functions, favor fast code, frame pointer omission
# /EHsc=C++ exceptions, /fp:fast=Enable reasonable FP optimizations, /GS-=No buffer security checks, /GL=whole program optimizations
# To include debug information add '/Z7' to CXXFLAGS and '/DEBUG' to LINKFLAGS
CXXFLAGS       = ['/nologo', '/Ox', '/fp:fast', '/arch:SSE2' ,'/D', 'WIN32', '/W3', '/EHsc', '/GS-', '/GL', '/MD', '/D', 'MTS_DEBUG', '/D', 'SINGLE_PRECISION', '/D', 'MTS_SSE', '/D', 'MTS_HAS_COHERENT_RT', '/D', '_CONSOLE', '/D', 'NDEBUG', '/openmp']
SHCXXFLAGS     = CXXFLAGS
TARGET_ARCH    = 'x86'
MSVC_VERSION   = '10.0'
LINKFLAGS      = ['/nologo', '/SUBSYSTEM:CONSOLE', '/MACHINE:X86', '/FIXED:NO', '/OPT:REF', '/OPT:ICF', '/LTCG', '/NODEFAULTLIB:LIBCMT']
BASEINCLUDE    = ['#include', '#tools/windows/include']
BASELIB        = ['pthreadVCE2', 'msvcrt', 'ws2_32']
OEXRINCLUDE    = ['#tools/windows/include/OpenEXR']
OEXRFLAGS      = ['/D', 'OPENEXR_DLL']
OEXRLIB        = ['IlmImf', 'IlmThread', 'Iex', 'zlib1', 'Half']
BOOSTINCLUDE   = ['#tools/boost']
BOOSTLIB       = ['boost_system-vc100-mt-1_44', 'boost_filesystem-vc100-mt-1_44']
COLLADAINCLUDE = ['#tools/windows/include/colladadom', '#tools/windows/include/colladadom/1.4']
COLLADALIB     = ['libcollada14dom21']
XERCESLIB      = ['xerces-c_3']
PNGLIB         = ['libpng13']
JPEGLIB        = ['jpeg62']
GLLIB          = ['opengl32', 'glu32', 'glew32mx', 'gdi32', 'user32']
GLFLAGS        = ['/D', 'GLEW_MX']
BASELIBDIR     = ['#tools/windows/lib32', '#tools/windows/lib32/vc100']
SHLIBPREFIX    = 'lib'
SHLIBSUFFIX    = '.dll'
PROGSUFFIX     = '.exe'
